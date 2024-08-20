# =============================================================================
#
# Hydrogeological Modflow NWT model built in Flopy.
# Simulates the unconfined groundwater in a complex topography in a steady-state.
# It takes groundwater recharge as main input via the MODFLOW RCH package. 
# Simulations are performed in parallel over a wide range of geomorphoclimatic conditions.
# The outflows are simulated using the MODFLOW DRN package. Drainage is distributed 
# on the entire ground surface as to simulate seepage faces (unknown a priori) naturally occurring
# as a result of the flowfield in the subsurface . 
# Input must include the DTM of the catchment , as well as a binary raster (CR) that acts as a mask  for the basin 
# and identifies the cells to which the drains are assigned (in this case, all of them).
# The final module performs the tracking of particles distributed uniformly on the surface where recharge
# is applied. It then calculates the residence time and the flowpath lengths of the particles and the age of 
# the streamflow (i.e. the outflow at the catchment outlet composed by quick runoff and groundwater contributions).
# 
# itmuni in the dis package defines the units of time , while lenuni defines the spatial units.
# 
#
# See "Morphological and hydrogeological controls of groundwater flows and water age distribution in mountain aquifers and streams"
# by A.Betterle and A.Bellin (Water Resources Reserch) for further details. 
# =============================================================================



import flopy
import flopy.utils.modpathfile as mpf
import os
import shutil
import sys
import flopy.utils.binaryfile as bf
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import multiprocessing
from pathlib import Path

###################################################################
#Set input Files and Directories ##################################
###################################################################

## Set the directories of the modflow exe
exeMODFLOW = Path("/**/**/MODFLOW-NWT_64.exe")
exeMODPATH = Path("/**/**/mpath7.exe")


#Set the paths to the input raster files
demPath = Path("/***/***/***.tif")          
crPath =  Path("/***/***/***.tif")   




#############################################################################
########### Parameters for the simulations - all combinations are explored ##
#############################################################################

# note:  simulations assume an hydraulic conductivity at the ground surface  K_G = 1e-7m/s ( assigned in terms of the ln-conductivity MEAN_Y= ln(K_G) = -16.12 )
# note:  R values here refer to dimensional potential recharges (m/year) and correspond to adimensional recharge R/K_G of [0.001	0.002154435	0.004641589	0.01	0.021544347	0.046415888	0.1	0.215443469	0.464158883	0.63419584	0.824454592	1
] 

# dimensional potential recharge r (m/year)
R = [0.0031536,	0.006794225, 0.014637715, 0.031536,	0.067942252, 0.146377145, 0.31536, 0.679422524,	1.463771455, 2.0, 2.6, 3.1536] 

# formation depth  D_i (m)
IMP_DEPTH =[10, 100, 1000]                                     

# mean and variance of the ln-hydraulic conductivity at the ground surface (ln K_G). Corresponding to K_G = 1e-7m/s
MEAN_Y = [-16.12]                                                 
VAR_Y = [0]

# Topographic stretch factor T_f (-)
TOPOGRAPHY_FACTOR = [0.25, 1, 4]                               

# Exponential depth-decay rate of the hydraulic conductivity (m-1)
ALPHA = [0.0001, 0.001, 0.01]                             
    
# Fixed porosity 
porosity=0.1                          
      
# Number of layers used in the numerical discetization             
layer_number=100                             

# Maximum number of particles assigned for each uppermost saturated cell during the particle tracking
num_particle_per_cell_max=4  

# Below this recharge the maximum number of particles for each cell is reduced to 2  to reduce the computational burden    
R_particle_split=0.679       
 
# Number of bins used to compute the histograms of the groundwater age
number_bins_gw_age_hist=25          



##########################
######### Main function ##
##########################

def MODELLO(crData, demData,IMP_DEPTH,TOPOGRAPHY_FACTOR,ALPHA,MEAN_Y,VAR_Y,R):     
    
        root_path= os.path.dirname(os.path.abspath(__file__))
        os.chdir(root_path)
        path= os.path.join(root_path,'output_parallelizzato',str(R))    # a folder named 'output_parallelizzato' must exist in the same directory as this script
        
        #deletes the folder if exists and creates a new one
        if os.path.exists(path):
         shutil.rmtree(path)
         
        os.mkdir(path)
        os.chdir(path)
        
        
        ###############################
        #Initialize Modflow Nwt solver 
        ###############################
        
        modelname = "model_parallelizzato" 
        modelpath = path
        
        
        # Initialization of the solver
        mf1 = flopy.modflow.Modflow(modelname, exe_name= exeMODFLOW, version="mfnwt", model_ws=modelpath)
        nwt = flopy.modflow.ModflowNwt(mf1 , maxiterout=15000,  maxitinner=10000, mxiterxmd = 10000, headtol=0.001, fluxtol=R/50*3600*24, linmeth=1, stoptol=1e-10, hclosexmd =1e-3, dbdtheta = 0.5, backflag=1, msdr=25, thickfact=1e-04)

    


        #########################
        # spatial discretization
        #########################

        nrow = demDs.RasterYSize
        ncol = demDs.RasterXSize
        delr = geot[1]
        delc = abs(geot[5])
        
        
        demData_stretched= TOPOGRAPHY_FACTOR * demData + 1
        crData[crData<0]=0
        demData_stretched[demData_stretched<0]=0
        
        
        ztop = np.ones((nrow , ncol)) *   np.max(demData_stretched[crData>0])       
        zbot = np.ones((nrow , ncol)) *   np.min(demData_stretched[crData>0]) - IMP_DEPTH
        nlay = layer_number
        delv = (ztop - zbot) / nlay
        botm = np.linspace(ztop, zbot, nlay + 1)
         
          
    
        ###############################
        #definition of flow packages ##
        ###############################
        
        #creates a homogeneus 3-d arrya as big a the domain with the surface hydraulic conductivity
        hk = np.zeros((nlay,nrow,ncol),dtype=np.float64)
        hk= np.exp(np.sqrt(VAR_Y)* hk + MEAN_Y)  
        hk= hk*3600*24  # goes to m/day
        
        #########################################################################
        # reduces the conductivity with depth. The array "reduction_factor_Ks" contain the factors applied to hk, which decreases exponentially. The lower limit is necessary to avoid numeric instabilities
        #######################################################################
        reduction_factor_Ks = np.ones(hk.shape, dtype=np.float32)
        for idx1 in range(nrow):    
         for idx2 in range(ncol):
          for idx3 in range(nlay):   
           if demData_stretched[idx1,idx2] >= botm[idx3+1,idx1,idx2] and demData_stretched[idx1,idx2]>0:    
            reduction_factor_Ks[idx3,idx1,idx2] = 0.001 + (1-0.001) * np.exp(- ALPHA * (demData_stretched[idx1,idx2] - ( botm[idx3,idx1,idx2] + botm[idx3+1,idx1,idx2] )/2  )  )

        # applies the reduction to the initial homogeneous hydraulic conductivity
        hk= np.multiply(hk,reduction_factor_Ks)
        
       ######################################################################################## 
        
        laytyp=np.ones((nlay), dtype=int)
        
        
        # Variables for the DIS package
        dis = flopy.modflow.ModflowDis(mf1, nlay,nrow,ncol,delr=delr,delc=delc,top=ztop,botm=botm[1:],itmuni=4) # <------itmuni=1: seconds ; itmuni = 4: days  
        
        # Variables for the BAS package
        iboundData = np.zeros(demData.shape, dtype=np.int32)
        iboundData[crData > 0 ] = 1
        
        
        #initial conditions for the first try of the solver  - alternative options can be uncommented
        
        strt= demData_stretched * 0.5 +100
        #strt= zbot + IMP_DEPTH[iter_1] + (ztop - zbot - IMP_DEPTH[iter_1]) * R[iter_5]/365/86400 / np.mean(hk) 
        #strt= demData_stretched * R[iter_5]/365/86400 / np.mean(hk) + 1
        #strt= demData_stretched * 0.5 + 200
        #strt= demData_stretched * 1/(1 + np.exp(5-10*R[iter_5]/365/86400 / np.mean(hk))) + 100
        #strt= zbot + Imp_depth - 10
        #strt= demData_stretched
        
        
        
        ######################################
        # Add BAS package to the MODFLOW model
        ######################################
        bas = flopy.modflow.ModflowBas(mf1,ibound=iboundData,strt=strt, hnoflo=-2.0E+020)   
        
        ######################################
        # Add UPW package to the MODFLOW model
        ######################################
        upw = flopy.modflow.ModflowUpw(mf1, laytyp = laytyp, hk = hk, ipakcb=53, hdry = -9999 , iphdry = 1) # <----!!!!! hdry = -1 , iphdry = 1)   #  <---- IMPORTANT! defines how the cells that dry because of the lowering water table are managed 
        
        ###################################################
        #Add the recharge package (RCH) to the MODFLOW model
        ###################################################
        rch_array = np.zeros((nrow, ncol), dtype=np.float32)     
        rch_array[crData>0]=R/365      #m/day
        rch_data = {0: rch_array}
        rch = flopy.modflow.ModflowRch(mf1, nrchop=3, rech =rch_data)

        

        ###################################################
        #Add the drain package (DRN) to the MODFLOW model 
        #Drainage conditions (potential seepage faces) are assigned on the ground surface
        ###################################################
        sorgenti = np.zeros(demData.shape, dtype=np.int32)
        sorgenti[crData >0 ] = 1
        lista = []
        for i in range(sorgenti.shape[0]):
            for q in range(sorgenti.shape[1]):
                
                for j in range(nlay):
                 if   demData_stretched[i,q] < botm[j,i,q] and demData_stretched[i,q] > botm[j+1,i,q] and sorgenti[i,q]>0: 
                  w=j
                  layer_pc[i,q]=j  #assagns to each cell the layer to where the ground surface lays
                
                
                if sorgenti[i,q] == 1:
                    C=1 *3600*24   
                    lista.append([w,i,q,demData_stretched[i,q], C ]) 
        rivDrn = {0:lista}
        
        drn = flopy.modflow.ModflowDrn(mf1,ipakcb=53, stress_period_data=rivDrn, filenames=None)
        
        

        flopy.modflow.ModflowOc(mf1, stress_period_data={(0, 0): ['save head', 'save budget',  'print head']})
        
        
        
        #########################################
        ####### writes files and run simulation##
        ##########################################
        
        #Write input files -> write file with extensions
        mf1.write_input()
        
        #run model -> gives the solution
        mf1.run_model()
    
    
        #####################
        #read output files##
        #####################
        
        fname = os.path.join(modelpath, modelname + ".hds")
        hds = bf.HeadFile(fname)
        times = hds.get_times()
        
        fname = os.path.join(modelpath, modelname + ".cbc")
        cbb = bf.CellBudgetFile(fname)
        kstpkper_list = cbb.get_kstpkper()
        
        drain_fluxes_3D= cbb.get_data(kstpkper=(0,0), text='DRAIN',full3D=True)
        drain_fluxes_2D= np.sum(drain_fluxes_3D[0], axis=0)
        drain_fluxes_2D=drain_fluxes_2D/(delc*delr)          # <------- groes from L^3/T to L/T (from discharge to flux)
        drain_fluxes_2D[crData==0]=np.NaN
        drain_fluxes_minus_recharge_2D=drain_fluxes_2D + rch_array 
        

        #direct runoff, is the rechage that does not enter the deep groundwater and becomes quick runoff 
        drain_fluxes_directrunoff_2D=np.copy(drain_fluxes_2D)
        drain_fluxes_directrunoff_2D[drain_fluxes_directrunoff_2D<-R/365]=-R/365  

        
        #extracts heads 
        fname = os.path.join(modelpath, 'model_parallelizzato.hds')
        hdobj = flopy.utils.HeadFile(fname)
        head = hdobj.get_data()
        head_0=head[-1,:,:]  #<-------- 
        
        #groundwater volume (km^3)
        GW_vol = np.sum(head>-9999) * delc * delr * delv[0,0] * porosity / 1000**3
    
    
    
        #######################################################################
        #######################################################################
        #### particle tracking with modpath 7 #################################
        #######################################################################
        #######################################################################
        

        ########################################################################
        # creates uniformely distributed particles -- tracking FORWARD #########
        ########################################################################
                
        plocs = []
        pids  = []
        localx= []
        localy= []
        localz= []
        particle_count=0
        
        #assigns a number of particles proportional to the actual recharge -  particles are assigned at an height corresponding to the hydraulic head   
        for idx1 in range(nrow):     
         for idx2 in range(ncol):   
           if crData[idx1,idx2] >0 and drain_fluxes_minus_recharge_2D[idx1,idx2]>0:   
              
               R_m_day=R/365
              
               if R>=R_particle_split:
                   
                   num_particle_per_cell_max_case= 4
                   
                   if  (drain_fluxes_minus_recharge_2D[idx1,idx2]  >  4/5 * R_m_day) & (crData[idx1,idx2]==1):
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                 
                       localx.append(0.25)  # relative position within each cell where each particle is released 
                       localx.append(0.25)
                       localx.append(0.75)
                       localx.append(0.75)
                   
                       localy.append(0.25)
                       localy.append(0.75)
                       localy.append(0.25)
                       localy.append(0.75)
                   
                       localz.append(0.95)  # 0 is the bottom of the cell 1 is the top 
                       localz.append(0.95)
                       localz.append(0.95)
                       localz.append(0.95)
                   
                       pids.append(particle_count)              
                       particle_count+=1  
                       pids.append(particle_count)
                       particle_count+=1 
                       pids.append(particle_count)              
                       particle_count+=1  
                       pids.append(particle_count)
                       particle_count+=1

                   if  (drain_fluxes_minus_recharge_2D[idx1,idx2]  >  3/5 * R_m_day)  &  (drain_fluxes_minus_recharge_2D[idx1,idx2]  <  4/5 * R_m_day)  & (crData[idx1,idx2]==1) :
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                 
                       localx.append(0.33)   # relative position within each cell where each particle is released 
                       localx.append(0.66)
                       localx.append(0.5)
           
                       localy.append(0.33)
                       localy.append(0.33)
                       localy.append(0.66)
  
                       localz.append(0.95)  # 0 is the bottom of the cell 1 is the top
                       localz.append(0.95)
                       localz.append(0.95)
                   
                       pids.append(particle_count)              
                       particle_count+=1  
                       pids.append(particle_count)
                       particle_count+=1 
                       pids.append(particle_count)              
                       particle_count+=1  

                   if  (drain_fluxes_minus_recharge_2D[idx1,idx2]  >  2/5 * R_m_day)  &  (drain_fluxes_minus_recharge_2D[idx1,idx2]  <  3/5 * R_m_day)  & (crData[idx1,idx2]==1) :
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                  
                       localx.append(0.33)  # relative position within each cell where each particle is released 
                       localx.append(0.66)
           
                       localy.append(0.33)
                       localy.append(0.66)
  
                       localz.append(0.95) # 0 is the bottom of the cell 1 is the top
                       localz.append(0.95)
                   
                       pids.append(particle_count)              
                       particle_count+=1  
                       pids.append(particle_count)
                       particle_count+=1 
                       
                   if  (drain_fluxes_minus_recharge_2D[idx1,idx2]  >  1/5 * R_m_day)  & ( drain_fluxes_minus_recharge_2D[idx1,idx2]  <  2/5 * R_m_day)  & (crData[idx1,idx2]==1):
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                  
                       localx.append(0.5)  
           
                       localy.append(0.5)
  
                       localz.append(0.95)  
                   
                       pids.append(particle_count)              
                       particle_count+=1                    
                       
                       
                       
                       
               elif R<R_particle_split:
                   
                   num_particle_per_cell_max_case= 2
                   
                   if  (drain_fluxes_minus_recharge_2D[idx1,idx2]  >  2/3 * R_m_day)  & (crData[idx1,idx2]==1):
                       
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
    

                       localx.append(0.25)  
                       localx.append(0.75)
    
                       localy.append(0.25)
                       localy.append(0.75)
        
                       localz.append(0.95)  
                       localz.append(0.95)
                   
                       pids.append(particle_count)              
                       particle_count+=1  
                       pids.append(particle_count)
                       particle_count+=1 

                   if  (drain_fluxes_minus_recharge_2D[idx1,idx2]  >  1/3 * R_m_day)  &   (drain_fluxes_minus_recharge_2D[idx1,idx2]  <  2/3 * R_m_day) & (crData[idx1,idx2]==1):
                       
                       plocs.append((layer_pc[idx1, idx2], idx1, idx2))
                  
                       localx.append(0.5)  
    
                       localy.append(0.5)
        
                       localz.append(0.95)  
                   
                       pids.append(particle_count)              
                       particle_count+=1  


           
        part0 = flopy.modpath.ParticleData(plocs, drape=1, structured=True, particleids=pids, localx=localx, localy=localy, localz=localz)  # drape=1: puts the particle in the first active layer below the one assigned in "plocs". drape=0 places the particle in the layer specificed by "plocs", if the layer is not active the particles are eliminated
        pg0   = flopy.modpath.ParticleGroup(particlegroupname='PG1', particledata=part0,filename='ex01a.pg1.sloc')
      
        particlegroups = [pg0]
        
        
        # default iface for MODFLOW-2005 and MODFLOW 6
        defaultiface = {'RECHARGE': 6, 'ET': 6, 'DRN': 6}     #<---- defines where conditions are applied (0: whole cell, 1-5: lower and lateral faces, 6:upper face) - if seepage faces are simulated it is important to assign the drainage on the upper face 
        defaultiface6 = {'RCH': 6, 'EVT': 6, 'DRN': 6}
         
        
        # create modpath files
        
        mp = flopy.modpath.Modpath7(modelname=modelname + '_mp', flowmodel=mf1, exe_name = exeMODPATH, model_ws=modelpath)
        
        mpbas = flopy.modpath.Modpath7Bas(mp, porosity=porosity,defaultiface=defaultiface)
        
        mpsim = flopy.modpath.Modpath7Sim(mp, simulationtype  ='combined',
                                          trackingdirection   ='forward',     
                                          weaksinkoption      ='stop_at',
                                          weaksourceoption    ='stop_at',
                                          budgetoutputoption  ='summary',
                                          budgetcellnumbers   =None,
                                          traceparticledata   =None,
                                          referencetime       =[0, 0, 0.],
                                          stoptimeoption      ='extend',
                                          timepointdata       =[500, 1000.],
                                          zonedataoption      ='off', 
                                          zones               =None,
                                          particlegroups      =particlegroups)


        # write modpath datasets
        mp.write_input()
        
        # run modpath
        mp.run_model()
        
        
        #####################################################################
        # get pathline file (bulky file, to be commented if not necessary)###
        #####################################################################
        #import flopy.utils.modpathfile as mpf                             
        pthobj = flopy.utils.PathlineFile(modelname + '_mp'+'.mppth')
        p = pthobj.get_alldata()          # pathfile for all particles
        #p1 = pthobj.get_data(partid=1)   # pathfile for a specified particle
        
        
        #counts the particles which infiltrate (GW) and those that are excluded (SW)
        num_particles_GW= len(p)
        num_particles_SW= np.sum(crData>0) * num_particle_per_cell_max_case  - num_particles_GW 
        
        # Computes the length and the velocity of the particles based on the pathfile object
        flowpath_lengths=np.zeros( (num_particles_GW,), dtype=np.float32)
        flowpath_age_mean=np.zeros( (num_particles_GW,), dtype=np.float32)
        flowpath_velocities=[]
        flowpath_velocities_mean=[]
        
        gw_age=[]
        gw_age_weights=[]


        for i in range(num_particles_GW):                                                         # iterates over all particles
            coord_array         = np.array([p[i].x[:-1], p[i].y[:-1], p[i].z[:-1], p[i].time[:-1]]).T       # creates an nx3 array with the vertexes of the flowpath trajectories
            delta_length        = np.sum(np.sqrt(np.diff(coord_array[:,0:3], axis=0)**2), axis=1)        # computes the length of each flowpath
            flowpath_lengths[i] = np.sum(delta_length)
            delta_time          = np.diff(coord_array[:,3], axis=0)
            flowpath_velocities.append(delta_length/delta_time)
            

            flowpath_age_mean[i]= np.sum(  (coord_array[0:-1,3] * delta_length / flowpath_velocities[-1])  [np.isfinite(1/flowpath_velocities[-1])]   )    /   np.sum(  1/ flowpath_velocities[-1][np.isfinite(1/flowpath_velocities[-1])] * delta_length[np.isfinite(1/flowpath_velocities[-1])]  )
        
            #computes the age of the groundwater
            gw_age.append( coord_array[0:-1,3][np.isfinite(1/flowpath_velocities[-1])]  )
            gw_age_weights.append( (delta_length / flowpath_velocities[-1]) [np.isfinite(1/flowpath_velocities[-1])]  )
 

        ############################################
        #calculates the statistics of GW age #######
        ############################################
          
        #creates a vector of all the ages of the particles      
        gw_age_array=np.array([])
        gw_age_weights_array=np.array([])
        
        #unroll each list into a single array (each list is a flowpath)
        for i in range(len(gw_age)): 
            gw_age_array= np.append(gw_age_array, gw_age[i]) 
            gw_age_weights_array= np.append(gw_age_weights_array, gw_age_weights[i]) 
        
        #when there are no particles assigns a dummy value   
        if len(gw_age_array)==0: 
            gw_age_array=np.array([np.nan])
            gw_age_weights_array=np.array([np.nan])
            
            
        #defines the weighted moments
        def weighted_mean(var, wts):
            return np.average(var, weights=wts)
        def weighted_variance(var, wts):
            return np.average((var - weighted_mean(var, wts))**2, weights=wts)
        def weighted_skew(var, wts):
            return (np.average((var - weighted_mean(var, wts))**3, weights=wts) / weighted_variance(var, wts)**(1.5))
        def weighted_kurtosis(var, wts):
            return (np.average((var - weighted_mean(var, wts))**4, weights=wts) / weighted_variance(var, wts)**(2))
            
        gw_age_mean   =  weighted_mean(gw_age_array,gw_age_weights_array)
        gw_age_var    =  weighted_variance(gw_age_array,gw_age_weights_array)
        gw_age_skew   =  weighted_skew(gw_age_array,gw_age_weights_array)
        gw_age_kurt   =  weighted_kurtosis(gw_age_array,gw_age_weights_array)
            
        gw_age_moments= np.array([gw_age_mean, gw_age_var, gw_age_skew, gw_age_kurt])
        
        #computes histogram 
        hist_gw_age_natural=[]
        hist_gw_age_log=[]

        hist, bins, _ = plt.hist(gw_age_array, bins=number_bins_gw_age_hist, range= (  10, np.nansum([100, np.max(gw_age_array)])  ),density=True, weights=gw_age_weights_array)
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        hist_log, bins_log, _ =plt.hist(gw_age_array, bins=logbins, density=True, weights=gw_age_weights_array)
        bins_width=bins[1:]-bins[0:-1]
        bins_width_log=bins_log[1:]-bins_log[0:-1]
        hist_gw_age_natural.append([hist,  bins[:-1],  bins_width])
        hist_gw_age_log.append([hist_log,  logbins[:-1],  bins_width_log])

        

    
        
        # calculates the surface and subsurface fluxes in two different ways (either using the modpath particles or the modflow fluxes)
        Q_gw_particle= R * delr * delc / num_particle_per_cell_max_case  *  num_particles_GW  /365/24/3600  # m3/sec
        Q_sw_particle= R * delr * delc / num_particle_per_cell_max_case  *  num_particles_SW  /365/24/3600  # m3/sec
        Q_gw_Q_sw_ratio_particle=Q_gw_particle/Q_sw_particle
        
        Q_gw_fluxes= np.nansum(drain_fluxes_minus_recharge_2D[drain_fluxes_minus_recharge_2D>0]) * delr * delc/24/3600
        Q_sw_fluxes= -1*np.nansum(drain_fluxes_directrunoff_2D) * delr * delc /24/3600
        Q_gw_Q_sw_ratio_fluxes=Q_gw_fluxes/Q_sw_fluxes
        
        Q_tot=R /365/24/3600 * catchment_area *1000**2
        Q_gw_fluxes_normalized=Q_gw_fluxes/Q_tot
        Q_sw_fluxes_normalized=Q_sw_fluxes/Q_tot
        
        ###############################
        # get travel times  ###########
        ###############################
        
        endobj = flopy.utils.EndpointFile(modelname + '_mp'+'.mpend')   # gets all the traveltimes
        e = endobj.get_alldata()
        traveltime=e.time
        #e1 = endobj.get_data(partid=1)   # gets a single traveltime
        

    
        #########################################################################################
        ####### end particle tracking ###########################################################
        #########################################################################################

        return traveltime, flowpath_lengths, drain_fluxes_2D, drain_fluxes_minus_recharge_2D, drain_fluxes_directrunoff_2D, GW_vol, Q_gw_Q_sw_ratio_particle, Q_gw_Q_sw_ratio_fluxes, Q_tot, Q_gw_fluxes_normalized, Q_sw_fluxes_normalized,  gw_age_moments, hist_gw_age_natural, hist_gw_age_log   



############################################
#####################################################
# MAIN ##############################################
#####################################################
############################################



#Open and read raster files 
demDs =gdal.Open(demPath)
crDs = gdal.Open(crPath)
geot = crDs.GetGeoTransform() 

# Get data as arrays
demData_original = demDs.GetRasterBand(1).ReadAsArray()
crData = crDs.GetRasterBand(1).ReadAsArray()

demData=np.array(demData_original, copy=True)  

# shifts the DTM so that the outlet has elevation 0
demData[crData>0]=demData[crData>0]-np.min(demData[crData>0])  

# Get statistics
stats = demDs.GetRasterBand(1).GetStatistics(0,1) 
stats

catchment_area= np.sum(crData>0) * geot[1] * abs(geot[5]) / 1000**2  # in km^2


#########################################
## allocates the memory for the variables
########################################

run_count=0
risultati=[]

drain_fluxes_2D_ARRAY                = np.zeros( (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA),) + demData.shape, dtype=np.float32)                # fluxes from the seepage faces (drainages)
drain_fluxes_minus_recharge_2D_ARRAY = np.zeros( (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA),) + demData.shape, dtype=np.float32) # effective recharge
drain_fluxes_directrunoff_2D_ARRAY   = np.zeros( (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA),) + demData.shape, dtype=np.float32)   # quick runoff 

R_K_ratio                = np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
Q_gw_Q_sw_ratio_particle =np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
Q_gw_Q_sw_ratio_fluxes   = np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
Q_tot                    = np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
Q_gw_fluxes_normalized   = np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
Q_sw_fluxes_normalized   = np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
GW_volume                = np.zeros( ((1),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)


GW_age_moments = np.zeros( ((4),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)

GW_age_hist_natural =[]
GW_age_hist_log     =[]

max_particle_nubmer    = num_particle_per_cell_max * sum(sum(crData>0))                   
traveltime_ARRAY       = -9999*np.ones( ((max_particle_nubmer),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
streamflow_age_ARRAY   = -9999*np.ones( ((max_particle_nubmer),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)
flowpath_lengths_ARRAY = -9999*np.ones( ((max_particle_nubmer),+ (len(R)*len(MEAN_Y)*len(IMP_DEPTH)*len(TOPOGRAPHY_FACTOR)*len(ALPHA))), dtype=np.float32)        

layer_pc=np.zeros((crData.shape[0], crData.shape[1]),dtype=np.int32)  # layer where the ground suface lays


##################################################################
## MAIN continues: runs the simulations and builds the outputs ###
##################################################################

for iter_1 in range(len(IMP_DEPTH)):
    
 for iter_2 in range(len(TOPOGRAPHY_FACTOR)): 

  for iter_3 in range(len(ALPHA)):   

   for iter_4 in range(len(MEAN_Y)): 
  
    # initializes parallel computations
    def main():
      pool      = multiprocessing.Pool(len(R))  # sets the number of cores to be used in the parallel computations, in this case is equal to the assigned recharge values  - reduce the number if you have less cores 
      risultati = pool.starmap(MODELLO, [(crData, demData,IMP_DEPTH[iter_1],TOPOGRAPHY_FACTOR[iter_2],ALPHA[iter_3],MEAN_Y[iter_4],VAR_Y[iter_4],ricarica) for ricarica in R])

      pool.close()
      pool.join()

      return risultati

        
    # better protect your main function when you use multiprocessing
    if __name__ == '__main__':
     
    
    
     risultati = risultati + main()
    
     
    
#builds the arrays of the outputs starting from the lists in "risultati"

for i in range(len(risultati)):
 traveltime_ARRAY[0:len(risultati[i][0]),i]        =risultati[i][0]
 flowpath_lengths_ARRAY[0:len(risultati[i][1]),i]  =risultati[i][1]
 drain_fluxes_2D_ARRAY[i,:,:]                      =risultati[i][2]
 drain_fluxes_minus_recharge_2D_ARRAY[i,:,:]       =risultati[i][3]
 drain_fluxes_directrunoff_2D_ARRAY[i,:,:]         =risultati[i][4]
 GW_volume[0,i]                                    =risultati[i][5]
 Q_gw_Q_sw_ratio_particle[0,i]                     =risultati[i][6]
 Q_gw_Q_sw_ratio_fluxes[0,i]                       =risultati[i][7]   
 Q_tot[0,i]                                        =risultati[i][8]  
 Q_gw_fluxes_normalized[0,i]                       =risultati[i][9]  
 Q_sw_fluxes_normalized[0,i]                       =risultati[i][10]  
 
 GW_age_moments[:,i,None]=risultati[i][11][:,None] 

 GW_age_hist_natural.append(risultati[i][12][0])
 GW_age_hist_log.append(risultati[i][13][0])

        
##################################################################################
#######end model runs#############################################################
##################################################################################






###########################################################
##### uncomment if you want to save model output variables
###########################################################

# =============================================================================
# import shelve
# 
# #where you want to save the results
# filename='E:/.../.../shelve.out'
# my_shelf = shelve.open(filename,'n') # 'n' for new
# 
# my_shelf['R']          = globals()['R']
# my_shelf['IMP_DEPTH']  = globals()['IMP_DEPTH']
# my_shelf['MEAN_Y']     = globals()['MEAN_Y']
# my_shelf['VAR_Y']      = globals()['VAR_Y']
# my_shelf['TOPOGRAPHY_FACTOR'] = globals()['TOPOGRAPHY_FACTOR']
# my_shelf['ALPHA']      = globals()['ALPHA']
# 
# my_shelf['demData_original']                    = globals()['demData_original']
# my_shelf['drain_fluxes_2D_ARRAY']               = globals()['drain_fluxes_2D_ARRAY']
# my_shelf['drain_fluxes_directrunoff_2D_ARRAY']  = globals()['drain_fluxes_directrunoff_2D_ARRAY']
# my_shelf['drain_fluxes_minus_recharge_2D_ARRAY']= globals()['drain_fluxes_minus_recharge_2D_ARRAY']
# my_shelf['traveltime_ARRAY']                    = globals()['traveltime_ARRAY']
# my_shelf['streamflow_age_ARRAY']                = globals()['streamflow_age_ARRAY']
# my_shelf['flowpath_lengths_ARRAY']              = globals()['flowpath_lengths_ARRAY']
# my_shelf['R']                                   = globals()['R']
# my_shelf['ALPHA']                               = globals()['ALPHA']
# my_shelf['TOPOGRAPHY_FACTOR']                   = globals()['TOPOGRAPHY_FACTOR']
# my_shelf['IMP_DEPTH']                           = globals()['IMP_DEPTH']
# my_shelf['porosity']                            = globals()['porosity']
# my_shelf['GW_volume']                           = globals()['GW_volume']
# 
# my_shelf['catchment_area']           = globals()['catchment_area']
# my_shelf['Q_tot']                    = globals()['Q_tot']
# my_shelf['Q_gw_fluxes_normalized']   = globals()['Q_gw_fluxes_normalized']
# my_shelf['Q_sw_fluxes_normalized']   = globals()['Q_sw_fluxes_normalized']
# my_shelf['Q_gw_Q_sw_ratio_fluxes']   = globals()['Q_gw_Q_sw_ratio_fluxes']
# my_shelf['Q_gw_Q_sw_ratio_particle'] = globals()['Q_gw_Q_sw_ratio_particle']
# my_shelf['GW_volume']                = globals()['GW_volume']
# 
# my_shelf['GW_age_moments']      = globals()['GW_age_moments']
# my_shelf['GW_age_hist_natural'] = globals()['GW_age_hist_natural']
# my_shelf['GW_age_hist_log']     = globals()['GW_age_hist_log']
# 
# my_shelf.close()
# =============================================================================


####################################################################
##### uncomment if you want to load the saved model output variables
####################################################################

# =============================================================================
# import shelve
#
# filename='E:/.../.../shelve.out'
# my_shelf = shelve.open(filename)
# for key in my_shelf:
#     globals()[key]=my_shelf[key]
# my_shelf.close()
# =============================================================================
