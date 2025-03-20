import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import flopy as fp
from pathlib import Path
import rasterio



model_name = "maso_transiente"


# PARAMETERS
ks               =  1   # m/day
porosity         =  0.1
specific_yield   =  porosity
n_lay            =  50
aq_tickness      =  200
depth_bottom     =  100
drain_multiplier =  1
recharge         =  0.003     # m/day


# Transient simulation settings
nper   =  2            # number of stress periods
perlen = [1, 36000]     # length of each stress period - day
nstp   = [1, 50]       # number of time steps per period
tsmult = [1.0] * nper  # time step multiplier

number_timesteps_to_save = 50  # must be lower than nstp
total_timesteps          = np.sum(np.array(nstp))
time_steps_to_save       = (np.linspace(1,total_timesteps-1, number_timesteps_to_save ).astype(int)).tolist()



# location of MODFLOW 6 executable
path_modflow = Path("C:/DEV/Exe/mf6.exe")

# DTM 
path_dtm = Path("E:/Trento/modello_transiente/DTM/DEM_Maso_100m.tif")

# Define the directory path for input and output files
base_dir = Path("E:/Trento/modello_transiente/simultions/maso_transiente")


#if True use a DTM, otherwise use sample topography
if True:
    with rasterio.open(path_dtm) as src:
        nodata      = src.nodata
        transform   = src.transform
        pxl_size    = transform[0]
        dtm         = src.read(1)
        dtm[dtm==nodata] = np.nan
        catchment   = ~np.isnan(dtm)      
else:
    #SIMPLIFIED TEST DTM
    dim_i   = 51  ; dim_j   = 51
    slope_i = 0   ; slope_j = 0.05
    dtm= np.ones((dim_i,dim_j))
    i, j = np.mgrid[0:dtm.shape[0], 0:dtm.shape[1]]
    catchment = ~np.isnan(dtm)
    pxl_size = 100
    dtm = dtm * i * slope_i*pxl_size + dtm * np.abs( (dim_j)/2 - (j+0.5) ) * slope_j*pxl_size




n_row, n_col     = dtm.shape 
i_catch, j_catch = np.where(catchment==1)


C   =  ks * pxl_size**2 *drain_multiplier   # drain conductance

#shift to elevation = 0 at the outlet
dtm= dtm - np.nanmin(dtm) 
z_max = np.nanmax(dtm) 
z_min = np.nanmin(dtm) 


if False:
    plt.imshow(dtm)
    plt.colorbar()
    plt.show()
    plt.imshow(catchment)


# CREATE SIMULATION OBJECT
sim = fp.mf6.MFSimulation(sim_name =  model_name, 
                          exe_name =  path_modflow, 
                          version  =  "mf6", 
                          sim_ws   =  base_dir,
                          memory_print_option = 'SUMMARY'
                          )


#TEMPORAL DISCRETIZATION
tdis = fp.mf6.ModflowTdis(sim,
                          pname       =  "tdis", 
                          time_units  =  "DAYS", 
                          nper        =  nper,                # number of stress periods
                          perioddata=list(zip(perlen, nstp, tsmult))   # length of stress period, number of timesteps in the stress period, timestep multiplier
                          )


#SOLVER SETTINGS
ims = fp.mf6.ModflowIms(sim, 
                        pname             = "ims",
                        complexity        = "COMPLEX",
                        print_option      = "SUMMARY",
                        inner_dvclose     = 1e-3,
                        outer_dvclose     = 1e-2,
                        inner_maximum     = 10000,
                        outer_maximum     = 10000,
                        #under_relaxation  = 'SIMPLE', under_relaxation_gamma = 0.5
                        #under_relaxation  = 'DBD', under_relaxation_theta = 0.7, under_relaxation_kappa = 0.01, under_relaxation_momentum = 0.001
                        )

   


# CREATE GROUNDWATER FLOW MODEL OBJECT
gwf = fp.mf6.ModflowGwf(sim, 
                        modelname      = model_name,
                        model_nam_file = f"{model_name}.nam",
                        save_flows     = True,
                        newtonoptions  = 'under_relaxation'   # <-----------------
                        )


dz           = (z_max + depth_bottom) / n_lay
cell_bottom  = np.linspace( z_max - dz ,  - depth_bottom, n_lay )
top          = np.ones_like(dtm) * z_max 
layer_limits = np.insert(cell_bottom, 0, z_max )
aspect_ratio = dz/pxl_size

cell_bottoms_3d = np.tile(cell_bottom[:, np.newaxis, np.newaxis], (1, n_row, n_col))


#SPECIFY ACTIVE CELLS
idomain = np.zeros([n_lay, n_row, n_col], dtype= np.uint8)
#idomain[:, catchment]  = 1 



surface_layer        = np.zeros([n_lay, n_row, n_col], dtype= np.uint8)
surface_layer_number = np.zeros([n_row, n_col], dtype= np.uint32)


for i in range(len(i_catch)):
    for j in range(n_lay):
        if  cell_bottom[j]  <= dtm[i_catch[i],j_catch[i] ]  :
            idomain[j, i_catch[i], j_catch[i]]  =  1
            if cell_bottom[j] + dz  >= dtm[i_catch[i],j_catch[i] ] :
                surface_layer[j,i_catch[i],j_catch[i]] =  1    
                surface_layer_number[i_catch[i],j_catch[i]] =  j



#SPATIAL DISCRETIZATION
dis = fp.mf6.ModflowGwfdis(
                           gwf,
                           idomain      =  idomain,    # 1 active cell, 0 inactive cell
                           length_units = 'METERS',
                           nlay         =  n_lay,
                           nrow         =  n_row,
                           ncol         =  n_col,
                           delr         =  pxl_size,
                           delc         =  pxl_size,
                           top          =  top,
                           botm         =  cell_bottom
                           )



# DEFINE ICELLTYPE FOR EACH LAYER
icelltype   = np.copy(idomain)



#rewet_record = [('wetfct', 1, 'iwetit', 1, 'ihdwet', 0)] # <-----------------



#NODE PROPERTY
npf = fp.mf6.ModflowGwfnpf(gwf, 
                           icelltype   =  icelltype,  # 1for unconfined, 0 confined
                           #cvoptions    = ['DEWATERED'],              # <-----------------  se "newtonoptions" disattivata di sopra 
                           #wetdry = 1,  rewet_record=rewet_record,    # <-----------------  se "newtonoptions" disattivata di sopra 
                           k           =  ks,         #k_horizontal
                           k33         =  ks,         #k_vertical
                           save_flows  =  True,
                           save_specific_discharge = True,     
                           )


#INITIAL CONDITIONS
param_scale  = 0.5 #(-)
param_offset = +10 #(m)
h_surf       = dtm*param_scale + param_offset
h_start      = np.ones((n_lay, n_row, n_col))* h_surf

#set fixed head to cells above the water table
for i in range(len(i_catch)):
    for j in range(n_lay):
        if  cell_bottom[j] > h_surf[i_catch[i],j_catch[i] ] :
            h_start[ j,  i_catch[i], j_catch[i] ] =  cell_bottom[j]   # + 0.01   # <----------------- se "newtonoptions" disattivata di sopra 




ic = fp.mf6.ModflowGwfic(gwf,
                         pname  = "ic", 
                         strt   = h_start)





#RECHARGE
recharge_array = catchment * recharge

#recharge_array = surface_layer * recharge


Rch = {0:recharge_array,  1:recharge_array*0}
#Rch = {0:recharge}


rch = fp.mf6.ModflowGwfrcha(
                            gwf, 
                            pname        = "rch", 
                            readasarrays = True,
                            irch         = surface_layer_number, # 0,
                            recharge     = Rch,
                            )


#SEEPAGE with drains
drain=[]

for i in range(len(i_catch)):
    for j in range(n_lay):
        if  layer_limits[j] >= dtm[i_catch[i],j_catch[i] ] > layer_limits[j+1] :
            drain.append([ j,  i_catch[i], j_catch[i]  ,  dtm[i_catch[i], j_catch[i]],  C ])
            


Drn = {0:drain, 1:drain}  

drn = fp.mf6.ModflowGwfdrn(gwf,
                           pname              ="drn", 
                           stress_period_data = Drn
                           )
 



# STORAGE SETTINGS
sto = fp.mf6.ModflowGwfsto(gwf,
                           sy=specific_yield,
                           iconvert=1, 
                           steady_state={0: True},
                           transient   ={1: True})


#OUTPUT CONTROL
head_filerecord   =  f"{model_name}.hds"
budget_filerecord =  f"{model_name}.cbc"
printrecord       = [("HEAD", "LAST")]

#CREATE THE SAVERECORD OPTION
saverecord      = [('HEAD',   'STEPS', ts) for ts in time_steps_to_save]
saverecord.extend([('BUDGET', 'STEPS', ts) for ts in time_steps_to_save])

#create the saverecord option
#saverecord      = [('HEAD',   'ALL'), ('BUDGET', 'ALL')]
#saverecord      = [('HEAD','FIRST'), ('BUDGET','FIRST'),('HEAD','LAST'), ('BUDGET','LAST')]



oc = fp.mf6.ModflowGwfoc(
                         gwf,
                         saverecord        = saverecord,
                         head_filerecord   = head_filerecord,
                         budget_filerecord = budget_filerecord,
                         printrecord       = printrecord,
                         )
    
    
if True:
    #horizontal slice plot
    pmv = fp.plot.PlotMapView(model=gwf, layer=30)
    pmv.plot_bc('DRN', color='b')
    pmv.plot_grid(colors='silver', lw=0.5)
    
    
    #cross section
    plt.figure(figsize = (5,3) )
    pxs = fp.plot.PlotCrossSection (model = gwf , line= { "Row" : 25 })
    pxs.plot_bc('DRN', color='b')
    #pxs.plot_bc('WEL')
    pxs.plot_grid(colors='silver', lw=0.1)
    
    
    
    
    
###########################
###### RUN  ###############
###########################


sim.write_simulation()
success, buff = sim.run_simulation()

if not success:
    raise Exception("MODFLOW 6 did not terminate normally.")
    
    
#///////////////////////////////////////////////////////////////////////    
#POST PROCESSING         /////////////////////////////////////////////// 
#///////////////////////////////////////////////////////////////////////   
    

# GET THE HEAD OBJECT
fname     = base_dir/f'{model_name}.hds'
hdobj     = fp.utils.HeadFile(fname)
times_out = hdobj.times   
Q_timeseries = []


#time_out_to_extract_index = 5


for time_out_to_extract_index, time in enumerate(times_out):
        
    # GET THE HEAD VALUES FOR A TIME STEP AS 3D ARRAY  
    #time      = hdobj.times[time_out_to_extract_index]
    head_3d   = hdobj.get_data(totim=time)   
    head_3d[head_3d==1.e+30] = np.nan
    head_3d[head_3d<=cell_bottoms_3d]= -999
    
    # GET THE HEAD VALUES ALSO FOR THE PREVIOUS TIMESTEP
    if time_out_to_extract_index > 0:
       time_before        = hdobj.times[time_out_to_extract_index-1]
       head_3d_before     = hdobj.get_data(totim=time_before)   
       head_3d_before[head_3d_before==1.e+30] = np.nan
       head_3d_before[head_3d_before<=cell_bottoms_3d]= -999 
    else:
        head_3d_before=np.copy(head_3d)
    
    
    
    
    # MIN-MAX heads
    head_max = np.nanmax(head_3d)
    head_min = np.nanmin(head_3d[head_3d>-999])
    print('head_max: ', head_max)   
    print('head_min: ', head_min)   
    
    
    
    
    
    
    #computes storage and delta storage
    water_table_layer    = n_lay - np.sum(head_3d>-999, axis=0) -1  # the layer where the water table lays
    i_indices, j_indices = np.ogrid[:n_row, :n_col]  # Create open grid for x and y indices
    head_water_table     = head_3d[water_table_layer, i_indices, j_indices]  # head at the water table cell
    bottom_elevation     = np.min(cell_bottoms_3d, axis=0)
    storage              = np.nansum( (head_water_table - bottom_elevation) * porosity * pxl_size**2)
    
    water_table_layer_before    = n_lay - np.sum(head_3d_before>-999, axis=0) -1 
    head_water_table_before     = head_3d_before[water_table_layer_before, i_indices, j_indices]  
    storage_before              = np.nansum( (head_water_table_before - bottom_elevation) * porosity * pxl_size**2)
    
    Delta_storage =  storage - storage_before
    
    
    
    #extracts cell budget file
    cvv  =  fp.utils.CellBudgetFile(base_dir/f'{model_name}.cbc' )
    Q_rch  =  cvv.get_data(text='RCH',totim=time)[0]
    Q_drn  =  cvv.get_data(text='DRN',totim=time)[0]
    
    #creates the 3d budget array
    Q_rch_3d = cvv.get_data(idx=None, totim=time, text='RCH',  full3D=True)[0].data
    Q_drn_3d = cvv.get_data(idx=None, totim=time, text='DRN',  full3D=True)[0].data
    
    
    #flattens
    Q_drn_2d = np.abs( np.sum(Q_drn_3d, axis=0) )
    Q_rch_2d = np.abs( np.sum(Q_rch_3d, axis=0) )
    
    Q_rch_2d = np.sum(Q_rch_3d, axis=0)
    
    
    R = Q_rch['q'].sum()
    Q = Q_drn['q'].sum() 
    
    print(f'R: {R} m^3/d')
    print(f'Q: {Q} m^3/d')
    print(f'Delta_Storage: {Delta_storage} m^3')
    print(f'Relative Error: { (R+Q)/R * 100} %')
    
    Q_timeseries.append(Q)
    
    
    Q_difference =  Q_rch_2d  - Q_drn_2d
    Q_rch_2d_NET =  Q_difference * (Q_difference>0)
    Q_drn_2d_NET =  np.abs( Q_difference * (Q_difference<0) )
    Q_rch_2d_NET[catchment==0] = np.nan
    Q_drn_2d_NET[catchment==0] = np.nan
    
    
    Q_drn_min= np.nanmin(Q_drn_2d_NET[Q_drn_2d_NET>0])
    Q_drn_max= np.nanmax(Q_drn_2d_NET)
    
    if np.sum(Q_rch_2d_NET)>0:
        Q_rch_min= np.nanmin(Q_rch_2d_NET[Q_rch_2d_NET>0])
    else:
        Q_rch_min = 0 
        
    Q_rch_max= np.nanmax(Q_rch_2d_NET)
    
    
    
    #FLOW RATIO
    Q_bf    = np.nansum(Q_drn_2d_NET)  # equals np.nansum(Q_rch_2d_NET)
    Q_dr    = np.nansum(Q_rch_2d) - Q_bf
    Q_ratio = Q_bf/(Q_bf + Q_dr)
    
    
    
    
    
    #################
    ## PLOTS ########
    #################
    
    if False:
            
        contour_interval = 10 
        levels = np.arange(head_min, head_max, contour_interval)
        
        
        #CUSTOM COLORMAP
        base_cmap      = plt.cm.viridis # Create a new colormap with a base of gray
        base_colors    = base_cmap(np.arange(base_cmap.N))
        specific_color = np.array([0.8, 0.8, 0.8, 1.0]) 
        base_colors[0] = specific_color  
        custom_cmap    = mcolors.ListedColormap(base_colors)
        
        
        
        #PLOT CROSS SECTION - ROW direction
        row_cross= 30
        plt.figure(figsize = (5,3) )
        pxs = fp.plot.PlotCrossSection (model = gwf , line= { "Row" : row_cross } ) 
        plot = pxs.plot_array(head_3d,vmin=0-1, vmax=head_max,cmap=custom_cmap)
        pxs.plot_grid(colors='silver', lw=0.1)
        pxs.contour_array( head_3d, levels=levels,  linewidths = .5, colors = 'b')
        plt.title(rf'Heads in Row cross-section {row_cross} - $\Delta_{{contour}}$ = {contour_interval} m')
        cbar = plt.colorbar(plot, orientation='vertical')
        cbar.set_label('Head (m)')  # Label for the colorbar
        plt.show()
        
        #PLOT CROSS SECTION - ROW direction
        col_cross= 40
        plt.figure(figsize = (5,3) )
        pxs = fp.plot.PlotCrossSection (model = gwf , line= { "Column" : col_cross } ) 
        plot = pxs.plot_array(head_3d,vmin=0-1, vmax=head_max,cmap=custom_cmap)
        pxs.plot_grid(colors='silver', lw=0.1)
        pxs.contour_array( head_3d, levels=levels,  linewidths = .5, colors = 'b')
        plt.title(rf'Heads in Column cross-section {row_cross} - $\Delta_{{contour}}$ = {contour_interval} m')
        cbar = plt.colorbar(plot, orientation='vertical')
        cbar.set_label('Head (m)')  # Label for the colorbar
        plt.show()
        
        
        
        #PLOT HORIZONTAL CROSS-SECTION
        layer_cross = n_lay-1
        pmv = fp.plot.PlotMapView(model=gwf,  layer=layer_cross)
        qm  = pmv.plot_array(head_3d,vmin=head_min, vmax=head_max)
        plt.colorbar(qm, shrink = 1, label ='head (m)')
        cs = pmv.contour_array( head_3d, levels=levels,  linewidths = .5, colors = 'b')
        plt.title(rf'Heads in Layer {layer_cross} - $\Delta_{{contour}}$ = {contour_interval} m')
        plt.show()
        
        
        
        
        # =============================================================================
        # ### PLOT NET SEEPAGE  - [L]^3/[T]
        # base_cmap      = plt.cm.viridis # Create a new colormap with a base of gray
        # base_colors    = base_cmap(np.arange(base_cmap.N))
        # specific_color = np.array([0.8, 0.8, 0.8, 1.0]) 
        # base_colors[0] = specific_color  
        # custom_cmap = mcolors.ListedColormap(base_colors)
        # 
        # plt.imshow(Q_drn_2d_NET, interpolation='none', cmap=custom_cmap, vmin=0, vmax=Q_drn_max)
        # plt.colorbar()
        # plt.title('Seepage discharge $[L]^3/[T]$')
        # plt.show()
        # =============================================================================
        
        
        
        ### PLOT NET SEEPAGE  - [L]^3/[T]
        base_cmap      = plt.cm.viridis # Create a new colormap with a base of gray
        base_colors    = base_cmap(np.arange(base_cmap.N))
        specific_color = np.array([0.8, 0.8, 0.8, 1.0]) 
        base_colors[0] = specific_color  
        custom_cmap    = mcolors.ListedColormap(base_colors)
        
        pmv = fp.plot.PlotMapView(model=gwf)
        qm  = pmv.plot_array(Q_drn_2d_NET,vmin=0, vmax=Q_drn_max, cmap=custom_cmap)
        plt.colorbar(qm, shrink = 1, label ='Seepage $[L]^3/[T]$')
        plt.title('Seepage discharge $[L]^3/[T]$')
        plt.show()
        
        
        
        
        
        
        # =============================================================================
        # ### PLOT NET RECHARGE  - [L]^3/[T]
        # Q_min= np.nanmin(Q_rch_2d_NET[Q_rch_2d_NET>0])
        # Q_max= np.nanmax(Q_rch_2d_NET)
        # 
        # base_cmap      = plt.cm.viridis # Create a new colormap with a base of gray
        # base_colors    = base_cmap(np.arange(base_cmap.N))
        # specific_color = np.array([0.8, 0.8, 0.8, 1.0]) 
        # base_colors[0] = specific_color  
        # custom_cmap = mcolors.ListedColormap(base_colors)
        # 
        # plt.imshow(Q_rch_2d_NET, interpolation='none', cmap=custom_cmap, vmin=0, vmax=Q_max)
        # plt.colorbar()
        # plt.title('Recharge discharge $[L]^3/[T]$')
        # plt.show()
        # =============================================================================
        
        
        ### PLOT NET RECHARGE  - [L]^3/[T]
        base_cmap      = plt.cm.viridis # Create a new colormap with a base of gray
        base_colors    = base_cmap(np.arange(base_cmap.N))
        specific_color = np.array([0.8, 0.8, 0.8, 1.0]) 
        base_colors[0] = specific_color  
        custom_cmap    = mcolors.ListedColormap(base_colors)
        
        pmv = fp.plot.PlotMapView(model=gwf)
        qm  = pmv.plot_array(Q_rch_2d_NET,vmin=0, vmax=Q_rch_max, cmap=custom_cmap)
        plt.colorbar(qm, shrink = 1, label ='Recharge $[L]^3/[T]$')
        plt.title('Recharge  $[L]^3/[T]$')
        plt.show()
    
    
    
# PLOT RECESSION
plt.plot(times_out,np.abs(np.array(Q_timeseries))/ 3600/24 )
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('t (days)')
plt.ylabel('Q ($m^3/s$)')
