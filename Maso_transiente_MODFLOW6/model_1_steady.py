import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import flopy as fp
from pathlib import Path
import rasterio



model_name = "maso_transiente"


# PARAMETERS
ks               =  0.01   # m/day
porosity         =  0.1
n_lay            =  100
aq_tickness      =  200
depth_bottom     =  100
drain_multiplier =  1
recharge         =  0.003     # m/day




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
    dim_i   = 50  ; dim_j   = 50
    slope_i = 0   ; slope_j = 0.05
    dtm= np.ones((dim_i,dim_j))
    i, j = np.mgrid[0:dtm.shape[0], 0:dtm.shape[1]]
    catchment = ~np.isnan(dtm)
    pxl_size = 100
    dtm = dtm * i * slope_i*pxl_size + dtm * np.abs( (dim_j-1)/2 - j) * slope_j*pxl_size




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
                          nper        =  1,                # number of stress periods
                          perioddata  =  [(1.0, 1, 1.0)]   # length of stress period, number of timesteps in the stress period, timestep multiplier
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



# Define ICELLTYPE for each layer
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


Rch = {0:recharge_array}
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
            

Drn = {0:drain}  # assigns to stress period 0

drn = fp.mf6.ModflowGwfdrn(gwf,
                           pname              ="drn", 
                           stress_period_data = Drn
                           )
 





#OUTPUT CONTROL
head_filerecord   =  f"{model_name}.hds"
budget_filerecord =  f"{model_name}.cbc"
saverecord  = [("HEAD", "ALL"), ("BUDGET", "ALL")]
printrecord = [("HEAD", "LAST")]



oc = fp.mf6.ModflowGwfoc(
                         gwf,
                         saverecord        = saverecord,
                         head_filerecord   = head_filerecord,
                         budget_filerecord = budget_filerecord,
                         printrecord       = printrecord,
                         )
    
    
if True:
    #horizontal plot
    pmv = fp.plot.PlotMapView(model=gwf, layer=30)
    pmv.plot_bc('DRN', color='b')
    pmv.plot_grid(colors='silver', lw=0.5)
    
    
    #cross section
    plt.figure(figsize = (5,3) )
    pxs = fp.plot.PlotCrossSection (model = gwf , line= { "Row" : 50 })
    pxs.plot_bc('DRN', color='b')
    #pxs.plot_bc('WEL')
    pxs.plot_grid(colors='silver', lw=0.1)
    
    
    
    
    
######################
###### RUN  ###############
#################################
sim.write_simulation()
success, buff = sim.run_simulation()

if not success:
    raise Exception("MODFLOW 6 did not terminate normally.")
    
    
    
# Get the head values for the last time step as 3d array  
fname = base_dir/f'{model_name}.hds'
hdobj = fp.utils.HeadFile(fname)
head_3d     = hdobj.get_data(totim=hdobj.times[-1])   
head_3d[head_3d==1.e+30] = np.nan
head_3d[head_3d<=cell_bottoms_3d]= 0



# MIN-MAX heads
head_max = np.nanmax(head_3d)
head_min = np.nanmin(head_3d)
print('head_max: ', head_max)   
print('head_min: ', head_min)   





contour_interval = 10 
levels = np.arange(head_min, head_max, contour_interval)



#PLOT CROSS SECTION - ROW direction
row_cross= 50
plt.figure(figsize = (5,3) )
pxs = fp.plot.PlotCrossSection (model = gwf , line= { "Row" : row_cross } ) 
plot = pxs.plot_array(head_3d,vmin=head_min, vmax=head_max)
pxs.plot_grid(colors='silver', lw=0.1)
pxs.contour_array( head_3d, levels=levels,  linewidths = .5, colors = 'b')
plt.title(rf'Heads in Row cross-section {row_cross} - $\Delta_{{contour}}$ = {contour_interval} m')
cbar = plt.colorbar(plot, orientation='vertical')
cbar.set_label('Head (m)')  # Label for the colorbar
plt.show()


#PLOT CROSS SECTION - ROW direction
col_cross= 50
plt.figure(figsize = (5,3) )
pxs = fp.plot.PlotCrossSection (model = gwf , line= { "Column" : col_cross } ) 
plot = pxs.plot_array(head_3d,vmin=head_min, vmax=head_max)
pxs.plot_grid(colors='silver', lw=0.1)
pxs.contour_array( head_3d, levels=levels,  linewidths = .5, colors = 'b')
plt.title(rf'Heads in Column cross-section {row_cross} - $\Delta_{{contour}}$ = {contour_interval} m')
cbar = plt.colorbar(plot, orientation='vertical')
cbar.set_label('Head (m)')  # Label for the colorbar
plt.show()



#PLOT HORIZONTAL CROSS-SECTION
layer_cross = 99
pmv = fp.plot.PlotMapView(model=gwf,  layer=layer_cross)
qm  = pmv.plot_array(head_3d,vmin=head_min, vmax=head_max)
plt.colorbar(qm, shrink = 1, label ='head (m)')
cs = pmv.contour_array( head_3d, levels=levels,  linewidths = .5, colors = 'b')
plt.title(rf'Heads in Layer {layer_cross} - $\Delta_{{contour}}$ = {contour_interval} m')
plt.show()




#extracts cell budget file
cvv  =  fp.utils.CellBudgetFile(base_dir/f'{model_name}.cbc' )
Q_rch  =  cvv.get_data(text='RCH')[0]
Q_drn  =  cvv.get_data(text='DRN')[0]

#creates the 3d budget array
Q_rch_3d = cvv.get_data(idx=None, kstpkper=(0,0), totim=None, text='RCH',  full3D=True)[0].data
Q_drn_3d = cvv.get_data(idx=None, kstpkper=(0,0), totim=None, text='DRN',  full3D=True)[0].data


#flattens
Q_drn_2d = np.abs( np.sum(Q_drn_3d, axis=0) )
Q_rch_2d = np.abs( np.sum(Q_rch_3d, axis=0) )

Q_rch_2d = np.sum(Q_rch_3d, axis=0)


V_in = Q_rch['q'].sum()
V_out= Q_drn['q'].sum() 

print(f'V_in: {V_in} m^3/d')
print(f'V_out: {V_out} m^3/d')
print(f'Relative Error: { (V_in+V_out)/V_in * 100} %')



Q_difference =  Q_rch_2d  - Q_drn_2d
Q_rch_2d_NET =  Q_difference * (Q_difference>0)
Q_drn_2d_NET =  np.abs( Q_difference * (Q_difference<0) )
Q_rch_2d_NET[catchment==0] = np.nan
Q_drn_2d_NET[catchment==0] = np.nan


Q_drn_min= np.nanmin(Q_drn_2d_NET[Q_drn_2d_NET>0])
Q_drn_max= np.nanmax(Q_drn_2d_NET)

Q_rch_min= np.nanmin(Q_rch_2d_NET[Q_rch_2d_NET>0])
Q_rch_max= np.nanmax(Q_rch_2d_NET)


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
plt.title('Recharge discharge $[L]^3/[T]$')
plt.show()

