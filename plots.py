import matplotlib.pyplot as plt
import numpy as np


# the following values correspond to those used in the simulations

R = [0.0031536,	0.006794225, 0.014637715, 0.031536,	0.067942252, 0.146377145, 0.31536, 0.679422524,	1.463771455, 2.0, 2.6, 3.1536]  # dimensional recharge r (m/year). Given Ks=1e-7  m/s it corresponds to dimensionless recharges from 10^-3 to 1. 
IMP_DEPTH =[10, 100, 1000]           # formation depth D_i (m)       
MEAN_Y = [-16.12]                    # ln(K_G). Log conductivity at the ground surface           
VAR_Y = [0]                          # variance of the log conductivity      
TOPOGRAPHY_FACTOR = [0.25, 1, 4]     # topographic stretch factor T_f (-)       
ALPHA = [0.0001, 0.001, 0.01]        # exponential decay rate of the hydraulic conductivity with depth (m-1)

Ks=1e-7          # K_g (m/s)
n=0.1            #porosity (-)
delta_space= 100 # horizontal dimension of the model cells (m) 

R_K_ratio = np.array(R)/Ks /3600 / 24 /365



#import model outputs from shelve.out
import shelve
filename='***/***/shelve.out'
my_shelf = shelve.open(filename)
for key in my_shelf:
 globals()[key]=my_shelf[key]
my_shelf.close()



representative_length =  (catchment_area*1000**2)**0.5   # (m)


#%%################################
###################################
## traveltime and flowpath lengths
###################################
###################################


inizio=0 
fine=0

L_k=len(IMP_DEPTH)
L_i=len(TOPOGRAPHY_FACTOR)
L_j=len(ALPHA)
L_R=len(R_K_ratio)

k=0   # choses for which IMP_DEPTH the plot is made 


# !!! traveltime array corresponds to the age of the groundwater component of the hydrological response

#prepara i dati del baseflow
traveltime_ARRAY_no9999_no0_anni=np.copy(traveltime_ARRAY)
traveltime_ARRAY_no9999_no0_anni[ np.logical_or( traveltime_ARRAY==-9999, traveltime_ARRAY==0) ]=np.nan
traveltime_ARRAY_no9999_no0_anni=traveltime_ARRAY_no9999_no0_anni /365   # goes from days to years

flowpath_lengths_ARRAY_no9999_no0=np.copy(flowpath_lengths_ARRAY)
flowpath_lengths_ARRAY_no9999_no0[ np.logical_or( flowpath_lengths_ARRAY==-9999, flowpath_lengths_ARRAY==0) ]=np.nan
flowpath_lengths_ARRAY_no9999_no0=flowpath_lengths_ARRAY_no9999_no0/1000  # goes from m to km




#initializes empyt arrays for the quartiles (first 3 columns) and for the mean and stdev (last two columns)
traveltime_ARRAY_statistics= np.zeros( (len(flowpath_lengths_ARRAY_no9999_no0[0]),5), dtype=np.float32)
flowpath_lengths_ARRAY_statistics= np.zeros( (len(flowpath_lengths_ARRAY[0]),5), dtype=np.float32)

#calcola i quatili traveltime  flowpathlength
traveltime_ARRAY_statistics[:,0]=np.nanpercentile(traveltime_ARRAY_no9999_no0_anni,25, axis=0)     #quartiles
traveltime_ARRAY_statistics[:,1]=np.nanpercentile(traveltime_ARRAY_no9999_no0_anni,50, axis=0)
traveltime_ARRAY_statistics[:,2]=np.nanpercentile(traveltime_ARRAY_no9999_no0_anni,75, axis=0)
traveltime_ARRAY_statistics[:,3]=np.nanmean(traveltime_ARRAY_no9999_no0_anni, axis=0)              #mean
traveltime_ARRAY_statistics[:,4]=np.sqrt(np.nanvar(traveltime_ARRAY_no9999_no0_anni, axis=0))      #stdev

flowpath_lengths_ARRAY_statistics[:,0]=np.nanpercentile(flowpath_lengths_ARRAY_no9999_no0,25, axis=0)     #quartiles
flowpath_lengths_ARRAY_statistics[:,1]=np.nanpercentile(flowpath_lengths_ARRAY_no9999_no0,50, axis=0)
flowpath_lengths_ARRAY_statistics[:,2]=np.nanpercentile(flowpath_lengths_ARRAY_no9999_no0,75, axis=0)
flowpath_lengths_ARRAY_statistics[:,3]=np.nanmean(flowpath_lengths_ARRAY_no9999_no0, axis=0)              #mean
flowpath_lengths_ARRAY_statistics[:,4]=np.sqrt(np.var(flowpath_lengths_ARRAY_no9999_no0, axis=0))         #stdev


#####################################################################
#%% plot 1 traveltimes (axis scale options : plot, semilogx, semilogy, loglog)
#####################################################################
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine=  k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,1])
  axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,0], traveltime_ARRAY_statistics[inizio:fine,2], alpha=0.3)
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)
  
i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='time (years)', title= title, xlim=(0.001, 1), ylim=(10, np.nanmax(traveltime_ARRAY_statistics[:,1])))
    i+=1
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
    
# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes.pdf')

#####################################################################
#%% plot 1 NORMALIZED ver1 traveltimes (axis scale options : plot, semilogx, semilogy, loglog)
#####################################################################
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine=  k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,1]  *   (Ks * 3600*24*365    /  representative_length )    )
  axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,0] *   (Ks * 3600*24*365    /  representative_length ), traveltime_ARRAY_statistics[inizio:fine,2] *   (Ks * 3600*24*365    /  representative_length ), alpha=0.3)
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3] *   (Ks * 3600*24*365    /  representative_length ))
  #axs[i, j].fill_between(R_K_ratio, (traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2) *   (Ks * 3600*24*365    /  representative_length ), (traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2) *   (Ks * 3600*24*365    /  representative_length ), alpha=0.3)
  
i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='T (-)', title= title, xlim=(0.001, 1), ylim=(0,10))  # np.nanmax(traveltime_ARRAY_statistics[:,1]*   (Ks * 3600*24*365    /  representative_length )  )))
    i+=1
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
    
# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes.pdf')




######################################################################
#%% plot 2 flowpathlengths
######################################################################
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine=  k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,1])
  axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,0], flowpath_lengths_ARRAY_statistics[inizio:fine,2], alpha=0.3)
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='length (km)', title= title, xlim=(0.001, 1), ylim=(0.5, np.nanmax(flowpath_lengths_ARRAY_statistics[:,2])))
    i+=1

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/flowpath_lenths.pdf')

######################################################################
#%% plot 2 NORMALIZED flowpathlengths 
######################################################################
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine=  k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,1]  /  (representative_length /1000)  )
  axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,0]/  (representative_length /1000), flowpath_lengths_ARRAY_statistics[inizio:fine,2]/  (representative_length /1000), alpha=0.3)
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='L (-)', title= title, xlim=(0.001, 1), ylim=(0,2)) # np.nanmax(flowpath_lengths_ARRAY_statistics[:,2])))
    i+=1

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/flowpath_lenths.pdf')





###################################################################################################
#%% plot 3 & 4 : same as plot 1 and 2 but with all IMPDEPTHS (without quartiles)
###################################################################################################
k1=0
k2=1
k3=2

#plot traveltimes
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio1:fine1,1],color='blue')
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio2:fine2,1],color='blue')
  axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio3:fine3,1],color='red')

  
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='time (years)', xlim=(0.001, 1), ylim=(0, np.nanmax(traveltime_ARRAY_statistics[:,1])))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes_all.pdf')



#%%plot flowpathlengths
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio1:fine1,1],color='blue')
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio2:fine2,1],color='green')
  axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio3:fine3,1],color='red')
  
  
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='length (km)', xlim=(0.001, 1), ylim=(0.1, np.nanmax(flowpath_lengths_ARRAY_statistics[:,1])))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/flowpath_lenths_all.pdf')



###################################################################################################
#%% plot 3 & 4 NORMALIZED: same as plot 1 and 2 but with all IMPDEPTHS (without quartiles)
###################################################################################################
k1=0
k2=1
k3=2

#plot traveltimes
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  #axs[i, j].semilox(...)
  axs[i, j].loglog(R_K_ratio, traveltime_ARRAY_statistics[inizio1:fine1,1]*   (Ks * 3600*24*365    /  representative_length )  ,color='blue')
  axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio1:fine1,0] *   (Ks * 3600*24*365    /  representative_length ), traveltime_ARRAY_statistics[inizio1:fine1,2] *   (Ks * 3600*24*365    /  representative_length ), alpha=0.3 ,color='blue')
  
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio2:fine2,1]*   (Ks * 3600*24*365    /  representative_length ),color='blue')
  
  axs[i, j].loglog(R_K_ratio, traveltime_ARRAY_statistics[inizio3:fine3,1]*   (Ks * 3600*24*365    /  representative_length )  ,color='red')
  axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio3:fine3,0] *   (Ks * 3600*24*365    /  representative_length ), traveltime_ARRAY_statistics[inizio3:fine3,2] *   (Ks * 3600*24*365    /  representative_length ), alpha=0.3 ,color='red')
  
  
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='T (-)', xlim=(0.001, 1), ylim=(0.01,10 ))  #np.nanmax(traveltime_ARRAY_statistics[:,1]*   (Ks * 3600*24*365    /  representative_length )   )))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['$D_i=10m$', '$D_i=1000m$'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes_all_loglog_quartiles.pdf')



#%%plot flowpathlengths
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  #axs[i, j].semilox(...)
  axs[i, j].loglog(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio1:fine1,1]/  (representative_length /1000)  ,color='blue')
  axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio1:fine1,0] /  (representative_length /1000), flowpath_lengths_ARRAY_statistics[inizio1:fine1,2] /  (representative_length /1000), alpha=0.3 ,color='blue')
  

  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio2:fine2,1]/  (representative_length /1000)  ,color='green')

  axs[i, j].loglog(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio3:fine3,1]/  (representative_length /1000)  ,color='red')
  axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio3:fine3,0] /  (representative_length /1000), flowpath_lengths_ARRAY_statistics[inizio3:fine3,2] /  (representative_length /1000), alpha=0.3 ,color='red')

  
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='L (-)', xlim=(0.001, 1), ylim=(0.01,10)) # np.nanmax(flowpath_lengths_ARRAY_statistics[:,1]/  (representative_length /1000)  )))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['$D_i=10m$', '$D_i=1000m$'], fontsize=7)
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/flowpath_lenths_all_loglog_quartiles.pdf')



#%%################################
###################################
## Q_sw Q_gw ratio
###################################
###################################

#normalized Q_gw and Q_sw
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, Q_sw_fluxes_normalized[0,inizio1:fine1],color='blue', linestyle='dashed')
  axs[i, j].semilogx(R_K_ratio, Q_gw_fluxes_normalized[0,inizio1:fine1],color='blue')
  #axs[i, j].semilogx(R_K_ratio, Q_sw_fluxes_normalized[0,inizio2:fine2],color='red')
  #axs[i, j].semilogx(R_K_ratio, Q_gw_fluxes_normalized[0,inizio2:fine2],color='magenta')
  axs[i, j].semilogx(R_K_ratio, Q_sw_fluxes_normalized[0,inizio3:fine3],color='red', linestyle='dashed')
  axs[i, j].semilogx(R_K_ratio, Q_gw_fluxes_normalized[0,inizio3:fine3],color='red')
  
  
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='Q/Qtot', xlim=(0.001, 1), ylim=(-0.1, 1.1))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
ax.legend(['Qdr (ID=0m)', 'Qbf (ID=0m)', 'Qdr (ID=1000m)', 'Qbf (ID=1000m)'], fontsize=7)    
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/Q_sw_Q_gw ratio.pdf')




#%% Q_gw  Q_sw ratio
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, Q_gw_Q_sw_ratio_particle[0,inizio1:fine1],color='blue')
  #axs[i, j].semilogx(R_K_ratio, Q_gw_Q_sw_ratio_particle[0,inizio2:fine2],color='green')
  axs[i, j].semilogx(R_K_ratio, Q_gw_Q_sw_ratio_particle[0,inizio3:fine3],color='red')
  
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='Q_bf/Q_dr', xlim=(0, 1), ylim=(0.1, 1000))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/GW_SW_fluxes_ratio.pdf')




#total Q  (è unico, stesso medesimo per tutti i casi)
fig, axs = plt.subplots(1,1)
axs = plt.loglog(R_K_ratio, Q_tot[0,inizio1:fine1],color='blue')
plt.xlim(0,1), plt.ylim(0,np.nanmax(Q_tot)), plt.xlabel('R/K (-)'), plt.ylabel('Q_tot (m^3/s)')

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/Q_tot.pdf')

#%%################################
###################################
## GW_volume  - not normalized 
###################################
###################################

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, (  GW_volume[0,inizio1:fine1]  - ( catchment_area * n * IMP_DEPTH[k1] /1000    ) ),color='blue')
  #axs[i, j].loglog(R_K_ratio, (  GW_volume[0,inizio2:fine2]  - ( catchment_area * n * IMP_DEPTH[k2] /1000  ) ) ,color='green')
  axs[i, j].semilogx(R_K_ratio, (  GW_volume[0,inizio3:fine3]  - ( catchment_area * n * IMP_DEPTH[k3] /1000    ) ),color='red')
  
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='GW_vol (km^3)', xlim=(0.001, 1), ylim=(1, np.nanmax(GW_volume)))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/GW_vol.pdf')


#%%################################
###################################
## GW_volume  - normalized
###################################
###################################

#calculates the volume of land above the lowest point of the dtm
dem_shifted=np.copy(demData_original)
dem_shifted[dem_shifted==-99999]=np.nan
dem_shifted=dem_shifted-np.nanmin(dem_shifted)
vol= np.nansum(catchment_area * n * np.nanmean(dem_shifted/1000))



fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].semilogx(R_K_ratio, (  GW_volume[0,inizio1:fine1]  - ( catchment_area * n * IMP_DEPTH[k1] /1000    ) ) / (vol * TOPOGRAPHY_FACTOR[i] ) ,color='blue')
  #axs[i, j].semilogx(R_K_ratio, (  GW_volume[0,inizio2:fine2]  - ( catchment_area * n * IMP_DEPTH[k2] /1000  ) ) / (vol * TOPOGRAPHY_FACTOR[i] ) ,color='green')
  axs[i, j].semilogx(R_K_ratio, (  GW_volume[0,inizio3:fine3]  - ( catchment_area * n * IMP_DEPTH[k3] /1000    ) ) / (vol * TOPOGRAPHY_FACTOR[i] ) ,color='red')
  
  #axs[i, j].semilogx(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, flowpath_lengths_ARRAY_statistics[inizio:fine,3]+flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, flowpath_lengths_ARRAY_statistics[inizio:fine,3]-flowpath_lengths_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='GW_vol (-)', xlim=(0.001, 1), ylim=(0., 1.1))  #, ylim=(np.nanmin(GW_volume), np.nanmax(GW_volume)))
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/GW_vol_normalized.pdf')





#%%################################
###################################
## Upwelling fluxes  - normalized 
###################################
###################################

k=2    # selects the impdepth level
kk=2   # selects the alpha value
kkk=[2,5,8] # selects the R/K ratio values
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), 3)

for i in range(len(TOPOGRAPHY_FACTOR)):
  for w in range(len(kkk)):
      
   daplottare=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  kk*(L_R) + kkk[w]
   
# =============================================================================
#    #gw outflow normalized by R_K_ratio
#    im= axs[i, w].imshow(drain_fluxes_minus_recharge_2D_ARRAY[daplottare,:,:]/ R_K_ratio[kkk[w]], interpolation='none',vmin=-0.05, vmax=-0,  cmap='viridis')
# =============================================================================
   
   #gw outflow normalized by R
   im= axs[i, w].imshow(drain_fluxes_minus_recharge_2D_ARRAY[daplottare,:,:]/ (R_K_ratio[kkk[w]] * (Ks*3600*24) ), interpolation='none',vmin=-10, vmax=-0,  cmap='viridis')
   


for ax in axs.flat:
    ax.set(xlabel='x (-)', ylabel='y (-)')
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    ax.set(xticklabels=[]), ax.set(yticklabels=[])
    ax.tick_params(left = False), ax.tick_params(bottom = False)
    
cbar = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.89)
#cbar.set_ticks([-0.05,-0.04,-0.03,-0.02,-0.01,-0])
cbar.set_ticks([-10,-8,-6,-4,-2,-0])
#cbar.set_label('Flux / (R/K)')
cbar.set_label('Flux / R (-)')

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/baseflow_flux.pdf')

#%%################################
###################################
## Direct runoff
###################################
###################################

k=2    # selects the impdepth level
kk=2   # selects the alpha value
kkk=[2,5,8] # selects the R/K ratio values
fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), 3)

for i in range(len(TOPOGRAPHY_FACTOR)):
  for w in range(len(kkk)):
      
   daplottare=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  kk*(L_R) + kkk[w]
   
# =============================================================================
#    #directrunoff normalized by R_K_ratio
#    im= axs[i, w].imshow(drain_fluxes_directrunoff_2D_ARRAY[daplottare,:,:] / R_K_ratio[kkk[w]], interpolation='none',vmin=-0.05, vmax=-0,  cmap='viridis')
# =============================================================================
   
   #directrunoff normalized by R
   im= axs[i, w].imshow(drain_fluxes_directrunoff_2D_ARRAY[daplottare,:,:] / (R_K_ratio[kkk[w]] * (Ks*3600*24) ) , interpolation='none',vmin=-1, vmax=-0,  cmap='viridis')
   

for ax in axs.flat:
    ax.set(xlabel='x (-)', ylabel='y (-)')
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    ax.set(xticklabels=[]), ax.set(yticklabels=[])
    ax.tick_params(left = False), ax.tick_params(bottom = False)
    
cbar = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.89)
#cbar.set_ticks([-0.05,-0.04,-0.03,-0.02,-0.01,-0])
cbar.set_ticks([-1,-0.8,-0.6,-0.4,-0.2,-0])
cbar.set_label('Flux / R (-)')

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/direct_runoff_flux.pdf')


#%%################################
###################################
## GW_ AGE moments   !!!!
###################################
###################################

#%% single distribution
#età delgroundwater come grafico a barre a partire dalle dimensioni in GW_age_hist_natural e GW_age_hist_log
scenario = 22
 # quale risultato plottare
plt.subplot(211)
plt.bar(GW_age_hist_natural[scenario][1], GW_age_hist_natural[scenario][0], align='edge', width=GW_age_hist_natural[scenario][2])
plt.subplot(212)
plt.bar(GW_age_hist_log[scenario][1], GW_age_hist_log[scenario][0], align='edge', width=GW_age_hist_log[scenario][2])
plt.xscale('log')
plt.xlim([1, 10e7])
plt.show()


#%% all



fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  
  
  
# =============================================================================
#   # median gw age
#   axs[i, j].loglog(R_K_ratio, GW_age_moments[0,inizio1:fine1]            ,color='blue')
#   axs[i, j].loglog(R_K_ratio, GW_age_moments[0,inizio3:fine3]            ,color='red')
#   for ax in axs.flat:
#     ax.set(xlabel='R/K (-)', ylabel='GW_age median', xlim=(0.001,1 ) , ylim=(np.nanmin(GW_age_moments[0,:]), np.nanmax(GW_age_moments[0,:])))
# =============================================================================
   
    
# =============================================================================
#   # cv
#   axs[i, j].semilogx(R_K_ratio, np.sqrt( GW_age_moments[1,inizio1:fine1] ) /  GW_age_moments[0,inizio1:fine1]            ,color='blue')
#   axs[i, j].semilogx(R_K_ratio, np.sqrt(GW_age_moments[1,inizio3:fine3])   /GW_age_moments[0,inizio3:fine3]            ,color='red')
#   for ax in axs.flat:
#     ax.set(xlabel='R/K (-)', ylabel='GW_age cv', xlim=(0.001,1 ) , ylim=(0.90, 1.2)  )
# =============================================================================
      
   
    
# =============================================================================
#   #gw age normalized ver 1
#   axs[i, j].semilogx(R_K_ratio, GW_age_moments[0,inizio1:fine1] *   (Ks * 3600*24*365    /  representative_length )            ,color='blue')
#   axs[i, j].semilogx(R_K_ratio, GW_age_moments[0,inizio3:fine3] *   (Ks * 3600*24*365    /  representative_length )            ,color='red')
#   for ax in axs.flat:
#     ax.set(xlabel='R/K (-)', ylabel='GW_age mean (-)', xlim=(0.001,1 ), ylim=(1e1,1e4))
# =============================================================================
   
    
# =============================================================================
#   #gw age normalized ver 2
#   axs[i, j].semilogx(R_K_ratio* (representative_length/(TOPOGRAPHY_FACTOR[i]*representative_hight)  )**2, GW_age_moments[0,inizio1:fine1] *   (Ks * 3600*24*365    /  representative_length )            ,color='blue')
#   axs[i, j].semilogx(R_K_ratio* (representative_length/(TOPOGRAPHY_FACTOR[i]*representative_hight)  )**2, GW_age_moments[0,inizio3:fine3] *   (Ks * 3600*24*365    /  representative_length )            ,color='red')
#   for ax in axs.flat:
#     ax.set(xlabel='(R/K)(L/H)$^2$ (-)', ylabel='GW_age mean (-)', xlim=(0.01,100 ), ylim=(1e1,1e4))    
# =============================================================================
   
    
  #ratio of the baseflow traveltimve and the groundwater age
  times_ratio =  traveltime_ARRAY_statistics[:,3]/  GW_age_moments[0,:]  
  times_ratio[np.isnan(times_ratio)] =  np.nan
  
  axs[i, j].semilogx(R_K_ratio,  times_ratio[inizio1:fine1] , color= 'blue')
  axs[i, j].semilogx(R_K_ratio,  times_ratio[inizio3:fine3] , color= 'red')
  for ax in axs.flat:
    ax.set(xlabel='R/K (-)', ylabel='$A_{bf}/A_{gw}$', xlim=(0.001,1 ),  ylim=(np.nanmin(times_ratio[times_ratio>0]),np.nanmax(times_ratio)) )
   



    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/GW_age_BF_age_means_ratio.pdf')




#%%################################
###################################
## STREAMFLOW age (considering also the direct runoff component with age 0)
###################################
###################################

#creates the array with streamflow age by addin zeros in a number proportional to the directrunoff flux to the traveltime array

age_streamflow =   np.copy(traveltime_ARRAY_no9999_no0_anni)

number_baseflow_particles =  np.sum(~np.isnan(age_streamflow[:,:]), axis=0)
number_directrunoff_particles =  (number_baseflow_particles   /   Q_gw_Q_sw_ratio_fluxes).astype(int)
max_directrunoff_particles= number_directrunoff_particles.max()

age_streamflow = np.vstack( (age_streamflow , np.zeros( ( max_directrunoff_particles, age_streamflow.shape[1]))*np.nan))

for i in range(age_streamflow.shape[1]):
    for j in range( number_directrunoff_particles[0,i] ):
        age_streamflow[ len(age_streamflow) - j -1, i] = 0
        
#calcola i quatili   
age_streamflow_statistics= np.zeros( (len(flowpath_lengths_ARRAY_no9999_no0[0]),5), dtype=np.float32)

age_streamflow_statistics[:,0]=np.nanpercentile(age_streamflow,25, axis=0)     #quartili
age_streamflow_statistics[:,1]=np.nanpercentile(age_streamflow,50, axis=0)
age_streamflow_statistics[:,2]=np.nanpercentile(age_streamflow,75, axis=0)
age_streamflow_statistics[:,3]=np.nanmean(age_streamflow, axis=0)              #media
age_streamflow_statistics[:,4]=np.sqrt(np.nanvar(age_streamflow, axis=0))      #deviazione statndard





#####################################################################
# comparison between the age of baseflow and the age of streamflow 
#####################################################################
k=2

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine=  k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #gw outflow age
  axs[i, j].loglog(R_K_ratio, age_streamflow_statistics[inizio:fine,3], color = 'blue')
  #axs[i, j].fill_between(R_K_ratio, age_streamflow_statistics[inizio:fine,0], age_streamflow_statistics[inizio:fine,2] , alpha=0.3)
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)
  
  axs[i, j].loglog(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3], color = 'red')
  
i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='time (years)', title= title, xlim=(0.001, 1), ylim=(10, 1e5))
    i+=1
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

ax.legend(['<TT> streamflow','<TT> baseflow'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/TT_baseflow_vs_streamflow.pdf')




#%%################################
###################################
## STREAMFLOW age - NORMALIZED ver1 (considerando anche l'età del direct runoff)
###################################
###################################


#####################################################################
# comparison between the age of baseflow and the age of streamflow 
#####################################################################
k=0

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio=k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine=  k*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #gw outflow age
  axs[i, j].loglog(R_K_ratio, age_streamflow_statistics[inizio:fine,3] *   (Ks * 3600*24*365    /  representative_length )   , color = 'blue')
  #axs[i, j].fill_between(R_K_ratio, age_streamflow_statistics[inizio:fine,0] *   (Ks * 3600*24*365    /  representative_length )   , age_streamflow_statistics[inizio:fine,2] , alpha=0.3)
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3] *   (Ks * 3600*24*365    /  representative_length )   )
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)
  
  axs[i, j].loglog(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3] *   (Ks * 3600*24*365    /  representative_length )  , color = 'red')
  
i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='T (-)', title= title, xlim=(0.001, 1), ylim=(0.01, 10))
    i+=1
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

ax.legend(['<TT> streamflow','<TT> baseflow'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/TT_baseflow_vs_streamflow.pdf')

#####################################################################
# ratio between the age of baseflow and the age of streamflow 
#####################################################################

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
plt.rc('axes', titlesize=10) #sets the fontsize of the title
TF_alpha=[]

for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  TF_alpha.append(fr'TF={TOPOGRAPHY_FACTOR[i]}, $\alpha$= {ALPHA[j]}')
  #ratio age streamflow tt and groundwater ouflows tt
  axs[i, j].semilogx(     R_K_ratio, (  age_streamflow_statistics[inizio1:fine1,3] *   (Ks * 3600*24*365    /  representative_length ) )   /  ( traveltime_ARRAY_statistics[inizio1:fine1,3] *   (Ks * 3600*24*365    /  representative_length ) )  , color = 'blue')
  axs[i, j].semilogx(     R_K_ratio, (  age_streamflow_statistics[inizio3:fine3,3] *   (Ks * 3600*24*365    /  representative_length ) )   /  ( traveltime_ARRAY_statistics[inizio3:fine3,3] *   (Ks * 3600*24*365    /  representative_length ) )  , color = 'red')

  
i=0
for ax in axs.flat:
    title= TF_alpha[i]
    ax.set(xlabel='R/K (-)', ylabel='$Ratio_{<TT>} (-)$', title= title, xlim=(0.001, 1), ylim=(0, 1.01))
    i+=1
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati
    
# set the spacing between subplots
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4,  hspace=0.4)

ax.legend(['imp_depth=0m', 'imp_depth=1000m'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/TT_baseflow_vs_streamflow_RATIO.pdf')


#%%##################################################################
# flowpath lenght vs residence time 
#####################################################################

k1=0
k2=1
k3=2

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].scatter(traveltime_ARRAY_no9999_no0_anni[:,inizio1+2], flowpath_lengths_ARRAY_no9999_no0[:,inizio1+2],color='blue',  marker = '.', s=0.001)
  axs[i, j].scatter(traveltime_ARRAY_no9999_no0_anni[:,fine1-1-3], flowpath_lengths_ARRAY_no9999_no0[:,fine1-1-3],color='red',  marker = '.', s=0.001)


  
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='time (years)', ylabel='distance (km)', xlim=(0.1, 1000000), ylim=(0.01, 100))
    ax.set_xscale("log")
    ax.set_yscale("log")
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 0.005', 'R/K = 0.5'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes_all.pdf')
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes_all.png', dpi=500, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None, metadata=None)


#%%##################################################################
# flowpath lenght vs residence time- NORMALIZED
#####################################################################

k1=0
k2=1
k3=2

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  #gw outflow age
  axs[i, j].scatter(traveltime_ARRAY_no9999_no0_anni[:,inizio1+2]  * (Ks * 3600*24*365    /  representative_length )   , flowpath_lengths_ARRAY_no9999_no0[:,inizio1+2] / (representative_length/1000)  ,color='blue',  marker = '.', s=0.001)
  axs[i, j].scatter(traveltime_ARRAY_no9999_no0_anni[:,fine1-1-3]  * (Ks * 3600*24*365    /  representative_length )   , flowpath_lengths_ARRAY_no9999_no0[:,fine1-1-3] / (representative_length/1000)  ,color='red',  marker = '.', s=0.001)


  
  #axs[i, j].semilogx(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3])
  #axs[i, j].fill_between(R_K_ratio, traveltime_ARRAY_statistics[inizio:fine,3]+traveltime_ARRAY_statistics[inizio:fine,4]/2, traveltime_ARRAY_statistics[inizio:fine,3]-traveltime_ARRAY_statistics[inizio:fine,4]/2, alpha=0.3)

for ax in axs.flat:
    ax.set(xlabel='T (-)', ylabel='distance (-)', xlim=(0.0001, 100), ylim=(0.01, 10))
    ax.set_xscale("log")
    ax.set_yscale("log")
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 0.005', 'R/K = 0.5'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes_all.pdf')
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/traveltimes_all.png', dpi=500, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None, metadata=None)




#%%density plot - single plot :
fig, ax = plt.subplots()
x_min=-1 ; x_max = 6 ; y_min=-2 ; y_max = 2  # in log10

ax.hexbin((traveltime_ARRAY_no9999_no0_anni[:,inizio3+2]), (flowpath_lengths_ARRAY_no9999_no0[:,inizio3+2]), gridsize=100, xscale='log', yscale='log', bins='log' ,marginals=False, cmap= 'viridis', extent=(x_min, x_max, y_min, y_max))

ax.set(xlim=(10**x_min, 10**x_max), ylim=(10**y_min, 10**y_max))

plt.show()


#%% CDF average flow velocity

k1=0
k2=1
k3=2

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     
     
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  

  
  velocity_1 = flowpath_lengths_ARRAY_no9999_no0[:,inizio1+2]  /  traveltime_ARRAY_no9999_no0_anni[:,inizio1+2]
  velocity_11 = flowpath_lengths_ARRAY_no9999_no0[:,fine1-1-3]  /  traveltime_ARRAY_no9999_no0_anni[:,fine1-1-3]
  velocity_2 = flowpath_lengths_ARRAY_no9999_no0[:,inizio3+2]  /  traveltime_ARRAY_no9999_no0_anni[:,inizio3+2]
  velocity_22 = flowpath_lengths_ARRAY_no9999_no0[:,fine3-1-3]  /  traveltime_ARRAY_no9999_no0_anni[:,fine3-1-3]
  
  binmin=0.0001
  binmax=0.1 # np.nanmax(velocity_1)
  numbins= 200
  bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins   )
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
  

  
  #axs[i, j].hist( velocity , bins=25, range= (0, np.nanmax(velocity)  ))   # linscale
  axs[i, j].hist( velocity_1 , bins=logbins, density=True, cumulative=True, histtype='step', color='blue')    #logscale
  axs[i, j].hist( velocity_11 , bins=logbins, density=True, cumulative=True, histtype='step', color='blue',linestyle='dashed')    #logscale

  axs[i, j].hist( velocity_2 , bins=logbins, density=True, cumulative=True, histtype='step', color='red')    #logscale
  axs[i, j].hist( velocity_22 , bins=logbins, density=True, cumulative=True, histtype='step', color='red',linestyle='dashed')    #logscale


for ax in axs.flat:
    ax.set(xlabel='velocity (km/years)', ylabel='CDF (-)', xlim=(binmin, binmax), ylim=(0.001, 1))
    ax.set_xscale("log")
    ax.set_yscale("log")
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 0.005, ID=0', 'R/K = 0.5, ID=0', 'R/K = 0.005, ID=1000', 'R/K = 0.5, ID=1000'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/flow_velocities.pdf')



#%% CDF average flow velocity  -   NORMALIZED

k1=0
k2=1
k3=2

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     
     
  inizio1=k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine1=  k1*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio2=k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine2=  k2*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  
  inizio3=k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R)    # definisce i ragne di cui plottare i valori
  fine3=  k3*(L_i*L_j*L_R) +  i*(L_j*L_R) +  j*(L_R) + L_R
  

  
  velocity_1  = flowpath_lengths_ARRAY_no9999_no0[:,inizio1+2]  /  traveltime_ARRAY_no9999_no0_anni[:,inizio1+2]  / ( Ks * 3600 * 24 * 365 / 1000)  # km/year
  velocity_11 = flowpath_lengths_ARRAY_no9999_no0[:,fine1-1-3]  /  traveltime_ARRAY_no9999_no0_anni[:,fine1-1-3]  / ( Ks * 3600 * 24 * 365 / 1000)
  velocity_2  = flowpath_lengths_ARRAY_no9999_no0[:,inizio3+2]  /  traveltime_ARRAY_no9999_no0_anni[:,inizio3+2]  / ( Ks * 3600 * 24 * 365 / 1000)
  velocity_22 = flowpath_lengths_ARRAY_no9999_no0[:,fine3-1-3]  /  traveltime_ARRAY_no9999_no0_anni[:,fine3-1-3]  / ( Ks * 3600 * 24 * 365 / 1000)
  
  binmin=0.01
  binmax=10 # np.nanmax(velocity_1)
  numbins= 200
  bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins   )
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
  

  
  #axs[i, j].hist( velocity , bins=25, range= (0, np.nanmax(velocity)  ))   # linscale
  axs[i, j].hist( velocity_1 , bins=logbins, density=True, cumulative=True, histtype='step', color='blue')    #logscale
  axs[i, j].hist( velocity_11 , bins=logbins, density=True, cumulative=True, histtype='step', color='blue',linestyle='dashed')    #logscale

  axs[i, j].hist( velocity_2 , bins=logbins, density=True, cumulative=True, histtype='step', color='red')    #logscale
  axs[i, j].hist( velocity_22 , bins=logbins, density=True, cumulative=True, histtype='step', color='red',linestyle='dashed')    #logscale


for ax in axs.flat:
    ax.set(xlabel='velocity (-)', ylabel='CDF (-)', xlim=(binmin, binmax), ylim=(0.001, 1))
    ax.set_xscale("log")
    ax.set_yscale("log")
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 0.005, ID=0', 'R/K = 0.5, ID=0', 'R/K = 0.005, ID=1000', 'R/K = 0.5, ID=1000'], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/flow_velocities.pdf')




#######################
#%% plot width function 
#######################


import rasterio

dem      = rasterio.open ( '***/***/DEM_Maso_100m.tif')
dem_data = dem.read()


src  = rasterio.open('E:/Trento/Paper_1/analisi_vario/gis/distance_to_outlet.tif')
data = src.read()
distance_linear = data[data>0]

plt.hist(distance_linear/1000,density=True, range=(0,19.5), bins=50, align='mid')
plt.xlabel('Distance (km)')
plt.ylabel('Frequency')
plt.title('Width Function')

src.close()
dem.close()
#plt.savefig('width_function.pdf', format='pdf')



##############
#%% plot selected traveltimes
##############

ID    = 2 
TF    = 1
alpha = 0

pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 

rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
selected_scenarios= rapporto + pointer

binmin=0.0001
binmax=1000 # np.nanmax(velocity_1)
numbins= 40
bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins   )
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))


plt0 = plt.hist( traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[0]] * (Ks * 3600*24*365    /  representative_length ), bins=logbins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   
plt1 = plt.hist( traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[3]] * (Ks * 3600*24*365    /  representative_length ), bins=logbins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   
plt2 = plt.hist( traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[6]] * (Ks * 3600*24*365    /  representative_length ), bins=logbins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   
plt3 = plt.hist( traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[9]] * (Ks * 3600*24*365    /  representative_length ), bins=logbins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   


plt.scatter(plt3[1][0:numbins-1],plt3[0], c='blue') # plotting t, b separately 
plt.scatter(plt2[1][0:numbins-1],plt2[0], c='lime') # plotting t, b separately 
plt.scatter(plt1[1][0:numbins-1],plt1[0], c='orange') # plotting t, b separately 
plt.scatter(plt0[1][0:numbins-1],plt0[0] , c='red') # plotting t, a separately 
plt.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001',], fontsize=7)
plt.xlabel('traveltimes(-)') #, xlim=(0.1, binmax))#, ylim=(0.000001, 0.25))
plt.ylabel('PDF (-)') #, xlim=(0.1, binmax))#, ylim=(0.000001, 0.25))
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.00001, binmax)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/aaa_prove_distribuzioni_ALLINONE_traveltimes.pdf')

plt.show()


#%% all traveltimes with gamma fit (to normalize multiply by (Ks * 3600*24*365    / (n* representative_length )  
import numpy as np
from scipy.stats import gamma


GAMMA_PARAMS = []

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     

  ID    = 2 
  TF    = i  
  alpha = j

  pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 
  rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
  selected_scenarios= rapporto + pointer
     
    
     
  binmin=0.0001
  binmax=10000 # np.nanmax(velocity_1)
  numbins= 40
  bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins   )
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))


  dati0 = traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[0]]*(Ks * 3600*24*365    /  (n * representative_length) )
  dati1 = traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[3]]*(Ks * 3600*24*365    /  (n * representative_length) )
  dati2 = traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[6]]*(Ks * 3600*24*365    /  (n * representative_length) )
  dati3 = traveltime_ARRAY_no9999_no0_anni[:,selected_scenarios[9]]*(Ks * 3600*24*365    /  (n * representative_length) )

  #removes nan  
  dati0 = dati0[~np.isnan(dati0)]
  dati1 = dati1[~np.isnan(dati1)]
  dati2 = dati2[~np.isnan(dati2)]
  dati3 = dati3[~np.isnan(dati3)]
  

  #fit the gamma 
  shape2, loc2, scale2 = gamma.fit(dati2, floc =  0)
  
  shape0, loc0, scale0 = gamma.fit(dati0, floc= loc2 , scale = scale2 )
  shape1, loc1, scale1 = gamma.fit(dati1, floc= loc2 , scale = scale2 )
  shape3, loc3, scale3 = gamma.fit(dati3, floc= loc2 , scale = scale2 )

  GAMMA_PARAMS.append([[shape0, scale0], [shape1, scale1], [shape2, scale2], [shape3, scale3] ]) 


  #generate numbers from the fitted distribution
  random_numbers0 = gamma.rvs(shape0, loc=loc0, scale=scale0, size=1000000)
  random_numbers1 = gamma.rvs(shape1, loc=loc1, scale=scale1, size=1000000)
  random_numbers2 = gamma.rvs(shape2, loc=loc2, scale=scale2, size=1000000)
  random_numbers3 = gamma.rvs(shape3, loc=loc3, scale=scale3, size=1000000)
  
  
# crea gli istogrammi e poi fa il plot  
  plt0 = plt.hist(dati0 , bins=logbins, density=True, cumulative=False, histtype='step', color='blue'  , alpha = 0)    #logscale
  plt1 = plt.hist(dati1 , bins=logbins, density=True, cumulative=False, histtype='step', color='lime'  , alpha = 0)    #logscale
  plt2 = plt.hist(dati2 , bins=logbins, density=True, cumulative=False, histtype='step', color='orange', alpha = 0)    #logscale
  plt3 = plt.hist(dati3 , bins=logbins, density=True, cumulative=False, histtype='step', color='red'   , alpha = 0)    #logscale

  fit0 = plt.hist(random_numbers0, bins=logbins, density=True, cumulative=False, histtype='step', color='blue'  , alpha = 0)
  fit1 = plt.hist(random_numbers1, bins=logbins, density=True, cumulative=False, histtype='step', color='lime'  , alpha = 0)  
  fit2 = plt.hist(random_numbers2, bins=logbins, density=True, cumulative=False, histtype='step', color='orange', alpha = 0)
  fit3 = plt.hist(random_numbers3, bins=logbins, density=True, cumulative=False, histtype='step', color='red'   , alpha = 0)

  
  axs[i, j].scatter(plt3[1][0:numbins-1],plt3[0], color = 'blue'  , s=5,  facecolors='none') # togliere "s=5" se si plotta con plot invece che con scatter
  axs[i, j].scatter(plt2[1][0:numbins-1],plt2[0], color = 'lime'  , s=5,  facecolors='none') 
  axs[i, j].scatter(plt1[1][0:numbins-1],plt1[0], color = 'orange', s=5,  facecolors='none') 
  axs[i, j].scatter(plt0[1][0:numbins-1],plt0[0], color = 'red'   , s=5,  facecolors='none') 

  axs[i, j].plot(fit3[1][0:numbins-1],fit3[0], color = 'blue'  ,  linewidth=1  ) # togliere "s=5" se si plotta con plot invece che con scatter
  axs[i, j].plot(fit2[1][0:numbins-1],fit2[0], color = 'lime'  ,  linewidth=1  ) 
  axs[i, j].plot(fit1[1][0:numbins-1],fit1[0], color = 'orange',  linewidth=1) 
  axs[i, j].plot(fit0[1][0:numbins-1],fit0[0], color = 'red'   ,  linewidth=1  ) 



for ax in axs.flat:
    #ax.set(xlabel='T (-)', ylabel='PDF (-)', xlim=(binmin, binmax), ylim=(0.000001, 1000))
    ax.set(xlabel='T (-)', ylabel='PDF (-)', xlim=(0.0001, 10000), ylim=(0.000001, 1000))
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.set_xticks([binmin,binmax])
    ax.set_xticks([0.0001,0.01,1,100,10000])
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001'], fontsize=5)

GAMMA_PARAMS_array=np.array(GAMMA_PARAMS)
GAMMA_PARAMS_array_reshaped= np.reshape(GAMMA_PARAMS_array, (GAMMA_PARAMS_array.shape[0]*GAMMA_PARAMS_array.shape[1] ,GAMMA_PARAMS_array.shape[2] ))
 #plt.savefig('E:/Post_TN/Risultati/plots_python/plots/distribuzioni_ALL_traveltimes_FIT_Di_1000m.pdf')



#################
#%% flowpath length
#################

ID    = 2 
TF    = 1
alpha = 0

pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 

rapporto=np.array([0,1,2,3,4,5,6,7,8,10])
selected_scenarios= rapporto + pointer

binmin=0.0001
binmax=30 # np.nanmax(velocity_1)
numbins= 50
bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins   )
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))


fig, axs = plt.subplots(len(selected_scenarios), 2,  figsize=(5, 20) )


for i in range(len(selected_scenarios)):

       
  #axs[i, j].hist( velocity , bins=25, range= (0, np.nanmax(velocity)  ))   # linscale

  axs[i,0].hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[i]] , bins=bins, density=True, cumulative=False, histtype='step', color='blue')    #logscale
  
  axs[i,1].hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[i]] , bins=bins, density=True, cumulative=True, histtype='step', color='blue')    #logscale


  
for ax in axs.flat:
    ax.set(xlabel='flowpath length (-)', ylabel=' (-)', xlim=(0.1, binmax))#, ylim=(0.000001, 0.25))
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend([''], fontsize=7)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/aaa_prove_distribuzioni_flowpathlengths.pdf')


#%% as above but with multiple plots
ID    = 0 
TF    = 2
alpha = 2

pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 
rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
selected_scenarios= rapporto + pointer

plt0 = plt.hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[0]] , bins=bins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   
plt1 = plt.hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[3]] , bins=bins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   
plt2 = plt.hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[6]] , bins=bins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   
plt3 = plt.hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[9]] , bins=bins, density=True, cumulative=False, histtype='step', color='blue', alpha = 0)   


plt.plot(plt3[1][0:numbins-1],plt3[0], 'blue') # plotting t, b separately 
plt.plot(plt2[1][0:numbins-1],plt2[0], 'lime') # plotting t, b separately 
plt.plot(plt1[1][0:numbins-1],plt1[0], 'orange') # plotting t, b separately 
plt.plot(plt0[1][0:numbins-1],plt0[0] , 'red') # plotting t, a separately 
plt.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001'], fontsize=7)
plt.xlabel('flowpath length (km)') #, xlim=(0.1, binmax))#, ylim=(0.000001, 0.25))
plt.ylabel('probability (-)') #, xlim=(0.1, binmax))#, ylim=(0.000001, 0.25))

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/aaa_prove_distribuzioni_ALLINONE_flowpathlengths.pdf')

plt.show()


#%% flowpathlength tutte insieme   (to normalize DIVIDE by  (representative_length/1000)   

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     

  ID    = 2 
  TF    = i  
  alpha = j

  pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 
  rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
  selected_scenarios= rapporto + pointer
     
    
     
  binmin=0.001
  binmax=3 # np.nanmax(velocity_1)
  numbins= 40
  bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins   )
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))


  
# alternatively plots directly with hist
  axs[i, j].hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[9]]/ (representative_length/1000) , bins=bins, density=True, cumulative=False, histtype='step', color='blue' )    #logscale
  axs[i, j].hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[6]]/ (representative_length/1000) , bins=bins, density=True, cumulative=False, histtype='step', color='lime'  )    #logscale
  axs[i, j].hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[3]]/ (representative_length/1000) , bins=bins, density=True, cumulative=False, histtype='step', color='orange')    #logscale
  axs[i, j].hist( flowpath_lengths_ARRAY_no9999_no0[:,selected_scenarios[0]]/ (representative_length/1000) , bins=bins, density=True, cumulative=False, histtype='step', color='red'   )    #logscale




for ax in axs.flat:
    ax.set(xlabel='Flowpath lengths (-)', ylabel='PDF (-)', xlim=(binmin, binmax), ylim=(0, 6))
    #ax.set_xscale("log")
    #☻ax.set_yscale("log")
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001'], fontsize=5)


#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/aaa_prove_distribuzioni_ALL_flowpathlenghts.pdf')




###################################################
#%% distribution of normalized fluxes  - fluxes in drain_fluxes_minus_recharge_2D_ARRAY are in m/year
####################################################


#linearizes otward fluxes 
linearized_gw_outflows = np.reshape (drain_fluxes_minus_recharge_2D_ARRAY[0,:,:], ( drain_fluxes_minus_recharge_2D_ARRAY[0,:,:].size, 1 )  )

for i in range(1,np.shape(drain_fluxes_minus_recharge_2D_ARRAY)[0]):
    
   linearized_gw_outflows=  np.append( linearized_gw_outflows, np.reshape (drain_fluxes_minus_recharge_2D_ARRAY[i,:,:], ( drain_fluxes_minus_recharge_2D_ARRAY[0,:,:].size, 1 )  ), axis = 1 )


#positive values are recharge (deleted)
linearized_gw_outflows[linearized_gw_outflows>0] = np.nan

linearized_gw_outflows=np.abs(linearized_gw_outflows)







binmin   = np.nanpercentile(linearized_gw_outflows,0.01) /(Ks*3600*24*365) # np.nanmin(linearized_gw_outflows) # np.nanpercentile(linearized_gw_outflows,1)
binmax   = np.nanpercentile(linearized_gw_outflows,99.99)/(Ks*3600*24*365)  # np.nanmax(linearized_gw_outflows) # np.nanpercentile(linearized_gw_outflows,99)
numbins  = 50

bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))




fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     

  ID    = 2 
  TF    = i  
  alpha = j

  pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 
  rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
  selected_scenarios= rapporto + pointer
  


# crea gli istogrammi e poi fa il plot  
  plt0 = plt.hist( linearized_gw_outflows[:,selected_scenarios[0]]/(Ks*3600*24*365) , bins=logbins, density=True, cumulative=False, histtype='step', color='blue'  , alpha = 0)    #logscale
  plt1 = plt.hist( linearized_gw_outflows[:,selected_scenarios[3]]/(Ks*3600*24*365) , bins=logbins, density=True, cumulative=False, histtype='step', color='lime'  , alpha = 0)    #logscale
  plt2 = plt.hist( linearized_gw_outflows[:,selected_scenarios[6]]/(Ks*3600*24*365) , bins=logbins, density=True, cumulative=False, histtype='step', color='orange', alpha = 0)    #logscale
  plt3 = plt.hist( linearized_gw_outflows[:,selected_scenarios[9]]/(Ks*3600*24*365) , bins=logbins, density=True, cumulative=False, histtype='step', color='red'   , alpha = 0)    #logscale

  axs[i, j].scatter(plt3[1][0:numbins-1],plt3[0], color = 'blue'  , s=5) # togliere "s=5" se si plotta con plot invece che con scatter
  axs[i, j].scatter(plt2[1][0:numbins-1],plt2[0], color = 'lime'  , s=5) 
  axs[i, j].scatter(plt1[1][0:numbins-1],plt1[0], color = 'orange', s=5) 
  axs[i, j].scatter(plt0[1][0:numbins-1],plt0[0], color = 'red'   , s=5) 



for ax in axs.flat:
    ax.set(xlabel='GW outflow (-)', ylabel='PDF (-)', xlim=(binmin, binmax), ylim=(0.01, 100000))
    ax.set_xscale("log")
    ax.set_yscale("log")
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

ax.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001'], fontsize=5)


#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/Distribuzione_GW_outflows.pdf')







###################################################
#%%## drainage density as a function of Q 
###################################################

def func_powerlaw(x, m, c):
    return  c*x**m 


binmin   = 0.1   # np.nanmin(linearized_gw_outflows) # np.nanpercentile(linearized_gw_outflows,1)
binmax   = 1000  # np.nanmax(linearized_gw_outflows) # np.nanpercentile(linearized_gw_outflows,99)
numbins  = 50

bins= np.arange( start = binmin, stop = binmax, step = (binmax-binmin)/numbins)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))



fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     

  ID    = 2 
  TF    = i  
  alpha = j

  pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 
  rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
  selected_scenarios= rapporto + pointer
  
  stream_length = np.sum( linearized_gw_outflows[:,selected_scenarios]>0, axis = 0)
  
  #stream_length = stream_length * 100 / 1000 *0.47  # goes from pxls to km (0.47*dx average length of a segment inside a 1x1 pixel of side dx)
  
  
  drainage_density = stream_length *(0.1**2) /catchment_area
  
  axs[i, j].scatter(R_K_ratio , drainage_density ,   color = 'blue'  , s=10)
  axs[i, j].grid(color = 'gray', linestyle = '--', linewidth = 0.5)

  if i == 0 or (i ==1 and j==2):
    popt, pcov = curve_fit(func_powerlaw, R_K_ratio[0:6], drainage_density[0:6] , p0 = (0,100))
  else:
    popt, pcov = curve_fit(func_powerlaw, R_K_ratio, drainage_density , p0 = (0,100))
  
  axs[i, j].plot(R_K_ratio, func_powerlaw(R_K_ratio, *popt), color="red", linewidth=1, linestyle="--")
  axs[i, j].text(0.005, 0.0002, f'$y(x)={round(popt[1],2)}x**{round(popt[0],2)}$', fontsize=7)


for ax in axs.flat:
    
    ax.set(xlabel='$Q_{sf}/(K_G A)$ (-)', ylabel='$D_d$ (-)', xlim=(np.min(R_K_ratio), np.max(R_K_ratio)) , ylim=(0.0001, 1))
    ax.set_xscale("log")
    ax.set_yscale("log")

    
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

#ax.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001'], fontsize=5)



GW_outflow_ALL = np.nansum( linearized_gw_outflows, axis = 0)
GW_outflow_ALL_reshaped = np.reshape ( GW_outflow_ALL, (27,12))

stream_length_ALL = np.nansum( linearized_gw_outflows>0  * 100 / 1000 *0.47, axis = 0)
stream_length_ALL_reshaped = np.reshape ( stream_length_ALL, (27,12))

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/drainage_density.pdf')


##########################################################################################################

a = [0.01, 0.1, 1]

Q_adim = np.logspace(-3, 1, 1000)


Dd     = a[i]* Q_adim**(2/3)    

plt.plot(Q_adim, a[0]* Q_adim**(2/3), 'r') # plotting t, a separately 
plt.plot(Q_adim, a[1]* Q_adim**(2/3), 'b') # plotting t, b separately 
plt.plot(Q_adim, a[2]* Q_adim**(2/3), 'g') # plotting t, c separately 
plt.xscale("log")
plt.yscale("log")
ax = plt.gca()
ax.set(xlabel='$Q_{sf}/(K_G A)$ (-)', ylabel='$D_d$ (-)', xlim=(np.min(R_K_ratio), np.max(R_K_ratio)) , ylim=(0.0001, 1))
plt.show()

plt.plot(Q_adim, 2/3*a[0]* Q_adim**(2/3-1), 'r') # plotting t, a separately 
plt.plot(Q_adim, 2/3*a[1]* Q_adim**(2/3-1), 'b') # plotting t, b separately 
plt.plot(Q_adim, 2/3*a[2]* Q_adim**(2/3-1), 'g') # plotting t, c separately 
plt.xscale("log")
plt.yscale("log")
ax = plt.gca()
ax.set_xlim([0.001, 1])
ax.set_ylim([0.01, 100])
plt.show()


###### alternative visualization - drainage density vs dimensionless Q 

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np


gam = 0.7

#A= np.arange(0.1,10,0.1)
A= np.logspace(-2, 1, 20)

# setup the normalization and the colormap
normalize = mcolors.LogNorm(vmin=A.min(), vmax=A.max())
colormap = cm.jet

# plot
for a in A:
    plt.plot(Q_adim,  a* Q_adim**gam , color=colormap(normalize(a)) , linewidth=2)

# setup the colorbar
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(A)
plt.colorbar(scalarmappaple, label= '$\lambda$ (-)')

ax = plt.gca()
ax.set(xlabel='$Q_{sf}/(K_G A)$ (-)', ylabel='$D_d$ (-)', xlim=(np.min(R_K_ratio), np.max(R_K_ratio)) , ylim=(0.0001, 1))
plt.xscale("log")
plt.yscale("log")
# show the figure
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/drainage_density_colorscale.pdf')

plt.show()





######  alternative visualization - derivative of the drainage density vs dimensionless Q 

# setup the normalization and the colormap
normalize = mcolors.LogNorm(vmin=A.min(), vmax=A.max())
colormap = cm.jet

# plot
for a in A:
    plt.plot(Q_adim,  gam*a* Q_adim**(gam-1) , color=colormap(normalize(a)) , linewidth=2)

# setup the colorbar
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(A)
plt.colorbar(scalarmappaple, label= '$\lambda$ (-)')

ax = plt.gca()
ax.set(xlabel='$Q_{sf}/(K_G A)$ (-)', ylabel='$D^I_d$ (-)', xlim=(np.min(R_K_ratio), np.max(R_K_ratio)) , ylim=(0.01, 100))
plt.xscale("log")
plt.yscale("log")
#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/drainage_density_derivative_colorscale.pdf')
# show the figure
plt.show()






##################################################################################
### as above but instead of using Q_sf uses Q_bf (baseflow rather than streamflow)
###################################################################################

import matplotlib.colors as colors

fig, axs = plt.subplots(len(TOPOGRAPHY_FACTOR), len(ALPHA))
for i in range(len(TOPOGRAPHY_FACTOR)):
 for j in range(len(ALPHA)):
     

  ID    = 2 
  TF    = i  
  alpha = j

  pointer = ID * (12 * 3 * 3) + TF * ( 12 * 3 ) + alpha * (12) 
  rapporto=np.array([0,1,2,3,4,5,6,7,8,9,10,11])
  selected_scenarios= rapporto + pointer
  
  stream_length = np.sum( linearized_gw_outflows[:,selected_scenarios]>0, axis = 0)
  
  #stream_length = stream_length * 100 / 1000 *0.47  # goes from pxls to km (0.47*dx average length of a segment inside a 1x1 pixel of side dx)
  
  
  drainage_density = stream_length *(0.1**2) /catchment_area
  
  
  Q_gw_fluxes_normalized_select = Q_gw_fluxes_normalized[:,selected_scenarios]



  
  axs[i, j].scatter(Q_gw_fluxes_normalized_select , drainage_density , c = R_K_ratio  , s=10, cmap = 'jet', norm = colors.LogNorm()) #    color = 'blue'  , s=10)
  axs[i, j].grid(color = 'gray', linestyle = '--', linewidth = 0.5)

for ax in axs.flat:
    
    ax.set(xlabel='$Q_{bf}/(R_{pot} A)$ (-)', ylabel='$D_d$ (-)', xlim=(0.0, 1) , ylim=(0.000, 1))
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    
    
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    
    #ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')   # fa i plot quadrati

#ax.legend(['R/K = 1', 'R/K = 0.1', 'R/K = 0.01', 'R/K = 0.001'], fontsize=5)

#plt.savefig('E:/Post_TN/Risultati/plots_python/plots/drainage_baseflow.pdf')



























