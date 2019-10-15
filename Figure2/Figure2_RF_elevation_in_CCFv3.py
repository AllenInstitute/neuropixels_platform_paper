# note: use 'get_data_from_warehouse' to generate the 'df' object

plt.figure(417811)
plt.clf()

areas = ('VISp','VISl','VISrl','VISal','VISpm','VISam')

plt.figure(11111)
plt.clf()

np.random.seed(19)

total_units = 0

for area_idx, area in enumerate(areas):
    
    selection = (df.structure_acronym == area) & \
                (df.on_screen_rf < 0.01) & \
                (df.dorsal_ventral_ccf_coordinate > 0) & \
                (df.area_rf < 1500) 
    
    azi = df[selection].azimuth_rf
    elev = df[selection].elevation_rf
    size = df[selection].area_rf / 40
    dv = -df[selection].dorsal_ventral_ccf_coordinate
    ml = df[selection].left_right_ccf_coordinate
    ap = df[selection].anterior_posterior_ccf_coordinate
    ori = df[selection].pref_ori_dg
    
    total_units += np.sum(selection)
    
    plt.scatter(-ml/1000,-ap/1000,c=elev,s=size,cmap='jet', alpha=0.4, edgecolors='none')#, vmin=1, vmax=7)
    plt.ylim([-11,-6])
    plt.xlim([-10.5,-6.5])
    
plt.title(total_units)

# %%

areas = ('LGd','LP')

plt.figure(11112)
plt.clf()

np.random.seed(27)

for area_idx, area in enumerate(areas):
    
    total_units = 0
    
    selection = (df.structure_acronym == area) & \
                (df.p_value_rf < 0.001) & \
                (df.dorsal_ventral_ccf_coordinate > 0) & \
                (df.area_rf < 800) 
    
    azi = df[selection].azimuth_rf.values
    elev = df[selection].elevation_rf.values
    size = df[selection].area_rf.values / 20
    dv = -df[selection].dorsal_ventral_ccf_coordinate.values
    ml = df[selection].left_right_ccf_coordinate.values
    ap = df[selection].anterior_posterior_ccf_coordinate.values
    ori = df[selection].pref_ori_dg.values
    
    order = np.random.permutation(len(ori))
    
    plt.subplot(2,2,area_idx+1)
    plt.scatter(ml[order]/1000,dv[order]/1000,c=elev[order],s=size[order],cmap='jet', alpha=0.4, edgecolors='none', vmin=-40, vmax=60)
    
    plt.xlim([6.75,8.5])
    plt.ylim([-4,-2.25])
    
    plt.subplot(2,2,area_idx+3)
    plt.scatter(ap[order]/1000,dv[order]/1000,c=elev[order],s=size[order],cmap='jet', alpha=0.4, edgecolors='none', vmin=-35, vmax=60)
    plt.title('AP/DV ')
    
    plt.xlim([6.5,9.0])
    plt.ylim([-4.5,-2])
    
    plt.title(np.sum(selection))
