# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

df = pd.read_csv(os.path.join(os.getcwd(), 'data', 'unit_table.csv'), low_memory=False)

# %%

areas = ('VISp','VISl','VISrl','VISal','VISpm','VISam','LGd','LP') 
num_units_per_area = np.zeros((len(df.specimen_id.unique()), len(areas)))
num_overlap_per_area = np.zeros((len(df.specimen_id.unique()), len(areas), len(areas)))

min_units = 0

num_per_mouse = []

for mouse_idx, mouse in enumerate(df.specimen_id.unique()):
    
    selection1 =  (df.specimen_id == mouse) & \
                  (df.amplitude_cutoff < 0.1) & \
                  (df.presence_ratio > 0.95) & \
                  (df.isi_violations < 0.5) & \
                  (df.quality == 'good')
                  
    num_per_mouse.append(np.sum(selection1))
    
    for area_idx, area in enumerate(areas):
        
        selection2 = (df.ecephys_structure_acronym == area)
        
        if np.sum(selection1 * selection2) > min_units:
                    
            num_units_per_area[mouse_idx, area_idx] = np.sum(selection1 & selection2)
        
simultaneous_areas = np.sum(num_units_per_area[:,:].astype('bool'),1)
did_record = np.sum(num_units_per_area.astype('bool'),0)

# %%

num_per_probe = []
has_ccf = 0

for probe_idx, probe in enumerate(df.ecephys_probe_id.unique()):
    
    selection = (df.ecephys_probe_id == probe) & \
                    (df.amplitude_cutoff < 0.1) & \
                    (df.presence_ratio > 0.95) & \
                    (df.isi_violations < 0.5)
    
    num_per_probe.append(np.sum(selection))
    
    if np.sum(df[selection]['anterior_posterior_ccf_coordinate']) > 0:
        
        has_ccf += 1

print('Total units per probe:')
print('Mean: ' + str(np.mean(num_per_probe)) + '+/-' + str(np.std(num_per_probe)))
print()

print('Total units per mouse:')
print('Mean: ' + str(np.mean(num_per_mouse)) + '+/-' + str(np.std(num_per_mouse)))
print()

# %%

unit_count = num_units_per_area.flatten()
unit_count[unit_count == 0] = np.nan
unit_mean = np.nanmean(unit_count)
unit_std = np.nanstd(unit_count)

print('Total units per area:')
print('Mean: ' + str(unit_mean) + '+/-' + str(unit_std))
print()
print('Total areas per recording:')
print('Mean: ' + str(np.mean(simultaneous_areas)) + '+/-' + str(np.std(simultaneous_areas)))

# %%
plt.figure(4171)
plt.clf()

plt.subplot2grid((1,6),(0,4),colspan=2,rowspan=1)

common_area_names = ['V1','LM','RL','AL','PM','AM','LGd','LP'] #,'CA1', 'CA3', 'DG', 'SUB']

h,b = np.histogram(simultaneous_areas, bins=np.arange(2,11)) #, color='gray')
plt.bar(b[:-1],h,width=1.0,color='gray', alpha=0.5)
plt.xlabel('Number of simultaneous visual areas')
plt.ylabel('Number of experiments')
plt.ylim([0,22])
plt.yticks(ticks=np.arange(0,22,5))

M = np.mean(simultaneous_areas)
plt.plot([M, M],[0,25],'--k')

ax = plt.gca()
[ax.spines[loc].set_visible(False) for loc in ['right', 'top']]   

ds = []

for i in range(len(common_area_names)):
    N = num_units_per_area[:,i]
    N = N[N > 0]
    ds.append(N)
    
plt.subplot2grid((1,6),(0,0),colspan=4,rowspan=1)

print('Total visual area units: ' + str(np.sum(np.concatenate(ds))))

ax = plt.gca()

box_parts = plt.boxplot(ds, patch_artist=True, whis=[5,95])
plt.xticks(ticks=np.arange(len(common_area_names))+1, labels=common_area_names, fontsize=12)
plt.ylabel('Number of units per experiment')

for i, pc in enumerate(box_parts['boxes']):
    pc.set_facecolor(get_color_palette(areas[i],'seaborn'))
    pc.set_edgecolor('k')
    pc.set_alpha(0.9)
    
for i, pc in enumerate(box_parts['medians']):
    pc.set_color('k')

plt.xlim([-0.5,9.5])
plt.ylim([0,180])

ax = plt.gca()
[ax.spines[loc].set_visible(False) for loc in ['right', 'top']]   

plt.tight_layout()

# %%


