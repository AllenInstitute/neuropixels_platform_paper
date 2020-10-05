import pandas as pd

import os

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import (
    EcephysSession, 
    removed_unused_stimulus_presentation_columns
)
from allensdk.brain_observatory.ecephys.stimulus_analysis.receptive_field_mapping import ReceptiveFieldMapping


# %%

cache_dir = '/mnt/hdd0/cache_dir_10_03' 

manifest_path = os.path.join(cache_dir, "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

import warnings
warnings.filterwarnings("ignore")

# %%
sessions = cache.get_sessions()

# %%
rfs = []
all_rf_unit_ids = []

for session_id, row in sessions.iterrows():
    
    print(session_id)
    
    nwb_path = cache_dir + '/session_' + str(session_id) + '/session_' + str(session_id) + '.nwb'
    
    session = EcephysSession.from_nwb_path(nwb_path, api_kwargs={
        "amplitude_cutoff_maximum": np.inf,
        "presence_ratio_minimum": -np.inf,
        "isi_violations_maximum": np.inf,
        "filter_by_validity": True
    })
    
    print('  Calculating RFs...')
    RFs = ReceptiveFieldMapping(session).receptive_fields
    
    rfs.append(RFs['spike_counts'].data)
    
    all_rf_unit_ids.append(RFs.unit_id.values)
    
# %%
    
    
all_RFs = np.concatenate(rfs, axis=2)

unit_ids = np.concatenate(all_rf_unit_ids)

# %%

df = df.reset_index('unit_id')

# %%

areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')

import scipy.interpolate as interpolate
from scipy.ndimage.filters import gaussian_filter

plt.figure(4171471)
plt.clf()

plt.figure(4171479)
plt.clf()

x = np.arange(9)
y = np.arange(9)
xx, yy = np.meshgrid(x, y)
xnew = np.arange(0,9,0.5)
ynew = np.arange(0,9,0.5)

for area_idx, area in enumerate(areas):
    
    selection = (df.structure_acronym == area) & \
                (df.on_screen_rf < 0.01) & \
                (df.area_rf < 2500) & \
                (df.firing_rate_dg > 0.1) & \
                (df.snr > 1) 
                
    ok = df[selection].index.values
    
    sub_rf = all_RFs[:,:,np.isin(unit_ids, ok)]
    
    new_unit_ids = unit_ids[np.isin(unit_ids, ok)]
    
    aligned_rf = np.empty((42,42,len(ok)))
    aligned_rf[:] = np.nan
    
    for i in range(sub_rf.shape[2]):
        
        rf = sub_rf[:,:,i]
        
        if not np.isnan(df.loc[new_unit_ids[i]].azimuth_rf):
            coords_y = (np.arange(19) + 12 - (df.loc[new_unit_ids[i]].azimuth_rf/10-1)*2+8).astype('int') # + 12 ).astype('int') # - 4).astype('int')
            coords_x = (np.arange(19) + 4 + (df.loc[new_unit_ids[i]].elevation_rf/10+3)*2).astype('int') # + 12 - .astype('int')(df.loc[new_unit_ids[i]].elevation_rf/10+3)*2+8
     
            f = interpolate.interp2d(x, y, rf, kind='linear')
            rf_interp = f(xnew,ynew)
            
            aligned_rf[coords_x[0]:coords_x[-1],coords_y[0]:coords_y[-1],i] = rf_interp / np.nanmax(rf_interp)
        
        
    plt.figure(4171471)
    plt.subplot(2,4,area_idx+1)
    
    M = gaussian_filter(np.nanmean(aligned_rf,2),1)
    
    plt.imshow(1-M[10:30,10:30], vmin=0.1, vmax=0.75, cmap='gray')
    
    m = np.nanmin(M)
    mm = np.nanmax(M)
    midpoint = (m + mm)/2
    
    plt.contour(M[10:30,10:30], levels=(m, midpoint,mm), colors=['r',get_color_palette(area, 'seaborn'),'b']) #, vmin=0.1, vmax=0.75, cmap='gray')
    plt.axis('off')
    

    plt.title(area)
    
    plt.figure(4171479)
    plt.contour(M[10:30,10:30], levels=(m, midpoint,mm), colors=['r',get_color_palette(area, 'seaborn'),'b']) #, vmin=0.1, vmax=0.75, cmap='gray')
    plt.xlim([0,18])
    plt.ylim([0,18])
    plt.axis('square')
    plt.grid('on')
    plt.xlim([3,16])
    plt.ylim([3,16])
    plt.xticks(ticks=np.arange(4,18,2),labels=np.arange(-30,40,10))
    plt.yticks(ticks=np.arange(4,18,2),labels=np.arange(-30,40,10))
    
    # %%
    