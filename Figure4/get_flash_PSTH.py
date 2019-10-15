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
from allensdk.brain_observatory.ecephys.stimulus_analysis.flashes import Flashes

# %%

cache_dir = '/mnt/hdd0/cache_dir_10_03' 

manifest_path = os.path.join(cache_dir, "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

import warnings
warnings.filterwarnings("ignore")

# %%
sessions = cache.get_sessions()

# %%

psths = []
all_unit_ids = []
all_delays = []

plt.figure(1910)
plt.clf()

for session_id, row in sessions.iterrows():
    
    print(session_id)
    
    nwb_path = cache_dir + '/session_' + str(session_id) + '/session_' + str(session_id) + '.nwb'
    
    session = EcephysSession.from_nwb_path(nwb_path, api_kwargs={
        "amplitude_cutoff_maximum": np.inf,
        "presence_ratio_minimum": -np.inf,
        "isi_violations_maximum": np.inf,
        "filter_by_validity": True
    })
    
    fl = Flashes(session)
  
    print('  loading session...')
    unit_ids = session.units.index.values
    
    print('  calculating psth...')
    dataset = session.presentationwise_spike_counts(np.arange(0,0.25,0.001),
                                                    fl.stim_table.index.values,
                                                    session.units.index.values)
    print('  calculating metrics...')
    ttf = [fl._get_time_to_first_spike(unit, fl._get_preferred_condition(unit))
                                         for unit in unit_ids]
    
    psths.append(dataset.data)
    all_unit_ids.append(unit_ids)
    all_delays.append(np.array(ttf))
    
    h,b = np.histogram(np.array(ttf), bins= np.linspace(0,0.2,30))
    plt.plot(b[:-1],h)
    plt.show()

# %%
    
total_psth = np.concatenate(psths, axis=2)

# %%
unit_ids = np.concatenate(all_unit_ids)
ttf = np.concatenate(all_delays)

# %%

    selection = (df.structure_acronym == area) #& \
    
    num_per_area[area_idx] = np.sum(selection)
    
    selection &= (df.on_screen_rf < 0.01) #& \
    
    num_with_rfs[area_idx] = np.sum(selection)
    
    selection &= (df.area_rf < 2500)
    selection &= (df.snr > 1)
    selection &= (df.firing_rate_dg > 0.1)

# %%

plt.figure(421711)
plt.clf()


areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')


from scipy.ndimage.filters import gaussian_filter1d

num_units = np.zeros((len(areas),))

for area_idx, area in enumerate(areas):
    
    selection = (df.structure_acronym == area) & \
                (df.on_screen_rf < 0.01) & \
                (df.time_to_first_spike_fl < 0.1) & \
                (df.area_rf < 2500) & \
                (df.firing_rate_dg > 0.1) & \
                (df.snr > 1)

    index_values = np.isin(unit_ids, df[selection].index.values)
                
    num_units[area_idx] = np.sum(selection)

    psth = np.mean(np.mean(total_psth[:,:,index_values],0),1) # average over trials
    
    psth = psth - np.mean(psth[:20])
    psth = gaussian_filter1d(psth,2)
    psth = psth[0:120] / 0.001
    
    plt.plot(psth, color=get_color_palette(area, 'seaborn'))
    
plt.legend(areas)
    
    

# %%