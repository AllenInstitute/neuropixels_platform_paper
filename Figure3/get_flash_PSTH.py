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
from allensdk.brain_observatory.ecephys.stimulus_analysis import Flashes
from allensdk.brain_observatory.ecephys.stimulus_analysis.drifting_gratings import DriftingGratings

# %%

cache_dir = '/mnt/nvme0/ecephys_cache_dir' 

manifest_path = os.path.join(cache_dir, "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

import warnings
warnings.filterwarnings("ignore")

# %%
sessions = cache.get_session_table()

# %%

all_unit_ids = []
all_delays = []
all_delays_baseline = []
all_delays_dg = []

for session_id, row in sessions.iloc[3:].iterrows():
    
    print(session_id)
    
    nwb_path = cache_dir + '/session_' + str(session_id) + '/session_' + str(session_id) + '.nwb'
    
    session = EcephysSession.from_nwb_path(nwb_path, api_kwargs={
        "amplitude_cutoff_maximum": np.inf,
        "presence_ratio_minimum": -np.inf,
        "isi_violations_maximum": np.inf,
        "filter_by_validity": True
    })
    
    stop
    
    fl = Flashes(session)
    
    stop
    
    # %%
    
    v1_units = session.units[session.units.ecephys_structure_acronym == 'VISp'].index.values
    # %%
    ratio = fl._get_on_off_ratio(v1_units[1])
    print(ratio)
    
    # %%
    dg = DriftingGratings(session)
    
    print('  calculating time to first spike...')
    unit_ids = session.units.index.values
    
    ttf_baseline = calc_time_to_first_spike(
                        session.presentationwise_spike_counts(np.arange(1.53,1.7,0.001),
                                                    fl.stim_table.index.values,
                                                    session.units.index.values,
                                                    binarize=True
                                                    ))
    
    ttf_actual_fl = calc_time_to_first_spike(
                        session.presentationwise_spike_counts(np.arange(0.03,0.2,0.001),
                                                    fl.stim_table.index.values,
                                                    session.units.index.values,
                                                    binarize=True
                                                    ))
    
    ttf_actual_dg = calc_time_to_first_spike(
                        session.presentationwise_spike_counts(np.arange(0.03,0.2,0.001),
                                                    dg.stim_table.index.values,
                                                    session.units.index.values,
                                                    binarize=True
                                                    ))
    
    all_unit_ids.append(unit_ids)
    all_delays.append(ttf_actual_fl)
    all_delays_baseline.append(ttf_baseline)
    all_delays_dg.append(ttf_actual_dg)
    
    #stop
    
# %%
    
def calc_time_to_first_spike(data):
    
 b = data.where(data > 0)
 b[:,-1,:] = 1
 
 a = np.nanargmin(b, axis=1)
 a = a.astype('float')
 
 a[a == np.max(a)] = np.nan

 c = np.nanmedian(a, 0)

 return c
    

# %%
unit_ids = np.concatenate(all_unit_ids)
tfs_fl = np.concatenate(all_delays)
tfs_baseline = np.concatenate(all_delays_baseline)
tfs_dg = np.concatenate(all_delays_dg)

# %%

ttf0 = pd.read_csv('/home/joshs/GitHub/neuropixels_platform_paper/data/time_to_first_spike.csv')

# %%

first_spike = pd.DataFrame(data={'ecephys_unit_id' : unit_ids,
                                'time_to_first_spike_control' : tfs_baseline,
                                'time_to_first_spike_new' : tfs_fl,
                                'time_to_first_spike_dg' : tfs_dg})

first_spike = first_spike.merge(ttf0, on='ecephys_unit_id')

# %%

plt.plot(first_spike.time_to_first_spike_new, first_spike.time_to_first_spike_dg, '.k', alpha=0.1)

#%%
first_spike.to_csv('/home/joshs/GitHub/neuropixels_platform_paper/data/time_to_first_spike_20200312.csv')

# %%

plt.figure(191)
plt.clf()

plt.rcParams.update({'font.size': 13})

areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')

HS = [-0.515, -0.357, -0.093, -0.059, 0.105, 0.152,0.327, 0.441]

Y1 = []
Y2 = []
Y3 = []
S1 = []
S2 = []
S3 = []
R1 = []

def get_bootstrap_95ci(M, measure_of_central_tendency, N=1000):
    n = int(len(M)/2)
    est = np.zeros((N,))
    for i in range(N):
        boot = M[np.random.permutation(len(M))[:n]]
        est[i] = measure_of_central_tendency(boot)
        
    return np.percentile(est,97.5) - np.nanmean(est)

measure_of_central_tendency = np.nanmean

from scipy.stats import linregress, pearsonr, spearmanr

for area_idx, area in enumerate(areas):
    
    selection = (df.ecephys_structure_acronym == area) #& \
    
    selection &= (df.on_screen_rf < 0.01) #& \
        
    selection &= (df.area_rf < 2500)
    selection &= (df.snr > 1)
    selection &= (df.firing_rate_dg > 0.1)
    
    selection &= (df.time_to_first_spike_fl < 0.1) 

    M1 = df[selection].time_to_first_spike_fl * 1000# * 1000 # + 30# * 1000 #+ 30
    M2 = df[selection].time_to_first_spike_control + 30
    M3 = df[selection].time_to_first_spike_dg + 30
    R = df[selection].firing_rate
    
    Y1.append(np.mean(M1))
    Y2.append(np.mean(M2))
    Y3.append(np.mean(M3))
    R1.append(np.mean(R.values))
    
    S1.append(get_bootstrap_95ci(M1.values, measure_of_central_tendency))
    S2.append(get_bootstrap_95ci(M2.values, measure_of_central_tendency))
    S3.append(get_bootstrap_95ci(M3.values, measure_of_central_tendency))
    
    
plt.subplot(1,3,1)
    
for k, area in enumerate(areas):
        
    plt.plot(HS[k], Y1[k], '.', color=get_color_palette(area, 'seaborn'))
    plt.errorbar(HS[k], Y1[k], yerr = S1[k], fmt='.',color=get_color_palette(area, 'seaborn'),alpha=0.4)
    
    plt.plot(HS[k], Y2[k], '.', color=get_color_palette(area, 'seaborn'))
    plt.errorbar(HS[k], Y2[k], yerr = S2[k], fmt='.',color=get_color_palette(area, 'seaborn'),alpha=0.4)
    
plt.xlabel('Anatomical hierarchy score')
plt.ylabel('Time to first spike (ms)')
    
def plot_best_fit(x, y, color='black', xlims=[4,12], text_loc = [10, 85]):

    slope,intercept,r,p,std = linregress(x,y)
    x2 = np.linspace(xlims[0],xlims[1],10)
    
    plt.plot(x2,x2*slope+intercept,'--', color=color, alpha=0.5)
    
    r_s,p_s = spearmanr(x,y)
    r_p,p_p = pearsonr(x,y)
    
    text =  '$r_P$ = ' + str(np.around(pow(r_p,1),2)) + '; $P_P$ = ' + str(np.around(p_p,6)) + '\n' + \
            '$r_S$ = ' + str(np.around(pow(r_s,1),2)) + '; $P_S$ = ' + str(np.around(p_s,6))

    plt.text(text_loc[0], text_loc[1],text,fontsize=11, color='black')

plot_best_fit(HS, Y1, xlims=[-0.75, 0.5], text_loc=[0.0, 65])
plot_best_fit(HS, Y2, xlims=[-0.75, 0.5], text_loc=[0.0, 85], color='gray')


plt.subplot(1,3,2)
for k, area in enumerate(areas):
        
    plt.plot(R1[k], Y1[k], '.', color=get_color_palette(area, 'seaborn'))
    plt.errorbar(R1[k], Y1[k], yerr = S1[k], fmt='.',color=get_color_palette(area, 'seaborn'),alpha=0.4)
    
    plt.plot(R1[k], Y2[k], '.', color=get_color_palette(area, 'seaborn'))
    plt.errorbar(R1[k], Y2[k], yerr = S2[k], fmt='.',color=get_color_palette(area, 'seaborn'),alpha=0.4)

plot_best_fit(R1, Y1, text_loc=[9,75])
plot_best_fit(R1, Y2, color='gray', text_loc=[9,95])

plt.xlabel('Mean firing rate')
plt.ylabel('Time to first spike (ms)')

plt.subplot(1,3,3)
for k, area in enumerate(areas):
        
    plt.plot(HS[k], Y3[k], '.', color=get_color_palette(area, 'seaborn'))
    plt.errorbar(HS[k], Y3[k], yerr = S3[k], fmt='.',color=get_color_palette(area, 'seaborn'),alpha=0.4)
    
plot_best_fit(HS, Y3, xlims=[-0.75, 0.5], text_loc=[0, 70])

plt.ylim([65, 92])

plt.xlabel('Anatomical hierarchy score')
plt.ylabel('Time to first spike (ms)')

# %%
    
    #plt.subplot(2,4,area_idx+1)
    plt.subplot(2,1,1)
    plt.plot(HS[area_idx], np.mean(M1), '.k')
    
    plt.subplot(2,1,2)
    plt.plot(HS[area_idx], np.mean(M2), '.r')
    #plt.plot(M1, M2, '.k')
    
   
   

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