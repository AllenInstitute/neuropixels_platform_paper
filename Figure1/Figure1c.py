import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import EcephysSession 
from allensdk.brain_observatory.ecephys.stimulus_analysis import ReceptiveFieldMapping
from allensdk.brain_observatory.ecephys.stimulus_analysis import DriftingGratings


cache_directory = "/mnt/nvme0/ecephys_cache_dir_2"

manifest_path = os.path.join(cache_directory, "manifest.json")

cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

# %%

# load session data
sessions = cache.get_session_table()
session_id = 756029989
session = cache.get_session_data(session_id)

# %%

# load LFP data
dg = DriftingGratings(session)
probe_id = session.probes.index.values[3]
lfp = session.get_lfp(probe_id)

# %%

plt.figure(1711, figsize=(9,12))
plt.clf()
    
trial_number = 248

areas = ('VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam', 'LGd', 'LP')

total_units = 0

freq1 = dg.stim_table.iloc[trial_number-1].temporal_frequency
freq2 = dg.stim_table.iloc[trial_number].temporal_frequency
freq3 = dg.stim_table.iloc[trial_number+1].temporal_frequency

t1 = np.linspace(0,2*pi*freq1*2,200)
t2 = np.linspace(0,2*pi*freq2*2,200)
t3 = np.linspace(0,2*pi*freq3*2,200)
z = np.zeros((100,))

a = np.concatenate((z,t1,z,t2,z,t3,z),axis=0)

num_plots = 13

plt.subplot(num_plots,1,1)
plt.plot(1-np.sin(a))
plt.xlim([0,len(a)])
plt.ylim([-1,3])
plt.axis('off')

for area_idx, area in enumerate(areas):
    
    units = session.units[(session.units.ecephys_structure_acronym == area)].index.values
                         
    total_units += len(units)

    t = np.arange(-4,6.0,0.001) 
    
    spike_data = session.presentationwise_spike_counts(
                bin_edges = t,
                stimulus_presentation_ids = dg.stim_table.index.values[trial_number],
                unit_ids = units
            )
    
    a = spike_data.sum(dim='time_relative_to_stimulus_onset')

    x, y = np.where(np.squeeze(spike_data.data))
    
    plt.subplot(num_plots,1,area_idx+2)
    
    plt.scatter(t[x],y,s=1, c='k',marker='|',alpha=0.4)
    
    plt.xlim([np.min(t),np.max(t)])
   
    plt.axis('off')
    
# %%

plt.subplot(num_plots,1,1)
plt.title('Unit count: ' + str(total_units))

plt.subplot(num_plots,1,11)

bounds = dg.stim_table.start_time.values[trial_number] + t

selection = (session.running_speed.start_time >= bounds[0]) & \
                  (session.running_speed.start_time < bounds[-1])

running_speeds = session.running_speed[selection].velocity
starts = session.running_speed[selection].start_time
ends = session.running_speed[selection].end_time

t = np.linspace(-4,6,len(running_speeds))

plt.plot(t, running_speeds, 'gray')
plt.plot([t[0],t[-1]],[0,0],'--g')
plt.ylim([-3,75])
plt.xlim([-4,6])

plt.subplot(num_plots,1,12)

eye_info = session.get_pupil_data()

selection = (eye_info.index.values >= bounds[0]) & \
            (eye_info.index.values < bounds[-1])

pupil_x = eye_info.pupil_center_x[selection]
pupil_y = eye_info.pupil_center_y[selection]
pupil_width = eye_info.pupil_width[selection]

t = np.linspace(-4,6,len(pupil_width))


plt.subplot(num_plots,1,12)
plt.plot(t, pupil_width, 'gray')
plt.xlim([-4,6])

lfp_plot = lfp.loc[dict(time=slice(np.min(bounds),np.max(bounds)))]

D = lfp_plot.data[:,40] * 1e6
plt.subplot(num_plots,1,10)

t = np.linspace(-4,6,len(D))

plt.plot(t,D,'-g',linewidth=0.8)

plt.xlim([np.min(t), np.max(t)])

plt.ylim([-1500,1500])


