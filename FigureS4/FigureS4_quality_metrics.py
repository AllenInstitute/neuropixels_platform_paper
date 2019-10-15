# plot distributions

# note: use 'get_data_from_warehouse' to generate the 'df' object

area_list = [['VISp', 'VISl', 'VISrl', 'VISam', 'VISpm', 'VIS', 'VISal','VISmma','VISmmp','VISli'],
             ['LGd','LD', 'LP', 'VPM', 'TH', 'MGm','MGv','MGd','PO','LGv','VL',
              'VPL','POL','Eth','PoT','PP','PIL','IntG','IGL','SGN','VPL','PF','RT'],
             ['CA1', 'CA2','CA3', 'DG', 'SUB', 'POST','PRE','ProS','HPF'],
             ['MB','SCig','SCiw','SCsg','SCzo','PPT','APN','NOT','MRN','OP','LT','RPF','CP']
             ]
             
area_names = ['cortex', 'thalamus', 'hippocampus','midbrain']

colors = ['#08858C','#FC6B6F', '#7ED04B','#FC9DFE']

plt.figure(1121471)
plt.clf()

metrics = [
        'firing_rate',
        'presence_ratio',
        'max_drift', 
        'amplitude',
        'spread',
        'duration',
        'isi_violations',
        'snr',
        'isolation_distance',
        'd_prime',
        'amplitude_cutoff',
        'nn_hit_rate'
        ]


labels = [
        'Overall firing rate (Hz)',
        'Presence ratio',
        'Maximum drift ($\mu$m)',
        'Waveform amplitude ($\mu$V)',
        'Waveform spread ($\mu$m)',
        'Waveform duration (ms)',
        'ISI violations',
        'SNR',
        'Isolation distance',
        "$d'$",
        'Amplitude cutoff',
        'Nearest-neighbors hit rate',
        
          ]



bins = [
       np.linspace(-3,2,100), 
        np.linspace(0,1,50),
        np.linspace(0,120,100), 
        np.linspace(0,500,100), 
        np.linspace(0,200,50),
        np.linspace(0,1.15,80),
        np.linspace(-5,3,100),
        np.linspace(0,8,100),
        np.linspace(0,200,100), 
        np.linspace(0,12,80),
        np.linspace(0,0.5,100),
        np.linspace(0.01,1,100)
        
        ]

use_log = [
        True,
        False,
        False, 
        False, 
        False,
        False,
        True,
        False,
        False,
        False,
        False,
        False
           ]

boundaries = [None,0.95, None, None, None, None, np.log10(0.5), None, None, None, 0.1, None ]

unit_count = np.zeros((len(area_list),))
max_values = np.zeros((len(area_list),))

N_rows = 4
N_cols = np.ceil(len(metrics)/N_rows)

for metric_idx, metric in enumerate(metrics):
        
    ax = plt.subplot(N_rows,N_cols,metric_idx+1)

    for area_idx, AREAS in enumerate(area_list):
        
        D = df[(df.structure_acronym.isin(AREAS)) & (df.quality == 'good')][metric]
        
        if metric == 'duration':
            D = D + np.random.rand(len(D))*0.02
        elif metric == 'spread':
            D = D + np.random.rand(len(D)) * 10
        
        if use_log[metric_idx]:
            h,b = np.histogram(np.log10(D+1e-4), bins=bins[metric_idx], density=True)
        else:
            h,b = np.histogram(D, bins=bins[metric_idx], density=True)
            
        if metric_idx == 0:
            unit_count[area_idx] = len(D)
            
        filter_win = 1
            
        if metric == 'isi_violations': # or metric == 'l_ratio':
            plt.bar(b[12]-1.2+area_idx*0.4,h[12],width=0.3,color=colors[area_idx])
            
            x = b[14:]
            y = gaussian_filter1d(h[13:],filter_win)
            max_values[area_idx] = h[12]
            
        elif metric == 'amplitude_cutoff': # or metric == 'l_ratio':
            plt.bar(0.53 + area_idx*0.03,h[-1],width=0.025,color=colors[area_idx])
            
            x = b[:-2]
            y = gaussian_filter1d(h[:-1],0.01)
            max_values[area_idx] = h[0]
 
        else:
            
            x = b[:-1]
            y = gaussian_filter1d(h,filter_win)
            max_values[area_idx] = np.max(y)
            
        plt.plot(x,y,color=colors[area_idx], linewidth=2.0)    
        
        if metric == 'isi_violations':
            plt.xticks(ticks=[-4.632,-3,-1,1],labels=['0','0.001','0.1','10'])
        elif metric == 'firing_rate':
            plt.xticks(ticks=[-2,-1,0,1],labels=['0.01','0.1','1','10'])
    
    plt.xlabel(labels[metric_idx])
    plt.gca().get_yaxis().set_visible(False)
    plt.ylim([0, np.max(max_values)*1.1])
    [ax.spines[loc].set_visible(False) for loc in ['right', 'top', 'left']]
    
    if boundaries[metric_idx] != None:
        plt.plot([boundaries[metric_idx], boundaries[metric_idx]], [0, np.max(max_values)*1.1], ':k')
        
plt.subplot(N_rows,N_cols, N_cols)

unit_count_string = [str(int(np.floor(uc/1e3))) + ',' + str(int(uc % 1000)).zfill(3) for uc in unit_count]

annotations = [area_names[i] + ' ($N$ = ' + unit_count_string[i] + ')' for i in range(len(unit_count))]
plt.legend(annotations)
        
plt.tight_layout()
            
# %%

# plot filtering steps


area_list = [['VISp', 'VISl', 'VISrl', 'VISam', 'VISpm', 'VIS', 'VISal','VISmma','VISmmp','VISli'],
             ['LGd','LD', 'LP', 'VPM', 'TH', 'MGm','MGv','MGd','PO','LGv','VL',
              'VPL','POL','Eth','PoT','PP','PIL','IntG','IGL','SGN','VPL','PF','RT'],
             ['CA1', 'CA2','CA3', 'DG', 'SUB', 'POST','PRE','ProS','HPF'],
             ['MB','SCig','SCiw','SCsg','SCzo','SCop','PPT','APN','NOT','MRN','OP','LT','RPF'],
             ['CP','ZI','grey','COAa',nan,'BMAa']
             ]
             
area_names = ['cortex', 'thalamus', 'hippocampus','midbrain','other','unlabeled']


main_areas = ['VISp','VISl','VISrl','VISam','VISpm','VISal','LP','LGd']

colors = ['#08858C','#FC6B6F', '#7ED04B','#FC9DFE','#5D5D5D','#C5C5C5']


values = np.zeros((5,6))

for area_idx, AREAS in enumerate(area_list):
    
     filter1 = all_metrics.structure_acronym.isin(AREAS)
     
     values[area_idx, 0] = np.sum(filter1)
     
     filter2 = df.quality == 'good'
     
     values[area_idx, 1] = np.sum(filter1 & filter2)
     
     filter3 = (df.on_screen_rf < 0.01)
               
     values[area_idx, 2] = np.sum(filter1 & filter2 & filter3)
     
     filter4 = (df.structure_acronym.isin(main_areas))
     
     values[area_idx, 3] = np.sum(filter1 & filter2 & filter3 & filter4)
     
     filter5 = (df.area_rf < 2500) & \
               (df.snr > 1) & \
               (df.firing_rate_dg > 0.1)
               
     values[area_idx, 4] = np.sum(filter1 & filter2 & filter3 & filter4 & filter5)
     
     filter6 = (df.time_to_first_spike_fl < 0.1)

     values[area_idx, 5] = np.sum(filter1 & filter2 & filter3 & filter4 & filter5 & filter6)
     
    
plt.figure(4711)
plt.clf()

totals = np.sum(values, 0)

levels = np.cumsum(values, 0)

for i in range(levels.shape[0]-1,-1,-1):
    
    plt.plot(levels[i,:],color=colors[i])
    plt.plot(levels[i,:],'.',color=colors[i])
    
for i in range(levels.shape[1]):
    plt.plot([i,6],[totals[i],totals[i]],':k',alpha=0.5)
    plt.text(6,totals[i]-1000,str(int(totals[i])))

plt.xlim([0,6.5])

plt.ylim([0,140000])
        

# %%

    
    