#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 17:04:32 2020

@author: joshs
"""

# %%

df = pd.read_csv(os.path.join(os.getcwd(), 'data', 'unit_table.csv'), low_memory=False)

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# %%

from scipy.ndimage.filters import gaussian_filter1d

areas = ['VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam']
area_names = ['V1', 'LM', 'RL', 'AL', 'PM', 'AM']

plt.figure(14914)
plt.clf()

for area_idx, area in enumerate(areas):
    
    plt.subplot(1,6,area_idx+1)
    
    sub_df = df[df.ecephys_structure_acronym == area]
    
    h,b = np.histogram(sub_df.cortical_depth, np.linspace(0.001,1.1,10), density=True)
    bin_width = np.mean(np.diff(b))
    plt.bar(b[:-1] + bin_width/2,h, width=bin_width, alpha=0.3)
    
    h,b = np.histogram(sub_df.cortical_depth, np.linspace(0.001,1.1,50), density=True)
    bin_width = np.mean(np.diff(b))
    plt.plot(b[:-1],gaussian_filter1d(h,3))
    plt.ylim([0,1.8])
    
    plt.xlabel('Cortical depth')
    
    plt.title(area_names[area_idx])
    
    ax = plt.gca()
    [ax.spines[loc].set_visible(False) for loc in ['right', 'top']]   

    
plt.tight_layout()

# %%

plt.rcParams.update({'font.size': 11})

plt.figure(14915)
plt.clf()

def convert_to_ms(value_in_s):
    return value_in_s*1000

def take_log(original_value):
    return np.log10(original_value)

def do_not_change(original_value):
    return original_value

def get_bootstrap_95ci(M, measure_of_central_tendency, N=1000):
    n = int(len(M)/2)
    est = np.zeros((N,))
    for i in range(N):
        boot = M[np.random.permutation(len(M))[:n]]
        est[i] = measure_of_central_tendency(boot)
        
    return np.percentile(est,97.5) - np.nanmean(est)

measure_of_central_tendency = np.nanmean

metrics = ['time_to_first_spike_fl', 'area_rf', 'mod_idx_dg', 'timescale_ac']
labels = ['Time to first spike (ms)', 'RF area ($deg^2$)', '$log_{10}$ modulation index', 'Autocorrelation timescale (ms)']
bins = [np.linspace(0,250,30), np.linspace(10,2000,30), np.linspace(-1.5,2,50), np.linspace(0,400,100)]
function_to_apply = [convert_to_ms, do_not_change, take_log, do_not_change]
y_lims = [[25,105], [0,1750], [-1.5,1.25], [0,100]]
y_line = [65, 600, 0.0, 50]

layers = [2,4,5,6]
layer_pts = [[[0.1,0.35],[0.35,0.55],[0.55,0.75],[0.75,1.0]],
             [[0.1,0.3],[0.3,0.45],[0.45,0.75],[0.75,1.0]],
             [[0.1,0.35],[0.35,0.5],[0.5,0.75],[0.75,1.0]],
             [[0.1,0.3],[0.3,0.4],[0.4,0.75],[0.75,1.0]],
             [[0.1,0.35],[0.35,0.45],[0.45,0.75],[0.75,1.0]],
             [[0.1,0.35],[0.35,0.47],[0.47,0.75],[0.75,1.0]]
             ]

HS = [-0.357, -0.093, -0.059, 0.152,0.327, 0.441]

mean_values = np.zeros((6,4,4))
ci_values = np.zeros((6,4,4))
total_n = np.zeros((6,4,4))

for area_idx, area in enumerate(areas):
    
    selection = (df.ecephys_structure_acronym == area) #& \
    
    selection &= (df.p_value_rf < 0.01) #& \
    selection &= (df.cortical_depth > 0.0)
    
    selection &= (df.area_rf < 2500)
    selection &= (df.snr > 1)
    selection &= (df.firing_rate_dg > 0.1)
    
    for metric_idx, metric in enumerate(metrics):
                
        plt.subplot(4, 6, area_idx+metric_idx * 6 + 1)
        
        if metric_idx == 0:
            selection &= (df.time_to_first_spike_fl < 0.1) 
            plt.title(area)
        elif metric_idx == 3:
            selection &= (df[metric] < 300)
            selection &= (df[metric] > 1)
            selection &= (df.spike_count_ac > 50)
            selection &= (df.err_ac < 20)
            #num_after_ac[area_idx] = np.sum(selection)

        M = function_to_apply[metric_idx](df[selection][metric].values) 
            
        if area_idx == 0:
            plt.ylabel(metric)
            
        if metric_idx == 2:
            plt.xlabel('Cortical depth')

        if metric_idx == 1:
            M = M + np.random.rand(len(M)) * 100
        
        #plt.scatter(df[selection].cortical_depth, M, c=df[selection].cortical_layer, s=4.0, alpha=0.6, cmap='Spectral')
        
        y = y_line[metric_idx]
        plt.plot([0.0,1.1],[y,y],':k',alpha=0.4)
        
        for layer_idx, layer in enumerate(layers):
            
            selection2 = (df.cortical_layer == layer)
            M = function_to_apply[metric_idx](df[selection & selection2][metric].values) 
            if metric_idx == 1:
                M = M + np.random.rand(len(M)) * 100
            
            m = np.mean(M)
            total_n[area_idx, metric_idx, layer_idx] = len(M)
            
            mean_values[area_idx, layer_idx, metric_idx] = m
            ci_values[area_idx, layer_idx, metric_idx] = get_bootstrap_95ci(M, measure_of_central_tendency)
            #plt.plot(layer_pts[area_idx][layer_idx],[m, m], '-k')
            
            s = np.around(np.mean(layer_pts[area_idx][layer_idx]),2)
            w = layer_pts[area_idx][layer_idx][1] - layer_pts[area_idx][layer_idx][0]
        
            plt.boxplot(M, positions = [s], widths=[w*0.75], sym='') #, notch=True)
            
        ax = plt.gca()
        [ax.spines[loc].set_visible(False) for loc in ['right', 'top']]  
        
        plt.ylim(y_lims[metric_idx])
        plt.xlim([0.0,1.0])
        #ax.set_xticks(ticks=np.arange(0.2,0.8,0.2))
        
plt.tight_layout()


print(np.sum(total_n,0))


plt.figure(1911)
plt.clf()

from scipy.stats import linregress, pearsonr, spearmanr

YYY = [55,500, 0.1, 27]
YLIMS = [[50,80],[450,1000],[-.4,0.4],[25,62]]

layer_names = ['L2/3', 'L4', 'L5', 'L6']

for metric_idx, metric in enumerate(metrics):
    
    for layer_idx, layer in enumerate(layers):
        
        plt.subplot(4,4,metric_idx * 4+layer_idx + 1)
        
        x = HS
        y = mean_values[:,layer_idx, metric_idx]
        s = ci_values[:,layer_idx, metric_idx]
        
        slope,intercept,r,p,std = linregress(x,y)
        x2 = np.linspace(-0.5,0.6,10)
        
        plt.plot(x2,x2*slope+intercept,'--k', alpha=0.5)
        
        r_s,p_s = spearmanr(x,y)
        r_p,p_p = pearsonr(x,y)
        
        text =  '$r_P$ = ' + str(np.around(pow(r_p,1),2)) + \
                 '; $P_P$ = ' + str(np.around(p_p,6)) + '\n' + \
                '$r_S$ = ' + str(np.around(pow(r_s,1),2)) + '; $P_S$ = ' + str(np.around(p_s,6))             #';
        plt.text(0.1,YYY[metric_idx],text,fontsize=11)
        
        ax = plt.gca()
        [ax.spines[loc].set_visible(False) for loc in ['right', 'top']]  
                
        
        for k, area in enumerate(areas):
        
            plt.plot(x[k], y[k], '.', color=get_color_palette(area, 'seaborn'))
            plt.errorbar(x[k], y[k], yerr = s[k], fmt='.',color=get_color_palette(area, 'seaborn'),alpha=0.4)
    
            
        plt.ylim(YLIMS[metric_idx])
        
        if metric_idx == 0:
            plt.title(layer_names[layer_idx])
        
        if layer_idx == 0:
            plt.ylabel(metric)
            
        if metric_idx == 3:
            plt.xlabel('Hierarchy score')
            
plt.tight_layout()
        
   
from scipy.stats import ranksums

significance_values = np.zeros((6,4))
difference_values = np.zeros((6,4))

for area_idx, area in enumerate(areas):
    
    #selection = (df.ecephys_structure_acronym == area) #& \
    
    selection = (df.p_value_rf < 0.01) #& \
    selection &= (df.cortical_depth > 0.0)
    
    selection &= (df.area_rf < 2500)
    selection &= (df.snr > 1)
    selection &= (df.firing_rate_dg > 0.1)
    
    for metric_idx, metric in enumerate(metrics):
                
        if metric_idx == 0:
            selection &= (df.time_to_first_spike_fl < 0.1) 
            plt.title(area)
        elif metric_idx == 3:
            selection &= (df[metric] < 300)
            selection &= (df[metric] > 1)
            selection &= (df.spike_count_ac > 50)
            selection &= (df.err_ac < 20)

  
        selection2 = (df.cortical_layer < 5)
        
        M_upper = function_to_apply[metric_idx](df[selection & selection2][metric].values) 
        
        selection2 = (df.cortical_layer >= 5)
        
        M_lower = function_to_apply[metric_idx](df[selection & selection2][metric].values) 
        
        difference_values[area_idx, metric_idx] = np.mean(M_upper) - np.mean(M_lower)
        
        stat, p = ranksums(M_upper, M_lower)
        
        significance_values[area_idx, metric_idx] = p
        