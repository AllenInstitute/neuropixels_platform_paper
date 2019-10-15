
from scipy.ndimage.filters import gaussian_filter1d
from sklearn.utils import resample    


# %%
code_directory = '/home/joshs/GitHub/neuropixels_platform_paper'

intrinsic_timescale = np.load(os.path.join(code_directory, 'data', '/intrinsic_timescale.npy'))

# %%

# %%
plt.figure(14781)
plt.clf()

# %%
areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')

color_palette= 'seaborn'

hierarchy_score = {'LGd' : -0.515,
                   'VISp' : -0.357,
                   'VISl' : -0.093,
                   'VISrl' : -0.059,
                   'LP' : 0.105,
                   'VISal' : 0.152,
                   'VISpm' : 0.327,
                   'VISam' : 0.441}

HS = [-0.515, -0.357, -0.093, -0.059, 0.105, 0.152,0.327, 0.441]

num_units = 0

def get_bootstrap_95ci(M, measure_of_central_tendency, N=1000):
    n = int(len(M)/2)
    est = np.zeros((N,))
    for i in range(N):
        boot = M[np.random.permutation(len(M))[:n]]
        est[i] = measure_of_central_tendency(boot)
        
    return np.percentile(est,97.5) - np.nanmean(est)

from scipy import stats

def convert_to_ms(value_in_s):
    return value_in_s*1000

def take_log(original_value):
    return np.log10(original_value)

def do_not_change(original_value):
    return original_value

measure_of_central_tendency = np.nanmean

np.random.seed(10)

num_per_area = np.zeros((8,))
num_with_rfs = np.zeros((8,))
num_after_filter = np.zeros((8,))
num_after_fl = np.zeros((8,))
mice_per_area = np.zeros((8,))

if True:
    metrics = ['time_to_first_spike_fl', 'area_rf', 'mod_idx_dg', 'autocorrelation_timescale']
    labels = ['Time to first spike (ms)', 'RF area ($deg^2$)', '$log_{10}$ modulation index', 'Autocorrelation timescale (ms)']
    bins = [np.linspace(25,120,30), np.linspace(10,2000,30), np.linspace(-1.5,2,50), np.linspace(0,400,100)]
    function_to_apply = [convert_to_ms, do_not_change, take_log, do_not_change]
    y_vals = [60, 520, 0.3, 60]

else:
    metrics = [ 'firing_rate', 'lifetime_sparseness_ns', 'g_osi_dg']
    labels = ['$log_{10}$ Firing rate', 'lifetime sparseness', 'OSI']
    bins = [np.linspace(-1,2), np.linspace(0,1), np.linspace(0,1)]
    function_to_apply = [take_log, do_not_change, do_not_change]
    y_vals = [0.75, 0.2, 0.1]

centers = np.zeros((8,len(metrics)))
errorbars = np.zeros((8,len(metrics)))

max_values = np.zeros((len(metrics),))

all_values = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}}

for area_idx, area in enumerate(areas):
    
    selection = (df.structure_acronym == area) #& \
    
    num_per_area[area_idx] = np.sum(selection)
    
    selection &= (df.on_screen_rf < 0.01) #& \
    
    num_with_rfs[area_idx] = np.sum(selection)
    
    selection &= (df.area_rf < 2500)
    selection &= (df.snr > 1)
    selection &= (df.firing_rate_dg > 0.1)

    num_after_filter[area_idx] = np.sum(selection)
                   
    mice_per_area[area_idx] = len(df[selection].specimen_id.unique())
    
    for metric_idx, metric in enumerate(metrics):
        
        if metric_idx == 0:
            selection &= (df.time_to_first_spike_fl < 0.1) 
            num_after_fl[area_idx] = np.sum(selection)

        if metric_idx != 3:
            M = function_to_apply[metric_idx](df[selection][metric].values) 
        else:
            M = intrinsic_timescale[area_idx]
        
        all_values[metric_idx][area] = M
        h, b = np.histogram(M, bins=bins[metric_idx], density=True)
        
        if metric_idx == 3:
            _filter = 6
        else:
            _filter = 1.5
    
        h_filt = gaussian_filter1d(h,_filter)
        
        max_values[metric_idx] = np.max([np.max(h_filt), max_values[metric_idx]])

        plt.subplot(len(metrics), 4, metric_idx*4+1)
        plt.plot(b[:-1],h_filt,color=get_color_palette(areas[area_idx], color_palette))
        plt.xlabel(labels[metric_idx])
        
        plt.subplot(len(metrics), 4, metric_idx*4+2)
        plt.plot(b[:-1],np.cumsum(h_filt),color=get_color_palette(areas[area_idx], color_palette))
        plt.xlabel(labels[metric_idx])
        
        centers[area_idx, metric_idx] = measure_of_central_tendency(M)
        errorbars[area_idx,  metric_idx] = get_bootstrap_95ci(M, measure_of_central_tendency) #np.nanstd(M) / np.sqrt(len(M))
        
print('TOTAL: ' + str(np.sum(num_per_area)))
    
from scipy.stats import linregress, pearsonr, spearmanr


x = HS
    
for i in range(len(metrics)):
    
    plt.subplot(len(metrics),4,i*4+3)
    y = centers[:,i]
    
    if i == 0:
        stop
        
    
    for k in range(8):
        plt.plot(x[k], centers[k,i],'.',color=get_color_palette(areas[k], color_palette))
        plt.errorbar(x[k], centers[k,i], yerr = errorbars[k,i], fmt='.',color=get_color_palette(areas[k], color_palette),alpha=0.8)

    slope,intercept,r,p,std = linregress(x,y)
    x2 = np.linspace(-0.75,0.5,10)
    
    plt.plot(x2,x2*slope+intercept,'--k', alpha=0.5)
    
    r_s,p_s = spearmanr(x,y)
    r_p,p_p = pearsonr(x,y)
    
    text =  '$r_P$ = ' + str(np.around(pow(r_p,1),2)) + '; $P_P$ = ' + str(np.around(p_p,6)) + '\n' + \
            '$r_S$ = ' + str(np.around(pow(r_s,1),2)) + '; $P_S$ = ' + str(np.around(p_s,6))
    
    plt.text(-0.30,y_vals[i],text,fontsize=8)
    plt.ylabel(labels[i])
        
for i in range(len(metrics)):
    for j in range(2):
        plt.subplot(len(metrics),4,i*4+1+j)
        ax = plt.gca()
        plt.gca().get_yaxis().set_visible(False)
        [ax.spines[loc].set_visible(False) for loc in ['right', 'top', 'left']]        
    
    plt.subplot(len(metrics),4,i*4+3)
    ax = plt.gca()
    [ax.spines[loc].set_visible(False) for loc in ['right', 'top']]   
    plt.xlim([-0.85,0.65])
        

from scipy.stats import ks_2samp, ranksums
from statsmodels.stats.multitest import multipletests

alpha = 0.05

common_names = ['LGN', 'V1', 'LM', 'RL', 'LP', 'AL', 'PM', 'AM']

for metric_idx, metric in enumerate(metrics):
    
    comparison_matrix = np.zeros((len(areas),len(areas))) #+ 1000
    
    for area_idx1, area1 in enumerate(areas):
        
        for area_idx2, area2 in enumerate(areas):
            
            if area_idx2 > area_idx1:
                
                v1 = all_values[metric_idx][area1]
                v2 = all_values[metric_idx][area2]
                
                z, p = ranksums(v1[np.invert(np.isnan(v1))],
                                   v2[np.invert(np.isnan(v2))] )
                comparison_matrix[area_idx1, area_idx2] = p + 1e-5
       
    p_values = comparison_matrix.flatten()
    ok_inds = np.where(p_values > 0)[0]
    inds = np.where(comparison_matrix > 0)
    indx = inds[0]
    indy = inds[1]
    
    reject, p_values_corrected, alphaSidak, alphacBonf = multipletests(p_values[ok_inds], alpha=alpha, method='fdr_bh')
            
    p_values_corrected2 = np.zeros((len(p_values),))
    p_values_corrected2[ok_inds] = p_values_corrected
    comparison_corrected = np.reshape(p_values_corrected2, comparison_matrix.shape)
    
    sig_thresh = np.log10(alpha)
    plot_range = 8
    
    plt.subplot(len(metrics),4,metric_idx*4+4)
    plt.imshow(np.log10(comparison_corrected),cmap='bone',vmin=-5,vmax=np.log10(0.05))
    
    plt.colorbar(fraction=0.026, pad=0.04)
    
    plt.xticks(ticks=np.arange(len(common_names)), labels=common_names)
    plt.yticks(ticks=np.arange(len(common_names)), labels=common_names)
    plt.ylim([-0.5,len(common_names)-0.5])
    plt.xlim([-0.5,len(common_names)-0.5])    
    
plt.tight_layout()