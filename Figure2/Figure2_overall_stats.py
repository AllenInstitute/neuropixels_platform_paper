from scipy.ndimage.filters import gaussian_filter1d
from scipy import stats

# note: use 'get_data_from_warehouse' to generate the 'df' object

areas = ('VISp','VISl','VISrl','VISal','VISpm','VISam','LGd', 'LP', 'HPC')

color_palette= 'seaborn'

key = 'specimen_id'

mice = df[key].unique()

values = np.zeros((len(mice), len(areas), 4))
values[:] = np.nan

np.random.seed(10)

total_units = 0

for mouse_idx, mouse in enumerate(mice):

    for area_idx, area in enumerate(areas):
        
        if area == 'HPC':
            
            selection = (df.structure_acronym.isin(['CA1','CA3','DG'])) & \
                        (df[key] == mouse)
        else:
            selection = (df.structure_acronym == area) & \
                        (df[key] == mouse)
                        
        values[mouse_idx, area_idx, 0] = np.median(df[selection].firing_rate)
        
        total_units += np.sum(selection)
    
        selection2 = selection * (df.on_screen_rf < 0.01)
        
        values[mouse_idx, area_idx, 1] = np.sum(selection2) / np.sum(selection)
        
        values[mouse_idx, area_idx, 2] = np.sum(selection2)
        values[mouse_idx, area_idx, 3] = np.sum(selection)
        
# %%
print('Percent significant RFs:')

M = np.sum(values[:,:8,2]/np.sum(values[:,:8,3]))
print(M*100)

plt.figure(15000)
plt.clf()

common_area_names = ['V1','LM','RL','AL','PM','AM', 'LGN', 'LP','HPC']

labels = ('Overall firing rate', 'Fraction of units with RFs') 

for i in range(2):
    
    plt.subplot(2,1,i+1)
    
    B = np.nanmean(values[:,:,i],0)
    C = np.nanstd(values[:,:,i],0)
    
    for j in range(len(common_area_names)):
            plt.bar(j, B[j], width=0.6, color=get_color_palette(areas[j], 'seaborn'))
        
    plt.errorbar(np.arange(len(common_area_names)), B, C, fmt='.',color='k')
        
    plt.xticks(ticks=np.arange(len(areas)), labels=common_area_names)
    
    plt.ylabel(labels[i])
    
plt.tight_layout()