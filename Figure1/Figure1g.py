from scipy.ndimage.filters import gaussian_filter1d
from scipy import stats

df = pd.read_csv(os.path.join(os.getcwd(), 'data', 'unit_table.csv'), low_memory=False)

# %%

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
            
            selection = (df.ecephys_structure_acronym.isin(['CA1','CA3','DG'])) & \
                        (df[key] == mouse)
        else:
            selection = (df.ecephys_structure_acronym == area) & \
                        (df[key] == mouse)
                        
        values[mouse_idx, area_idx, 0] = np.median(df[selection].firing_rate)
        
        total_units += np.sum(selection)
    
        selection2 = selection * (df.p_value_rf < 0.01)
        
        values[mouse_idx, area_idx, 1] = np.sum(selection2) / np.sum(selection)
        
        values[mouse_idx, area_idx, 2] = np.sum(selection2)
        values[mouse_idx, area_idx, 3] = np.sum(selection)
        
# %%

M = np.sum(values[:,:8,2]/np.sum(values[:,:8,3]))
print('Percentage with significant RFs: ' + str(M*100))

plt.figure(15000)
plt.clf()

common_area_names = ['V1','LM','RL','AL','PM','AM', 'LGN', 'LP','HPC']

i = 1

B = np.nanmean(values[:,:,i],0)
C = np.nanstd(values[:,:,i],0)

for j in range(len(common_area_names)):
    
    V = values[:,j,i]
    V = V[V>0]
    plt.scatter(np.random.rand(len(V))*0.25-0.125+j,
                V,
                edgecolors=get_color_palette(areas[j], 'seaborn'),
                facecolors='none',alpha=0.5)
    plt.bar(j, B[j], width=0.6, edgecolor=get_color_palette(areas[j], 'seaborn'), 
            facecolor='none')
    
plt.errorbar(np.arange(len(common_area_names)), B, C, fmt='.',color='k')
    
plt.xticks(ticks=np.arange(len(areas)), labels=common_area_names)

plt.ylabel('Fraction of units with RFs')
    
plt.tight_layout()

