import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy import signal


def autocorr(x):
    """
    x is 1 d
    """
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]

def autocorr2D(x):
    
    corr = signal.correlate(x, x, mode='same')
    return corr

def load_npz(filename):
    """
    load npz files with sparse matrix and dimension
    output dense matrix with the correct dim
    """
    npzfile = np.load(filename, allow_pickle=True) 
    sparse_matrix = npzfile['arr_0'][0]
    ndim=npzfile['arr_0'][1]

    new_matrix_2d = np.array(sparse_matrix.todense())
    new_matrix = new_matrix_2d.reshape(ndim)
    return new_matrix


def get_bootstrap_95ci(M, N=5000):
    n = int(len(M)/2)
    est = np.zeros((N,))
    for i in range(N):
        boot = M[np.random.permutation(len(M))[:n]]
        est[i] = np.nanmedian(boot)
        
    return np.percentile(est,97.5) - np.nanmedian(est)

def get_color(name='seaborn'):
	import cmocean  
	name='seaborn'
	if name=='cmocean':
	    hierarchy_colors = cmocean.cm.phase(np.arange(1.0,0.1,-0.124))           
	    color_bank = {'LGd' : hierarchy_colors[0],
	                  'VISp' : hierarchy_colors[1],
	                  'VISl' : hierarchy_colors[2],
	                  'VISrl' : hierarchy_colors[3],
	                  'LP' : hierarchy_colors[4],
	                  'VISal' : hierarchy_colors[5],
	                  'VISpm' : hierarchy_colors[6],
	                  'VISam' : hierarchy_colors[7],
	                 }

	if name=='seaborn':
	    colors = [[217,141,194],
	              [129,116,177],
	              [78,115,174],
	              [101,178,201],
	              [88,167,106],
	              [202,183,120],
	              [219,132,87],
	              [194,79,84]]

	    def scale_colors(color):
	        return [col/255. for col in color]

	    hierarchy_colors = [scale_colors(col) for col in colors]

	    color_bank = {
	            'VISp' : hierarchy_colors[1],
	            'VISl' : hierarchy_colors[2],
	            'VISal' : hierarchy_colors[5],
	            'VISrl' : hierarchy_colors[3],
	            'VISpm' : hierarchy_colors[6],
	            'VISam' : hierarchy_colors[7],


	            'DG' : '#A4A4A4',
	            'CA3' : '#6D6D6D',
	            'CA1' : '#5B5B5B',
	            'CA2' : '#5B5B5B',
	            'CA' : '#7ED04B',
	            'POST' : '#A4A4A4',
	            'SUB' : '#A4A4A4',
	            'HPC' : '#A4A4A4',

	            'LGd' : hierarchy_colors[0],
	            'LP' : hierarchy_colors[4]
	            }
	return color_bank

def define_areas_hierarchy_score():
	areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')
	#areas = ('VISp','VISl','VISrl','VISal','VISpm','VISam')
	HS = [-0.515, -0.357, -0.093, -0.059, 0.105, 0.152,0.327, 0.441]
	HSA = [-0.6329, -0.4209, -0.08555, -0.054969, 0.17871887, 0.0226128859, 0.112409304, 0.274742343]
	HS_C = [-0.357, -0.093, -0.059, 0.152,0.327, 0.441]
	return areas, HS, HSA, HS_C


areas_all = ['LGd', #LGN
         'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'VISmma', #VIS
         'LP','TH', # thalamus
         'MB','MGm', 'MGv', 'MGd', 'APN','Eth','POL','ProS' ,'NOT', 'VPM','RPF','SGN', 'PRE', 'POST', 'SUB', 'HPF', 'ZI', 'IntG', # others
         'DG','CA1', 'CA3', # hippo
         'none'   # not labeled
        ]

# select mouse IDs
df_exp=pd.read_csv('/Volumes/local1/work_allen/Ephys/resorted/experiment_table_2019-08-27.csv', index_col=0)
df_exp['stimulus_set'].unique()

mouseIDs = ['mouse'+str(i) for i in df_exp.index.values]
print(len(mouseIDs))

# process autocorrelation in each mouse each area
AMO=[]
SEP=[]
m_id=[]
sep_start=0
probes = []
for ii, mouseID in enumerate(mouseIDs):
    print(mouseID)
    # 1. load spikes
    basepath = '/Volumes/local1/work_allen/Ephys/resorted/'+mouseID
    if os.path.isfile('/Volumes/local1/work_allen/Ephys/resorted/'+mouseID+'/matrix/flash_all_units.npz'):
        spikes = load_npz(basepath+'/matrix/flash_all_units.npz')
        # add area without layer
        # with the new file, can't load with index_col
        df_old = pd.read_csv(basepath+'/matrix/'+mouseID+'_all_units_meta.csv')
        if 'b' in df_old.ccf.unique()[0]:
            df_old['ccf'] = df_old['ccf'].apply(lambda x: x[2:-1])
        #print(df_old.ccf.unique())
        
        df_old['areas_group']=np.zeros(len(df_old))
        for a in areas_all:
            df_old['areas_group'][df_old['ccf'].str.contains(a, case=False, na=False, regex=False)]=a
        # select units in cortex and LP and LGD
        idx=[]
        for probe in df_old.probe_id.unique():
            df_tmp = df_old[(df_old.probe_id==probe) & (df_old.snr>1) & (df_old.qc_amp<0.1)]
            areas = df_tmp.areas_group.unique()
            for a in areas:
                if 'VIS' in str(a) or 'LP' in str(a) or 'LGd' in str(a):
                    idx.append(df_tmp[df_tmp.areas_group==a].index.values)
        idx = [item for sublist in idx for item in sublist]         

        df=df_old.iloc[idx]
        df = df.reset_index().drop(['index'], axis=1)
        spikes = spikes[idx, :, :, :]
        
        FR = df.FR.values
        assert len(df)==spikes.shape[0]

        #constrain by RF on screen
        # load RF fit
        df_rf = pd.read_csv('/Volumes/local1/work_allen/Ephys/processed_data/RF_features/resorted/'+mouseID+'_rf_features.csv', index_col=0)    
        rf_index=[]
        for index, row in df.iterrows():
            probe=row['probe_id']
            unit=row['unit_id']
            df_rf_tmp=df_rf[df_rf.probe_id==probe]
            # RF fit exist
            if unit in df_rf_tmp.unit_id.values:
                pos = df_rf[(df_rf.probe_id==probe) & (df_rf.unit_id==unit)].as_matrix(columns = ('rf_center_x1','rf_center_y1'))
                pos = pos.astype(float)[0]
                # RF center on screen
                if pos[0]>0 and pos[0]<8 and pos[1]>0 and pos[1]<8:
                    rf_index.append(index)
        spikes_new=spikes[rf_index,:,:,:]
        df_new=df.iloc[rf_index]
        df_new = df_new.reset_index().drop(['index'], axis=1)

        probenames = df_new.probe_id.unique().astype(str)
        separations = [0]
        #separations[-1]=separations[-1]-1
        for probe in probenames:
            index = np.where(df_new.probe_id==probe)[0]
            separations = np.concatenate([separations, [index[-1]+1]],axis=0)

        ACCG=np.zeros((spikes_new.shape[0], spikes_new.shape[-1]))
        for i in range(spikes_new.shape[0]):
            x = spikes_new[i,1,:,:]
            tmp = autocorr2D(x)
            ACCG[i,:]=tmp.mean(0)
        
        print(ACCG.shape)
        
        if ii==0:
            AMO=ACCG
            probes = df_new.areas_group
            m_id.append([mouseID]*len(df_new))
        else:
            AMO=np.concatenate([AMO, ACCG], axis=0)
            probes = np.concatenate([probes, df_new.areas_group], axis=0)
            m_id.append([mouseID]*len(df_new))

        separations=separations+sep_start
        SEP.append(separations)
        sep_start=separations[-1]
m_id = [item for sublist in m_id for item in sublist] 

# plot population ACCG across mouse

areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')

m = []
for i, probe in enumerate(areas):
    #plt.plot(np.nanmean(AMO[np.where(probes==probe)[0],:], axis=0), c=color_bank[probe])
    m.append(max(np.nanmean(AMO[np.where(probes==probe)[0],:], axis=0)))   

for i, probe in enumerate(areas):
    plt.plot(np.nanmean(AMO[np.where(probes==probe)[0],:], axis=0)/m[i], c=color_bank[probe], label=probe)
plt.legend()
#plt.xlim([500, 1000])


