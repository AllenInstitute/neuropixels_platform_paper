import matplotlib.pyplot as plt
import numpy as np

def test_color_palette(name):
    
    areas = ['VISp','VISl','VISal','VISrl','VISpm','VISam','DG','CA3','CA1','POST','SUB','LGd','LP','LD','APN','MRN']
    
    plt.figure(478)
    plt.clf()
    
    for idx, area in enumerate(areas):
        val = np.random.rand(1) + 1
        plt.bar(idx, val, width=0.5, color=get_color_palette(area, name=name))

# %%

def get_color_palette(area, name='Allen CCF'):
    
    default_color = '#C9C9C9'
    
    if name == 'Steinmetz':
        
        palette = {'VISp' : '#1d2836',
                'VISl' : '#2d405a',
                'VISal' : '#4c6d9e',
                'VISrl' : '#3c587f',
                'VISpm' : '#56769d',
                'VISam' : '#7395bc',

                'DG' : '#432135',
                'CA3' : '#703457',
                'CA1' : '#9a4376',
                'CA' : '#9a4376',
                'POST' : '#bc568e',
                'SUB' : '#d26d9c',
                
                'LGd' : '#3c6636',
                'LP' : '#73ad6c',
                'LD' : '#31522f',
                
                'APN' : '#c15355',
                'MRN' : '#984445'
                }
        
    elif name == 'Allen CCF' :
        
        palette = {'VISp' : '#08858C',
                'VISl' : '#08858C',
                'VISal' : '#08858C',
                'VISrl' : '#009FAC',
                'VISpm' : '#08858C',
                'VISam' : '#08858C',
                
                'DG' : '#7ED04B',
                'CA3' : '#7ED04B',
                'CA1' : '#7ED04B',
                'CA' : '#7ED04B',
                'POST' : '#48C83C',
                'SUB' : '#4FC244',
                
                'LGd' : '#FF8084',
                'LP' : '#FF909F',
                'LD' : '#FF909F',
                
                'APN' : '#FF90FF',
                'MRN' : '#FF90FF'
                }
        
    elif name == 'Rainbow' :
        
        palette = {'VISp' : '#F6BB42',
                'VISl' : '#37BC9B',
                'VISal' : '#967ADC',
                'VISrl' : '#4A89DC',
                'VISpm' : '#E9573F',
                'VISam' : '#DA4453',
                
                'DG' : '#37BC9B',
                'CA3' : '#37BC9B',
                'CA1' : '#37BC9B',
                'CA' : '#7ED04B',
                'POST' : '#48CFAD',
                'SUB' : '#48CFAD',
                
                'LGd' : '#D770AD',
                'LP' : '#EC87C0',
                'LD' : '#EC87C0',
                
                'APN' : '#434A54',
                'MRN' : '#656D78'
                }
        
    elif name == 'cmocean':
                
        import cmocean
        
        hierarchy_colors = cmocean.cm.phase(np.arange(1.0,0.1,-0.124))
            
        palette = {
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
        
    elif name == 'seaborn':

        colors = [[217,141,194],
                  [129,116,177],
                  [78,115,174],
                  [101,178,201],
                  [88,167,106],
                  [202,183,120],
                  [219,132,87],
                  [194,79,84]]
        
        def scale_colors(color):
            return [col/255 for col in color]
        
        hierarchy_colors = [scale_colors(col) for col in colors]
        
        palette = {
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
        
    else:
        raise Error('No matching palette name')
        
    if area in palette.keys():
        
        return palette[area]
    
    else:
        return default_color
    
    