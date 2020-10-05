#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 13:27:45 2020

@author: joshs
"""

# %%%


import os

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d

from scipy.ndimage.filters import gaussian_filter1d
from scipy.stats import pearsonr
from scipy import signal
from scipy.optimize import curve_fit

from scipy.ndimage.filters import gaussian_filter1d

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.stimulus_analysis.flashes import Flashes

import warnings

# This code constructs the units data frame used for most of the figures

### PATH VARIABLES ###
cache_directory = '/mnt/nvme0/ecephys_cache_dir_2'
######################
# %%

manifest_path = os.path.join(cache_directory, "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

sessions = cache.get_session_table()

# %%

def calculate_autocorrelation_timescale(dataset, timespan = (0.04, 0.29), unit_id=None):
    
    t1 = timespan[0]
    t2 = timespan[1]
    
    if unit_id is None:
        spikes = dataset.sum(dim='unit_id').sel(
                            time_relative_to_stimulus_onset = slice(t1, t2)
                            ).data
    else:
        spikes = dataset.sel(
                              unit_id = units.index.values[unit_idx],
                              time_relative_to_stimulus_onset = 
                              slice(t1, t2)
                              ).data
        
    
    nbins = spikes.shape[1]

    ac_matrix = autocorr2D(spikes)

    ACCG = np.mean(ac_matrix, axis=0)
    ACCG = ACCG[nbins//2:]
    t = np.linspace(0,nbins/2*10,len(ACCG))

    params, error = fit_exp(t, ACCG)
    
    timescale = params[1]
    err = error[1]
    spike_count = np.sum(spikes)
    
    return timescale, err, spike_count
    

def calculate_intrinsic_timescale(dataset, timespan = (1.0, 2.0), unit_id=None):

    t1 = timespan[0]
    t2 = timespan[1]
    
    if unit_id is None:
        spikes = dataset.sum(dim='unit_id').sel(
                            time_relative_to_stimulus_onset = slice(t1, t2)
                            ).data
    else:
        spikes = dataset.sel(
                              unit_id = units.index.values[unit_idx],
                              time_relative_to_stimulus_onset = 
                              slice(t1, t2)
                              ).data
 
    nbins = spikes.shape[1]

    rsc_matrix, T = intrinsic_timescale(spikes)
  
    t = np.linspace(10,500,49)

    try:
        params, error = fit_exp(t, T[1:50])
    except ValueError:
        params = [np.nan, np.nan, np.nan]
        error = [np.nan, np.nan, np.nan]

    timescale = params[1]
    err = error[1]
    spike_count = np.sum(spikes)

    return timescale, err, spike_count

 
def intrinsic_timescale(data):

    nbins = data.shape[1]
    rsc_matrix = np.empty((nbins, nbins)) * 0

    for i in np.arange(nbins-1):
        for j in np.arange(i+1, nbins):
            r, p = pearsonr(data[:,i], data[:,j])
            rsc_matrix[i, j] = r
            
    T = np.zeros((nbins-1,))
    for i in range(nbins-1):
        T[i] = np.nanmean(np.diag(rsc_matrix, k=i+1))

    return rsc_matrix, T

def autocorr2D(x):
    
    ac_matrix = signal.correlate(x, x, mode='same')
    ac_matrix = np.delete(ac_matrix, 
                          [ac_matrix.shape[0] // 2], 
                          axis=0)
    
    return ac_matrix

exponential = lambda t, a, b, c: a * np.exp(-1 / b * t) + c
    
def fit_exp(t, y):
    
    params, pcov = curve_fit(exponential, 
                       t, y, p0 = (5, 20, 0.1), method ='trf',
                       bounds = ([0, 1, -np.inf], [np.inf, 1000, np.inf]),
                       maxfev = 1000000000)
    
    error = np.sqrt(np.diag(pcov))
    
    return params, error

def plot_exp(t, params):


    y = exponential(t, *params)

    plt.plot(t, y, '-',)
    
def plot_raster(data):
    
    x, y = np.where(data)
    y = np.random.rand(len(y)) 
    plt.scatter(y,x,s=1)


# %%

areas = ['LGd','VISp','VISrl','VISl','VISal','LP','VISpm','VISam']

overall_df = []
units_df = []

for session_id in sessions.index.values:
    
    print('Session: ' + str(session_id))

    session = cache.get_session_data(session_id, amplitude_cutoff_maximum = np.inf,
                                              isi_violations_maximum = np.inf,
                                              presence_ratio_minium = -np.inf)
    
    fl = Flashes(session)
    
    stim_table = fl.stim_table

    for area in areas:
        
        units = session.units[(session.units.ecephys_structure_acronym == area) & 
                              (session.units.p_value_rf < 0.01) & 
                              (session.units.area_rf < 2500) & 
                              (session.units.snr > 1) &
                              (session.units.firing_rate_dg > 0.1)]
        
        if len(units) > 0:
            
            print(' ' + area + ': ' + str(len(units)) + ' units')
        
            dataset = session.presentationwise_spike_counts(
                    bin_edges=np.arange(0, 2.01, 0.01),
                    stimulus_presentation_ids = stim_table.index.values,
                    unit_ids=units.index.values
                )
            
            if len(units) > 15:
                
                timescale_ac, err_ac, spike_count_ac = calculate_autocorrelation_timescale(dataset)
                timescale_it, err_it, spike_count_it = calculate_intrinsic_timescale(dataset)
            
                        
                overall_df.append(pd.DataFrame(data = {'area' : [area],
                                                       'session' : [session_id],
                                                       'timescale_ac' : [timescale_ac],
                                                       'err_ac': [err_ac],
                                                       'spike_count_ac' : [spike_count_ac],
                                                       'timescale_it' : [timescale_it],
                                                       'err_it' : [err_it],
                                                       'spike_count_it' : [spike_count_it],
                                                       'unit_count' : [len(units)]}))
                
            timescale_units_ac = np.zeros((len(units),))
            err_units_ac = np.zeros((len(units),))
            spike_count_units_ac = np.zeros((len(units),))
            timescale_units_it = np.zeros((len(units),))
            err_units_it = np.zeros((len(units),))
            spike_count_units_it = np.zeros((len(units),))
            
            print('   Looping through units...')
                
            for unit_idx, unit_id in enumerate(units.index.values):
                
                ts_ac, err_ac, sc_ac = \
                    calculate_autocorrelation_timescale(dataset, unit_id=unit_id)
                    
                ts_it, err_it, sc_it = \
                    calculate_intrinsic_timescale(dataset, unit_id=unit_id)
                    
                timescale_units_ac[unit_idx] = ts_ac
                err_units_ac[unit_idx] = err_ac
                spike_count_units_ac[unit_idx] = sc_ac
                timescale_units_it[unit_idx] = ts_it
                err_units_it[unit_idx] = err_it
                spike_count_units_it[unit_idx] = sc_it
                
            units_df.append(pd.DataFrame(data = {'area' : [area] * len(units),
                                                 'session' : [session_id] * len(units),
                                                 'unit_id' : units.index.values,
                                                 'timescale_ac' : timescale_units_ac,
                                                 'err_ac' : err_units_ac,
                                                 'spike_count_ac' : spike_count_units_ac,
                                                 'timescale_it' : timescale_units_it,
                                                 'err_it' : err_units_it,
                                                 'spike_count_it' : spike_count_units_it
                                                 }))
                
# %%

df = pd.concat(units_df)

# %%

df.to_csv('/mnt/nvme0/unit_tables/unit_table_autocorr_200327.csv')

# %%

df = pd.concat(overall_df)

# %%%

df.to_csv('/mnt/nvme0/unit_tables/area_table_autocorr_200327.csv')

# %%

df = pd.read_csv('/mnt/nvme0/unit_tables/unit_table_autocorr_200327.csv',
                 index_col =0)

# %%

from scipy import stats
from scipy.ndimage.filters import gaussian_filter1d

def get_bootstrap_95ci(M, measure_of_central_tendency, N=1000):
    n = int(len(M)/2)
    est = np.zeros((N,))
    for i in range(N):
        boot = M[np.random.permutation(len(M))[:n]]
        est[i] = measure_of_central_tendency(boot)
        
    return np.percentile(est,97.5) - np.nanmean(est)

def do_not_change(original_value):
    return original_value

measure_of_central_tendency = np.nanmean

areas = ('LGd','VISp','VISl','VISrl','LP','VISal','VISpm','VISam')

color_palette= 'seaborn'

HS = [-0.515, -0.357, -0.093, -0.059, 0.105, 0.152,0.327, 0.441]


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'font.size': 13})

# %%

from scipy.stats import linregress, spearmanr, pearsonr

metric = 'timescale_ac'

total_units = 0

plt.figure(19191)
plt.clf()

plt.subplot(3,2,1)

centers = np.zeros((8,))
errorbars = np.zeros((8,))

for area_idx, area in enumerate(areas):
    
    sub_df = df[(df.area == area) &
                (df[metric] < 300) &
                (df[metric] > 1) & 
                (df.spike_count_ac > 50) & 
                (df.err_ac < 20)]

    y = sub_df[metric].values
    
    centers[area_idx] = measure_of_central_tendency(y)
    errorbars[area_idx] = get_bootstrap_95ci(y, measure_of_central_tendency)
    
    h, b = np.histogram(y, bins=np.linspace(0,120,40), density=True)

    h_filt = gaussian_filter1d(h,1.5)
    
    max_value = np.max([np.max(h_filt), 0.025])

    plt.plot(b[:-1],h_filt,color=get_color_palette(areas[area_idx], color_palette))
    plt.xlabel('Response decay timescale')
    
    total_units += len(sub_df)
    
plt.ylim([0, max_value * 1.1])
#plt.gca().get_yaxis().set_visible(False)
ax = plt.gca()
[ax.spines[loc].set_visible(False) for loc in ['right', 'top']]        

plt.subplot(3,2,2)

x = HS
y = centers

for k in range(8):
    plt.plot(x[k], y[k],'.',color=get_color_palette(areas[k], color_palette))
    plt.errorbar(x[k], y[k], yerr = errorbars[k], fmt='.',color=get_color_palette(areas[k], color_palette),alpha=0.8)

slope,intercept,r,p,std = linregress(x,y)
x2 = np.linspace(-0.75,0.5,10)

plt.plot(x2,x2*slope+intercept,'--k', alpha=0.5)

r_s,p_s = spearmanr(x,y)
r_p,p_p = pearsonr(x,y)

text =  '$r_P$ = ' + str(np.around(pow(r_p,1),2)) + '; $P_P$ = ' + str(np.around(p_p,4)) + '\n' + \
        '$r_S$ = ' + str(np.around(pow(r_s,1),2)) + '; $P_S$ = ' + str(np.around(p_s,4))

plt.text(-0.0,36,text,fontsize=10)
plt.ylabel('Response decay timescale')
plt.xlabel('Anatomical hierarchy score')
ax = plt.gca()
[ax.spines[loc].set_visible(False) for loc in ['right', 'top']] 


# %%
metric = 'timescale_it'

total_units_it = 0

plt.subplot(3,2,3)

centers = np.zeros((8,))
errorbars = np.zeros((8,))

for area_idx, area in enumerate(areas):
    
    sub_df = df[(df.area == area) &
                (df[metric] < 400) &
                (df[metric] > 5) & 
                (df.spike_count_it > 100) & 
                (df.err_it < 100)]
    y = sub_df[metric]
    
    centers[area_idx] = measure_of_central_tendency(y)
    errorbars[area_idx] = get_bootstrap_95ci(y, measure_of_central_tendency)
    
    h, b = np.histogram(y, bins=np.linspace(0,200,40), density=True)

    h_filt = gaussian_filter1d(h,1.5)
    
    max_value = np.max(h_filt)

    plt.plot(b[:-1],h_filt,color=get_color_palette(areas[area_idx], color_palette))
    plt.xlabel('Intrinsic timescale')
    
    total_units_it += len(sub_df)
    
plt.ylim([0, max_value * 1.4])
#plt.gca().get_yaxis().set_visible(False)
ax = plt.gca()
[ax.spines[loc].set_visible(False) for loc in ['top', 'right']]        

plt.subplot(3,2,4)

x = HS
y = centers

for k in range(8):
    plt.plot(x[k], y[k],'.',color=get_color_palette(areas[k], color_palette))
    plt.errorbar(x[k], y[k], yerr = errorbars[k], fmt='.',color=get_color_palette(areas[k], color_palette),alpha=0.8)

slope,intercept,r,p,std = linregress(x,y)
x2 = np.linspace(-0.75,0.5,10)

plt.plot(x2,x2*slope+intercept,'--k', alpha=0.5)

r_s,p_s = spearmanr(x,y)
r_p,p_p = pearsonr(x,y)

text =  '$r_P$ = ' + str(np.around(pow(r_p,1),2)) + '; $P_P$ = ' + str(np.around(p_p,3)) + '\n' + \
        '$r_S$ = ' + str(np.around(pow(r_s,1),2)) + '; $P_S$ = ' + str(np.around(p_s,3))

plt.text(-0.4,60,text,fontsize=10)
plt.ylim([40,100])
plt.ylabel('Intrinsic timescale')
plt.xlabel('Anatomical hierarchy score')
ax = plt.gca()
[ax.spines[loc].set_visible(False) for loc in ['right', 'top']] 



# %%







    