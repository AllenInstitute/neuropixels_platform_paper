# -*- coding: utf-8 -*-
"""
Created on Fri Aug 9 10:23:56 2019

@author: Xiaoxuan Jia
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter


def fit_exp(rsc_time_matrix, color):
    
    intr = np.abs(rsc_time_matrix)
    

    tmp = np.nanmean(intr, axis=0)
    n=intr.shape[0]
    
    t = np.arange(len(tmp))[1:]
    y=gaussian_filter(np.nanmean(tmp, axis=0)[1:],0.78)
    
    p, amo = curve_fit(lambda t,a,b,c: a*np.exp(-1/b*t)+c,  t,  y,  p0=(-4, 2, 1), maxfev = 1000000000)

    a=p[0]
    b=p[1] # this is the intrinsic timescale
    c=p[2]
    y_std = np.nanstd(tmp, axis=0)[1:]/np.sqrt(n)

    return t, y, y_std, a, b, c