# -*- coding: utf-8 -*-

import os

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import (
    EcephysSession, 
    removed_unused_stimulus_presentation_columns
)

# This code constructs the units data frame used for Figure 1, Figure 3,
# and Extended Data Figures 1, 4, and 10

### PATH VARIABLES ##############################
cache_directory = '/mnt/nvme0/ecephys_cache_dir_2'
code_directory = '/home/joshs/GitHub/neuropixels_platform_paper'
###################################################

# %%

# 1. Download pre-computed response metrics from AllenSDK

manifest_path = os.path.join(cache_directory, "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

metrics = cache.get_unit_analysis_metrics_by_session_type('brain_observatory_1.1', filter_by_validity=False, amplitude_cutoff_maximum = np.inf,
                                                          presence_ratio_minimum = -np.inf,
                                                          isi_violations_maximum = np.inf)

metrics2 = cache.get_unit_analysis_metrics_by_session_type('functional_connectivity', filter_by_validity=False, amplitude_cutoff_maximum = np.inf,
                                                          presence_ratio_minimum = -np.inf,
                                                          isi_violations_maximum = np.inf)

all_metrics = pd.concat([metrics, metrics2], sort=True)

# %% 

# 2. Add time to first spike (not included in data release)

first_spike = pd.read_csv(os.path.join(code_directory, 'data', 'time_to_first_spike.csv'),index_col=0)

all_metrics = all_metrics.merge(first_spike, on='ecephys_unit_id', sort=False)

# %%

# 3. Add layer info and cortical depth (derived from CCF annotations)

layer_info = pd.read_csv(os.path.join(code_directory, 'data', 'layer_info.csv'),index_col=0)

all_metrics = all_metrics.merge(layer_info, on='ecephys_unit_id', sort=False)

# %%

# 4. Add intrinsic timescale and response decay timescale

timescale_metrics = pd.read_csv(os.path.join(code_directory, 'data', 'timescale_metrics.csv'),
                                             index_col=0, low_memory=False)

all_metrics = all_metrics.merge(timescale_metrics, on='ecephys_unit_id', how='outer')

# %%

# 5. Write to a CSV file for use by other scripts

all_metrics.to_csv(os.path.join(code_directory, 'data', 'unit_table.csv'))

# %%


