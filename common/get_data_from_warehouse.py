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

# This code constructs the units data frame used for most of the figures

### PATH VARIABLES ###
cache_directory = '/mnt/hdd0/cache_dir_10_03'
code_directory = '/home/joshs/GitHub/neuropixels_platform_paper'
######################

manifest_path = os.path.join(cache_directory, "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

metrics = cache.get_unit_analysis_metrics_by_session_type('brain_observatory_1.1', filter_by_validity=False, amplitude_cutoff_maximum = np.inf,
                                                          presence_ratio_minimum = -np.inf,
                                                          isi_violations_maximum = np.inf)

metrics2 = cache.get_unit_analysis_metrics_by_session_type('functional_connectivity', filter_by_validity=False, amplitude_cutoff_maximum = np.inf,
                                                          presence_ratio_minimum = -np.inf,
                                                          isi_violations_maximum = np.inf)

all_metrics = pd.concat([metrics, metrics2])

ttf_df = pd.read_csv(os.path.join(code_direcotry, data, 'time_to_first_spike.csv'),index_col=0)

df = all_metrics.merge(ttf_df, on='ecephys_unit_id')
