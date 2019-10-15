from allensdk.brain_observatory.ecephys.ecephys_session import EcephysSession
from allensdk.brain_observatory.ecephys.stimulus_analysis.drifting_gratings import DriftingGratings 
from allensdk.brain_observatory.ecephys.stimulus_analysis.receptive_field_mapping import ReceptiveFieldMapping    
from allensdk.brain_observatory.ecephys.stimulus_analysis.flashes import Flashes    

cache_directory = "/mnt/hdd0/cache_dir_10_08/"

manifest_path = os.path.join(cache_directory, "manifest.json")

cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

# %%

sessions = cache.get_sessions()

session_id = 756029989

# %%

session = cache.get_session_data(session_id)

dg = DriftingGratings(session)
rf = ReceptiveFieldMapping(session)
fl = Flashes(session, psth_resolution=0.01)

# %%

v = session.units[session.units['manual_structure_acronym'] == 'VISp']

# %%

unit_id = 42
UNIT = v.index.values[unit_id]

plt.figure(1)
plt.clf()
fl.plot_conditionwise_raster(UNIT)

plt.figure(2)
plt.clf()
fl.plot_response(UNIT)

# %%
plt.figure(3)
plt.clf()
rf.plot_rf(UNIT)

plt.figure(4)
plt.clf()
rf.plot_conditionwise_raster(UNIT)

# %%

plt.figure(5)
plt.clf()
dg.make_star_plot(UNIT)
plt.ylim([-5,5])
plt.axis('equal')

plt.figure(6)
plt.clf()
dg.plot_conditionwise_raster(UNIT)

