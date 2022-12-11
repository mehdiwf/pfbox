from simplot import *
import pfbox_simplot as pfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

c_dir = "../../water_c_code/looped_vdw/compar/2D_NX100_simu_VdW_compar_0"
rust_dir = "../src/testoutput/defaultdir_c_comparison"
# rust_dir = "../src/testoutput/defaultdir"

simdic_c = extract_simulation_info(c_dir)
simdic_rust  = pfp.extract_pfbox_sim_info(rust_dir)

dfc0 = get_df_at_percent_time(simdic_c, 0)
dfc_f = get_df_at_percent_time(simdic_c, 100)
dfr_f = pfp.get_df_at_percent_time(simdic_rust, 100)
dfr = simdic_rust['df']
dfc = simdic_c['df']

dfr_mf = dfr[dfr['time'] == 8999]

diff = dfc_f['vx_profile'] - dfr_mf["v_x"]

# pfp.show_plot_evolution(simdic_rust, "density")
pfp.show_plot_evolution(simdic_rust, "v_x")

df_vx = dfr[dfr['column'] == 25]
plt.plot(df_vx['time'], df_vx['v_x'], label="rust")

df_vx = dfc[dfc['j'] == 25]
plt.plot(df_vx['time'], df_vx['vx_profile'], label="c")
plt.legend()
plt.show()

# show_plot_evolution(simdic_c, "vx_profile")
# pfp.show_plot_evolution(simdic_rust, "v_x")

# plt.plot(dfr_f['column'], diff)
# plt.show()
# print(dfc['time'].max())
# print(dfr['time'].max())

# pfp.show_plot_evolution(simdic_rust, "density")
