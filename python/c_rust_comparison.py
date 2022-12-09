from simplot import *
import pfbox_simplot as pfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

c_dir = "../../water_c_code/looped_vdw/compar/2D_NX100_simu_VdW_compar_0/"
rust_dir = "../src/testoutput/defaultdir_c_comparison/"

simdic_c = extract_simulation_info(c_dir)
simdic_rust  = pfp.extract_pfbox_sim_info(rust_dir)

dfc0 = get_df_at_percent_time(simdic_c, 0)
dfc = get_df_at_percent_time(simdic_c, 100)
dfr = pfp.get_df_at_percent_time(simdic_rust, 100)

show_plot_evolution(simdic_c, "vx_profile")
pfp.show_plot_evolution(simdic_rust, "v_y")
pfp.show_plot_evolution(simdic_rust, "density")
