from simplot import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Température pression densité_coex_vapeur densité_coex_liquide densité_spinodale_vapeur densité_spinodale_liquide
thierry_cols = ["temp",
                "pressure",
                "dens_coex_gas",
                "dens_coex_liquid",
                "dens_spinodal_gas",
                "dens_spinodal_liquid"]

data_thierry_path = "/home/mehdi/workdir/dossiers/ilm/these/code_simulations/ressources/thierry/diagramme_VdW.dat"
# df_thierry = pd.read_csv(data_thierry_path, delimiter = ' ',
#                          header=None)
# df_thierry.columns = thierry_cols

data_dir_begin = "2D_NX100_simu_VdW_heatc_critic_"
# data_dir_list = get_directory_list("./critic_simu", data_dir_begin)
# data_dir_list = get_directory_list("../src/testoutput/", data_dir_begin)
# data_dir_test = data_dir_list[-1]

simdic = extract_simulation_info("../src/testoutput", prefix = "step_")

# to_plot = 'density_profile'
# to_plot = 'Pxx_profile'
# to_plot = 'Pyy_profile'
# to_plot = 'vx_profile'
to_plot = 'T_profile'
to_plot = 'density'

# gradient = np.linspace(0, 1, len(data_dir_list))
# # from red to blue
# color_tuples = [(c, 0., 1-c) for c in gradient]

initial_temp_list = []
gas_density_list = []
liquid_density_list = []
final_temp_list = []
final_temp_std_list = []

show_plot_evolution(simdic, to_plot, interval=1, save=False, space_index_column = 'column')

# for datadir in data_dir_list:
#     simdic = extract_simulation_info(datadir)
#     df_lasttime = get_df_at_percent_time(simdic, 100)
#     liquid_density = df_lasttime[df_lasttime['j'] == 25]['density_profile']
#     gas_density = df_lasttime[df_lasttime['j'] == 75]['density_profile']

#     final_temp_list.append(df_lasttime['T_profile'].mean())
#     final_temp_std_list.append(df_lasttime['T_profile'].std())
#     initial_temp_list.append(simdic['readme']['T0'])
#     gas_density_list.append(gas_density)
#     liquid_density_list.append(liquid_density)

# plt.scatter(df_thierry['dens_coex_gas'], df_thierry['temp'],
#             color='g',
#             marker='o',
#             label='Thierry VdW gas',
#             alpha = 0.1)
# plt.scatter(df_thierry['dens_coex_liquid'], df_thierry['temp'],
#             color='g',
#             marker='o',
#             label='Thierry VdW liquid',
#             alpha = 0.1)
# plt.errorbar(gas_density_list, final_temp_list,
#              yerr = final_temp_std_list,
#              marker = 'd',
#              ls='',
#              color='b',
#              label='gas with final T')
# plt.errorbar(liquid_density_list, final_temp_list,
#              yerr = final_temp_std_list,
#              marker = 'd',
#              ls='',
#              color='r',
#              label='liquid with final T')
# plt.xlabel('density')
# plt.ylabel('temperature')
# plt.legend()
# plt.savefig('critic_temperature_comparison.pdf')
# plt.show()

# ------------------- ALREADY DONE PLOTS

plot_convergence = False
plot_indiv_evolution = False
plot_last_time = False

if plot_indiv_evolution:
    show_plot_evolution(simdic, to_plot,
                        interval=1,
                        save=False)

if plot_convergence:
    for i_sim, data_dir in enumerate(data_dir_list):
        simdic_i = extract_simulation_info(data_dir)
        T0 = simdic_i['readme']['T0']
        time_list, differences_list = get_convergence_data(simdic_i, to_plot)
        plt.plot(time_list, differences_list,
                 label=f'simu {i_sim}, T = {T0}',
                 color = color_tuples[i_sim])
    plt.xlabel('time')
    plt.ylabel(f'{to_plot} convergence')
    plt.legend()
    # plt.savefig(figure_name)
    plt.show()        

# last time
if plot_last_time:
    for i_sim, data_dir in enumerate(data_dir_list):
        simu_dic = extract_simulation_info(data_dir)
        time_index_to_plot = simu_dic['profile_time_list'][:]
        T0 = simu_dic['readme']['T0']
        df_lasttime = get_df_at_percent_time(simu_dic, 100)
        
        # # ---------------- plotting last time 
        plt.plot(df_lasttime['j'],
                 df_lasttime[to_plot],
                 color=color_tuples[i_sim],
                 label=f'simu {i_sim}, T = {T0}')
    plt.xlabel('column index')
    plt.ylabel(to_plot)
    plt.legend()
    plt.savefig('test.pdf')
    plt.show()

