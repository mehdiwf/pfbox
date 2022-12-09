import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def reformat_dataframe(df):
    """if there is a '#' in the headers line, shifts the columns names to
    their right place (this deletes the two last columns of the
    dataframe, that are in general just NaN values)
    """
    columns_name = df.columns[1:]
    last_columns_name = df.columns[-1]
    df.drop(columns = last_columns_name, inplace=True)
    df.columns = columns_name
    return None

def get_ordered_filelist(directory, prefix):
    file_list = os.listdir(directory)
    files_list = list(filter(lambda x: x.startswith(prefix), file_list))
    
    # sort the list by the number after "profile_"
    nb = len(prefix)
    files_list.sort(key=lambda x: int(x[nb:]))
    # files_list_abspath = [f"{directory}/{fname}" for fname in files_list]
    return files_list

def concatenate_all_data(filenames, time_list):
    """time_list is the time associated with each filename in filenames
    (filename is the fullpath)
    """
    df = pd.DataFrame()
    for i, fname in enumerate(filenames):
        df_temp = pd.read_csv(fname, delimiter = " ", header=0)
        reformat_dataframe(df_temp)
        df_temp['time'] = time_list[i]
        df = pd.concat((df, df_temp))
    return df

# --------------------------
def get_simulation_data(data_directory, prefix = "step_"):
    """returns a df with all the profile data in it, and the
    profile_index list (time list)
    """
    profile_files = get_ordered_filelist(data_directory, prefix)
    profile_files_fullpath = [f"{data_directory}/{fname}" for fname in profile_files]
    profile_index = [int(fname[len(prefix):]) for fname in profile_files]
    df = concatenate_all_data(profile_files_fullpath, profile_index)
    return df, profile_index
# --------------------------

def extract_pfbox_sim_info(simulation_dir, prefix = "step_"):
    """returns a dictionnary with a df containing all the profiles data of
    the simulation, a list containing all the time accessibles for
    profiles data, and a dictionnary containing the readme info of the
    simulation:

    keys:
    'df'
    'profile_time_list'
    'readme'

    """    
    simu_info = {}
    readmefile = f"{simulation_dir}/README.txt"
    df, profile_index = get_simulation_data(f"{simulation_dir}/sim_data", prefix= prefix)
    simu_info['df'] = df
    simu_info['profile_time_list'] = profile_index
    # do something with readme :toho:
    return simu_info
    
def extract_simulation_info(simulation_dir,
                            readme = True, # unused :toho:
                            prefix = "profile_"
                            ):
    """returns a dictionnary with a df containing all the profiles data of
    the simulation, a list containing all the time accessibles for
    profiles data, and a dictionnary containing the readme info of the
    simulation:

    keys:
    'df'
    'profile_time_list'
    'readme'

    """
    simu_info = {}
    readmefile = f"{simulation_dir}/README"
    df, profile_index = get_simulation_data(simulation_dir, prefix= prefix)
    simu_info['df'] = df
    simu_info['profile_time_list'] = profile_index
    # if readme:
        # do something with readme :toho:
    return simu_info

def get_convergence_data(simulation_dictionnary, column):
    """using the output of extract_simulation_info function returns x and
    y of convergence data of the simulation for a particular df column

    the convergence data is the sum of the absolute value of the derivative for each index
    """
    differences_list = []
    time_list = []
    simu_time_index = simulation_dictionnary['profile_time_list']
    df = simulation_dictionnary['df']
    dt = simulation_dictionnary['readme']['dt']
    for i in range(len(simu_time_index)-1):
        t = simu_time_index[i]
        t_plus_dt = simu_time_index[i+1]
        df_t = df[df['time'] == t]
        df_t_plus_dt = df[df['time'] == t_plus_dt]
        # calculations of the derivative
        derivatives = (df_t[column] - df_t_plus_dt[column])/dt
        differences = np.sum(np.abs(derivatives))
        differences_list.append(differences)
        time_list.append(t)
    return time_list, differences_list

def get_directory_list(searching_directory, beginning_dir_name, sort=True):
    """get the list of files/directories beginning with
    'beginning_dir_name', and order them if they are in the form
    '{beginning_dir_name}XX' with XX being a number
    """
    files = os.listdir(searching_directory)
    data_dir_list = list(filter(lambda x: x.startswith(beginning_dir_name),
                                files))
    if sort:
        data_dir_list.sort(key=lambda x: int(x[len(beginning_dir_name):]))
    datadir_list_abspath = [f"{searching_directory}/{fname}" for fname in data_dir_list]
    return datadir_list_abspath
    
def get_df_at_percent_time(simulation_dic, percent):
    """
    percent = 100 -> get last time of simulation
    """
    df = simulation_dic['df']
    profile_time_index = simulation_dic['profile_time_list']
    time_array = np.array(profile_time_index)
    max_time = time_array[-1]
    corresponding_time = (percent/100)*max_time
    closest_time_arg = np.argmin(np.abs(time_array - corresponding_time))
    df_time = df[df['time'] == profile_time_index[closest_time_arg]]
    return df_time

def show_plot_evolution(simu_dic, column, interval=1, save=False,
                        space_index_column = 'column'):
    """
    save is the str of the pdf file if you want to save your plot
    """
    
    time_index_to_plot = simu_dic['profile_time_list'][::interval]
    df = simu_dic['df']
    
    gradient = np.linspace(0, 1, len(time_index_to_plot))
    # from red to blue
    color_tuples = [(1-c, 0., c) for c in gradient]
    
    # plotting the evolution of a parameter
    for i, time in enumerate(time_index_to_plot[1:-1]):
        df_time = df[df['time'] == time]
        plt.plot(df_time[space_index_column],
                 df_time[column], color = color_tuples[i])
    init_time = time_index_to_plot[0]
    df_time_0 = df[df['time'] == init_time]
    final_time = time_index_to_plot[-1]
    df_time_final = df[df['time'] == final_time]
    plt.plot(df_time_0[space_index_column],
             df_time_0[column], color = color_tuples[0],
             label=f"time = {init_time}")
    plt.plot(df_time_final[space_index_column],
             df_time_final[column], color = color_tuples[-1],
             label=f"time = {final_time}")
    plt.title(f"evolution of {column} for different times")
    plt.xlabel('column index')
    plt.ylabel(column)
    plt.legend()
    if save:
        plt.savefig(save)
    plt.show()
