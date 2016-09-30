import datetime
from itertools import chain
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from scipy import stats
import seaborn as sns

def read_data(filepath):
    
    """Read the multicultivator CSV data, formats NaN values, and returns a dataframe.
    
    This function has not much functionality, however, it simplifies 
    the process of reading the data into a functional dataframe.
    For now this is all it needs to do. If we encounter any variations
    in the default file, e.g. the note section spans mulitple rows, 
    this function needs to be adapted.
    """
         
    df = pd.read_csv(filepath, skiprows=4, index_col=False, na_values='Overflow')
    
    return df.fillna(0)

def read_layout(filepath):
    
    """Read the multicultivator CSV data and parses the note section.
    
    This function reads the multicultivator CSV file line by line,
    extracts the content of each vessel from the note section, 
    and returns a list with potential columns names for the dataframe 
    and a list of the vessel contents.
    """
    
    with open(filepath) as file:
        read = file.readlines()
    
    layout = read[1]
    layout = layout.split(':',1)[1][:-1]
    layout = re.split('[._*-,]', layout)
    
        
    names = chain.from_iterable((name+'_680',
                                 name+'_720')
                                 for name in layout)
    
    layout = chain.from_iterable((name,name) for name in layout)
    
    return ['Time','Temp']+list(names), list(layout)

def find_groups(lst):
    
    loc_dict = {}

    for strain in set(lst):
        count = 0
        index = []
        for i in range(lst.count(strain)):
            count = lst.index(strain, count)
            index.append(count+2)
            count += 1
        loc_dict[strain] = index
        
    return loc_dict

def format_date_stamp(file):
    
    with open(file) as f:
        read = f.readlines()
    
    date = read[3]
    date = date.split(':',1)[1][:-1]
    
    day = date.split(', ')[1]
    time = date.split(', ')[0]
    
    stamp = day + time
    timestamp = datetime.datetime.strptime(stamp, '%d %b %Y%H:%M:%S')
    
    return timestamp

def create_date(column, file, option=False):
    
    if option:
        column = pd.to_timedelta(column, unit='s')
        column = format_date_stamp(file) + column
    
    return column

def concat_measurements(df, *args):
    
    joined_df = df.iloc[:,:-1]
    
    for df in args:
        joined_df = pd.concat([joined_df, df.iloc[:,2:-1]], axis=1)
        
    return joined_df

def plot_raw(df, labels, time_unit='day'):
    
    time = {'second': [1, 50000], 'minute': [60, 1000],
            'hour': [60*60, 50], 'day': [60*60*24, 2]}
    
    fig, ax = plt.subplots(len(df.columns[2:-1])//8, len(labels)//(len(df.columns[2:-1])//8),
                           sharey='row', sharex='col', figsize=(18,4*len(df.columns[2:-1])//16))
    sns.despine()
    
    max_y_680 = np.around(df.iloc[:,2::2].max(axis=1).max(), decimals=2)
    max_y_720 = np.around(df.iloc[:,3::2].max(axis=1).max(), decimals=2)
    
    max_x = np.around(df['Time'].max()/time[time_unit][0])+time[time_unit][1]

    
    for sample, axis, label in zip(range(2,len(df.columns[:-1]),2),
                                   ax[::2].reshape(1,len(labels)//2)[0],
                                   labels[::2]):
        axis.plot(df['Time']/time[time_unit][0], df.iloc[:,[sample]],
                     color='skyblue')
        axis.set_title(label)
        axis.set_ylim(0, max_y_680)
        axis.set_xlim(0, max_x)
        axis.set_xticks([])
    
    for sample, axis in zip(range(2,len(df.columns[:-1]),2),
                            ax[1::2].reshape(1,len(labels)//2)[0]):
        axis.plot(df['Time']/time[time_unit][0], df.iloc[:,[sample+1]],
                     color='mediumseagreen')
        axis.set_ylim(0, max_y_720)
        axis.set_xlim(0, max_x)
        axis.set_xticks([])

    for line in range(0,len(df.columns[2:-1])//8,2):
        ax[line][7].text(max_x, max_y_680/2, '680 nm', fontsize=15)
        ax[line+1][7].text(max_x, max_y_720/2, '720 nm', fontsize=15)

    plt.savefig('raw.pdf', bbox_inches='tight')
    
def plot_groups(df, labels, time_unit='day', order=None, ci=0.95):
    
    groups = find_groups(labels)
    
    time = {'second': [1, 50000], 'minute': [60, 1000],
            'hour': [60*60, 50], 'day': [60*60*24, 2]}
    
    if not order:
        order = groups.keys()
    
    fig, ax = plt.subplots(2, len(groups), sharey='row', sharex='col', figsize=(10,5))
    sns.despine()

    max_y_680 = np.around(df.iloc[:,2::2].max(axis=1).max(), decimals=2)
    max_y_720 = np.around(df.iloc[:,3::2].max(axis=1).max(), decimals=2)
    
    max_x = np.around(df['Time'].max()/time[time_unit][0])+time[time_unit][1]

    for group, axis in zip(order, ax.reshape(6,1)[:3]):

        mean_680 = df.iloc[:,groups[group][::2]].mean(axis=1)
        std_680 = df.iloc[:,groups[group][::2]].std(axis=1)
        
        conf_int = stats.t.interval(ci, len(df)-1, loc=mean_680, scale=std_680/np.sqrt(len(df)))

        axis[0].plot(df['Time']/time[time_unit][0], mean_680, color='skyblue', linewidth=1)
        axis[0].fill_between(df['Time']/time[time_unit][0], conf_int[0], conf_int[1],
                             color='skyblue', alpha=0.6, linewidth=0)
        axis[0].set_ylim(0, max_y_680)
        axis[0].set_xlim(0, max_x)
        axis[0].set_title(group)

    ax[0][2].text(max_x, max_y_680/2, '680 nm', fontsize=12)

    for group, axis in zip(order, ax.reshape(6,1)[3:]):

        mean_720 = df.iloc[:,groups[group][1::2]].mean(axis=1)
        std_720 = df.iloc[:,groups[group][1::2]].std(axis=1)
        
        conf_int = stats.t.interval(ci, len(df)-1, loc=mean_720, scale=std_680/np.sqrt(len(df)))

        axis[0].plot(df['Time']/time[time_unit][0], mean_720, color='mediumseagreen', linewidth=1)
        axis[0].fill_between(df['Time']/time[time_unit][0], conf_int[0], conf_int[1],
                             color='mediumseagreen', alpha=0.6, linewidth=0)
        axis[0].set_ylim(0, max_y_720)
        axis[0].set_xlim(0, max_x)

    ax[1][2].text(max_x, max_y_720/2, '720 nm', fontsize=12)
    
    plt.savefig('group.pdf', bbox_inches='tight')
    
def growth_stats(df, labels, group, lower_bound, upper_bound, time_max, wave_length='720', time_unit='day'):
    
    groups = find_groups(labels)
    
    time = {'second': [1, 50000], 'minute': [60, 1000],
            'hour': [60*60, 50], 'day': [60*60*24, 2]}
    
    if wave_length == '680':
        x = 0
    elif wave_length == '720':
        x = 1
    else:
        raise ValueError("%s not allowed, please use '680' or '720'" % wave_length)
    
    part = df[(df.iloc[:,groups[group][x::2]].mean(axis=1) >= lower_bound) & 
              (df.iloc[:,groups[group][x::2]].mean(axis=1) <= upper_bound) & 
              (df['Time']/time[time_unit][0] <= time_max)]

    slope, intercept, r_value, p_value, std_err = stats.linregress(part['Time'],
                                                                   part.iloc[:,groups[group][x::2]].mean(axis=1))
    
    print('Slope:', slope)
    print('y-axis intercept:', intercept)
    print('R value:', r_value)
    print('R**2:', r_value**2)
    
    return slope, intercept, r_value, p_value, std_err
