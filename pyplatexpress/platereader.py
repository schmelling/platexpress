import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import string

def read_data(filepath):
    
    '''Read the plate reader CSV data and returns a dataframe.
    
    This function has not much functionality, however, it simplifies 
    the process of reading the data into a functional dataframe.
    For now this is all it needs to do. If we encounter any variations
    in the default file this function needs to be adapted.
    '''
    
    df = pd.read_csv(filepath, skiprows=60)
    df = df.iloc[:241,1:-1]
    df['Time'] = pd.to_datetime(df['Time'])
    df[list(df.columns[1:])] = df[list(df.columns[1:])].apply(pd.to_numeric)
    
    return df

def plot_96well(df):
    
    max_y = np.around(plate.iloc[:,2:].max(axis=1).max(), decimals=2)
    
    fig, ax = plt.subplots(8, 12, sharex='col', sharey='row', figsize=(20,12))
    sns.despine()
    
    for column, axis in zip(df.columns[2:], ax.reshape(96, 1)):
        axis[0].plot(df['Time'], df[column])
        axis[0].set_ylim(0,max_y)
        axis[0].set_title(column)
        
        axis[0].tick_params(labelsize=8)
        axis[0].set_xticks([])
        
    fig.tight_layout()
        
    plt.savefig('raw_96well.pdf', bbox_inches='tight')
    
def read_layout(filepath, skipped=False):
    
    layout = pd.read_csv(filepath)
    layout = layout.fillna('empty')
    names = list(layout.iloc[:,1:].as_matrix().reshape(1,96)[0])
    
    if skipped:
        names = [value for value in names if value != 'skip']
        
    return ['Time','Temp']+list(names), list(names)

def remove_empty(df, labels):
    
    labels = [value for value in labels if value != 'empty']

    return df.drop('empty', axis=1), labels

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

def skip_wells(df, filepath, wells):

    layout = pd.read_csv(filepath)
    
    if isinstance(wells, str):
        
        row = string.ascii_uppercase.index(wells[0])
        layout.iloc[row,int(wells[1])] = 'skip'
        
    elif isinstance(wells, list):
        
        for well in wells:
            
            row = string.ascii_uppercase.index(well[0])
            layout.iloc[row,int(well[1])] = 'skip'
            
    layout.to_csv(filepath.split('.')[0]+'_skip.csv', index=False)
    
    return df.drop(wells, axis=1)

def plot_groups(df, labels, order=None, error='std', ci=0.95):
    
    groups = find_groups(labels)
    
    groups = {k: v for k, v in groups.items() if 'blank' not in k}
    
    if not order:
        order = groups.keys()
    
    fig, ax = plt.subplots(1, len(groups), sharey='row', sharex='col', figsize=(15,4.5))
    sns.despine()

    for group, axis in zip(order, ax.reshape(len(groups),1)):

        mean = df.iloc[:,groups[group]].mean(axis=1)
        std = df.iloc[:,groups[group]].std(axis=1)
        
        conf_int = stats.t.interval(ci, len(df)-1, loc=mean, scale=std/np.sqrt(len(df)))

        axis[0].plot(range(len(df['Time'])), mean, color='seagreen', linewidth=1)
        
        if error == 'std':
            axis[0].fill_between(range(len(df['Time'])), mean-std, mean+std,
                                 color='mediumseagreen', alpha=0.5, linewidth=0)
        elif error == 'ci':
            axis[0].fill_between(range(len(df['Time'])), conf_int[0], conf_int[1],
                                 color='mediumseagreen', alpha=0.5, linewidth=0)

        axis[0].set_title(group)
    
    plt.savefig('group_96well.pdf', bbox_inches='tight')
