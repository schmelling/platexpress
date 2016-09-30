import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

    for column, axis in zip(df.columns[2:], ax.reshape(96, 1)):
        axis[0].plot(df['Time'], df[column])
        axis[0].set_ylim(0,max_y)
        axis[0].set_title(column)
 
        axis[0].tick_params(labelsize=8)
        axis[0].set_xticks([])
        
    plt.savefig('raw_96well.pdf')
