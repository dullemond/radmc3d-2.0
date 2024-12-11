import numpy as np
import pandas as pd

pd.set_option('display.max_rows', 1000)

def read_paths():
    df = pd.read_fwf('photon_paths.out',skiprows=2,header=None)
    df.columns  = ['iphot','ievent','x','y','z','nx','ny','nz']
    df['r']     = (df['x']**2+df['y']**2+df['z']**2)**0.5
    df['theta'] = np.arccos(df['z']/df['r'])
    rc          = (df['x']**2+df['y']**2)**0.5
    phi         = np.arctan2(df['y'],df['x'])
    phi[phi<0] += 2*np.pi
    df['phi']   = phi
    df.set_index(['iphot','ievent'],inplace=True,drop=True)
    return df
