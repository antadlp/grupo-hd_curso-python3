import numpy as np
import pandas as pd

def coulomb(qi, qj, ri, rj):
    f = 138.935
    rij = np.subtract(rj, ri)
    rij3 = np.power(np.linalg.norm(rij), 3)
    return f*np.divide(qi*qj, rij3)*rij
    
    
def multi_coulomb(qi, qj, ri, rj):
    
    f = 138.935
    rij = np.subtract(rj, ri)
    rij_dot = np.sum(np.multiply(rij, rij), axis=1)
    rij_norm = np.sqrt(rij_dot)
    rij3 = np.power(rij_norm, 3)
    
    _ = f*np.divide(np.multiply(qi, qj), rij3)
    
    return np.multiply(rij.T, _)
    
def get_coulomb2d(df): 

    dic = {}
    dic['fx'] = {}
    dic['fy'] = {}

    rjs = df[['x', 'y']].values
    qjs = df['e'].values
    fc = 0
    for idx in df.index:

        ri = np.array([df[['x', 'y']].loc[idx].values])
        qi = df['e'].loc[idx]

        rjs = df[['x', 'y']].drop(idx).values
        qjs = df['e'].drop(idx).values

        fc = multi_coulomb(qi, qjs, ri, rjs)

        dic['fx'][idx] = sum(fc[0])
        dic['fy'][idx] = sum(fc[1])


    return pd.DataFrame(dic)
    
    
    
    
    
    
    
    
    
    
    

