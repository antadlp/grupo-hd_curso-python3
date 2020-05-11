import pandas as pd
import numpy as np
from mendeleev import element
from fuerzas import *

def ranf():
    return np.random.uniform()

def conf_inicial_2d(**kwargs):
    
    N = kwargs.get('N', 20)
    dens = kwargs.get('dens', 1)
    iter_check = kwargs.get('iter_check', 10)
    tol = kwargs.get('tol', 0.7)
    seed = kwargs.get('seed', 0)
    path_out = kwargs.get('path_out', 'conf_inicial.pkl')
    elements = kwargs.get('elements', ['He'])
    fsize = kwargs.get('fsize', 250)
   
    np.random.seed(seed)
    elms = np.random.choice(elements, size=N)
    
    mu_e =  kwargs.get('mu_e', 0.7)
    sigma_e = kwargs.get('sigma_e', 0.3)
    dist_e = np.abs(np.random.normal(mu_e, sigma_e, N))
    charges = kwargs.get('charges', dist_e)

    box = np.power(np.divide(N, dens), np.divide(1, 2))
    box2 = box/2

    r = {}
    r['x'] = {}
    r['y'] = {}
    r['vx'] = {}
    r['vy'] = {}
    r['element'] = {}
    r['m'] = {}
    r['e'] = {}
    r['sizes'] = {}
    
    for i in range(N):
        r['x'][i] = box*ranf() - box2
        r['y'][i] = box*ranf() - box2
        r['vx'][i] = np.random.randn()
        r['vy'][i] = np.random.randn()
        el = elms[i]
        m = element(el).atomic_weight        
        r['element'][i] = el
        r['m'][i] = m
        r['e'][i] = charges[i]
        r['sizes'][i] = fsize*m


    for k in range(iter_check):
        dr = 0
        for i in range(N - 1):
            rxi = r['x'][i]
            ryi = r['y'][i]
            jlist = np.arange(i+1, N, 1)
            for j in jlist:
                rxij = rxi - r['x'][j]
                ryij = ryi - r['y'][j]
                rsij = rxij*rxij+ryij*ryij
                flag_rsij = False
                if (rsij < tol):
                    flag_rsij = True
                while(flag_rsij):
                    r['x'][j] = box*ranf() - box2
                    r['y'][j] = box*ranf() - box2
                    dr=dr+1
                    rxij = rxi - r['x'][j]
                    ryij = ryi - r['y'][j]
                    rsij = rxij*rxij+ryij*ryij

                    if (rsij >= 0.7):
                        flag_rsij = False
            
    dfr = pd.DataFrame(r)

    center = dfr[['x', 'y']].mean().values
    mv = np.array([box/2, box/2]) - center
    rdic = {}
    rdic['x'] = {}
    rdic['y'] = {}
    for idx in dfr.index:
        rdic['x'][idx] = dfr['x'].loc[idx] + mv[0]
        rdic['y'][idx] = dfr['y'].loc[idx] + mv[1]

    dfr.update(pd.DataFrame(rdic))
    dfr = dfr.join(get_coulomb2d(dfr))
   
    dic = {}
    dic['df'] = dfr
    dic['dr'] = dr
    dic['box'] = box
    dic['path_out'] = path_out
    
    dfr.to_pickle(path_out)
    
    return dic


