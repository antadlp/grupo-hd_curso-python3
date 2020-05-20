import pandas as pd
import numpy as np
from mendeleev import element
import sys
sys.path.append("/home/antadlp/Documents/grupo-hd_curso-python3/consola2")
import subprocess


def ranf():
    return np.random.uniform()

def conf_inicial_md(**kwargs):
    
    N = kwargs.get('N', 20)
    dens = kwargs.get('dens', 0.6)
    iter_check = kwargs.get('iter_check', 50)
    tol = kwargs.get('tol', 0.7)
    seed = kwargs.get('seed', 0)
    path_out = kwargs.get('path_out', 'conf_inicial.pkl')
    elements = kwargs.get('elements', ['He'])
   
    np.random.seed(seed)
    elms = np.random.choice(elements, size=N)
    
    box = np.power(np.divide(N, dens), np.divide(1, 2))
    box2 = box/2

    r = {}
    r['x'] = {}
    r['y'] = {}
    r['z'] = {}
    r['vx'] = {}
    r['vy'] = {}
    r['vz'] = {}
    r['element'] = {}
    r['m'] = {}
    
    for i in range(N):
        r['x'][i] = box*ranf() - box2
        r['y'][i] = box*ranf() - box2
        r['z'][i] = box*ranf() - box2
        el = elms[i]
        m = element(el).atomic_weight        
        r['element'][i] = el
        r['m'][i] = m
    
    contador_iter_check = 1
    while(contador_iter_check <= iter_check):
        # se acostumbra while cuando NO CONOCES EL NUMERO
        # DE PASOS en ciclos de convergencia
        #
        # Se acostumbra el for cuando SI CONOCES EL
        # numero de pasos
#    for k in range(iter_check):
        dr = 0
        for i in range(N - 1):
            rxi = r['x'][i]
            ryi = r['y'][i]
            rzi = r['z'][i]
            jlist = np.arange(i+1, N, 1)
            for j in jlist:
                rxij = rxi - r['x'][j]
                ryij = ryi - r['y'][j]
                rzij = rzi - r['z'][j]
                rsij = rxij*rxij+ryij*ryij+rzij*rzij
                flag_rsij = False
                if (rsij < tol):
                    flag_rsij = True
                while(flag_rsij):
                    r['x'][j] = box*ranf() - box2
                    r['y'][j] = box*ranf() - box2
                    r['z'][j] = box*ranf() - box2
                    dr=dr+1
                    rxij = rxi - r['x'][j]
                    ryij = ryi - r['y'][j]
                    rzij = rzi - r['z'][j]
                    rsij = rxij*rxij+ryij*ryij+rzij*rzij

                    if (rsij >= 0.7):
                        flag_rsij = False
        contador_iter_check+=1
            
    dfr = pd.DataFrame(r)
  
    dic = {}
    dic['df'] = dfr
    dic['dr'] = dr
    dic['box'] = box
    dic['path_out'] = path_out
    
    dfr.to_pickle(path_out)
    
    return dic
    
    

def crear_mdorg(**kwargs):
    
    PATH_COORDS = "coords.dat"
    PATH_COORDS = kwargs.get("path_coords", PATH_COORDS)
    
    PATH_ADDVEL_BIN = "./addvel.exe"
    PATH_ADDVEL_BIN = kwargs.get("path_addvel_bin", PATH_ADDVEL_BIN)

    MDORG_FILE = "mdorg.dat"
    MDORG_FILE = kwargs.get("mdorg_file", MDORG_FILE)
    
    INPUT_VEL_NAME = "addvelX.inp"
    INPUT_VEL_NAME = kwargs.get("input_vel_name", INPUT_VEL_NAME)
    
    NUM_PARTICULAS = 33
    NUM_PARTICULAS = kwargs.get("num_particulas", NUM_PARTICULAS)
    
    NON_FROZEN = NUM_PARTICULAS
    NON_FROZEN = kwargs.get("non_frozen", NON_FROZEN)
    
    TEMPERATURA = 0.8
    TEMPERATURA = kwargs.get("temperatura", TEMPERATURA)
    
    DENS = 0.3 
    DENS = kwargs.get("dens", DENS)

    dic = conf_inicial_md(N=NUM_PARTICULAS, dens=DENS, iter_check=15)
    df = dic['df']
    L = dic['box']
    L = np.round(L, decimals=5)

    f = open(PATH_COORDS, 'w')
    f.write("{}\n".format(len(df)))
    f.write("{:12.5f}{:12.5f}{:12.5f}\n".format(L, L, L))

    for idx in df.index:

        x = df['x'].loc[idx]
        y = df['y'].loc[idx]
        z = df['z'].loc[idx]

        f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(x, y, z))

    f.close()

    f = open(INPUT_VEL_NAME, 'w')
    f.write("{}\n{}\n{}\n{}\n{}".format(PATH_COORDS,\
                                        MDORG_FILE,\
                                        NUM_PARTICULAS,\
                                        NON_FROZEN,\
                                        TEMPERATURA))

    f.close()

    cmd = [PATH_ADDVEL_BIN, '<', INPUT_VEL_NAME]
    cmd = [str(item) for item in cmd]
    cmd = ' '.join(cmd)
    process = subprocess.run(cmd, check=True,\
                                shell=True,\
                        stdout=subprocess.PIPE,\
                        universal_newlines=True)
    
    
    return
    
    
