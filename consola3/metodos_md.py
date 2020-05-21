import pandas as pd
import numpy as np
import subprocess
import sys


def ranf():
    return np.random.uniform()
    
def proper_round(num, dec=0):
    num = str(num)[:str(num).index('.')+dec+2]
    if num[-1]>='5':
      a = num[:-2-(not dec)]       # integer part
      b = int(num[-2-(not dec)])+1 # decimal part
      return float(a)+b**(-dec+1) if a and b == 10 else float(a+str(b))
    return float(num[:-1])
    

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
        m = 1.0        
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
    
    
def config_md(**kwargs):
    
    MDORG_FILE = kwargs.get("mdorg_file", "mdorg.dat")
    
    f = open(MDORG_FILE, 'r')
    
    N = int(f.readline().split()[0])
    
    L = f.readline().split()
    Lx = float(L[0])
    Ly = float(L[1])
    Lz = float(L[2])
    
    rx = []
    ry = []
    rz = []
    
    vx = []
    vy = []
    vz = []
    
    for line in f:
        
        sp = line.split()
        
        rx.append(sp[0])
        ry.append(sp[1])
        rz.append(sp[2])
        
        vx.append(sp[3])
        vy.append(sp[4])
        vz.append(sp[5])
    
    f.close()
    
    dic = {}
    dic['rx'] = np.array(rx).astype(np.float)
    dic['ry'] = np.array(ry).astype(np.float)
    _z = np.array(rz).astype(np.float)
    dic['rz'] = _z - Lz/2.0
    
    dic['vx'] = np.array(vx).astype(np.float)
    dic['vy'] = np.array(vy).astype(np.float)
    dic['vz'] = np.array(vz).astype(np.float)
    
    dic['n'] = N
    
    dic['Lx'] = Lx
    dic['Ly'] = Ly
    dic['Lz'] = Lz

    
    return dic
    
def start_md(**kwargs):
    
    read_from_file = kwargs.get("read_from_file", True)
    n = kwargs.get("n", 100)
    L = kwargs.get("L", np.array([5.0, 5.0, 5.0]))
    
    path_input = kwargs.get("path_input", "pyMD.inp")
    
    conf0_file = kwargs.get("conf0_file", "mdorg.dat")
    output_conf = kwargs.get("output_conf", "md-w.dat")
    output_gr = kwargs.get("output_gr", "grmd.dat")
    num_blocks = kwargs.get("num_blocks", 10)
    steps_per_block = kwargs.get("steps_per_block", 100)
    time_step = kwargs.get("time_step", 0.005)
    temp_option = kwargs.get("temp_option", True)
    temperature = kwargs.get("temperature", 0.7)
    potential_cutoff = kwargs.get("potential_cutoff", 3.0)
    
   
    if read_from_file:
        
        f = open(path_input, 'r')
        for line in f:
            if "INITIAL CONF FILE" in line:
                conf0_file = line.split()[-1]
        f.close()
        
        f = open(path_input, 'r')
        for line in f:
            if "OUTPUT CONF FILE" in line:
                output_conf = line.split()[-1]
        f.close()
        
        f = open(path_input, 'r')
        for line in f:
            if "OUTPUT G(R) FILE" in line:
                output_gr = line.split()[-1]
        f.close()    
        
        f = open(path_input, 'r')
        for line in f:
            if "NUM BLOCKS" in line:
                num_blocks = int(line.split()[-1])
        f.close()  
        
        f = open(path_input, 'r')
        for line in f:
            if "STEPS PER BLOCK" in line:
                steps_per_block = int(line.split()[-1])
        f.close()              
        
        f = open(path_input, 'r')
        for line in f:
            if "TIME STEP" in line:
                time_step = float(line.split()[-1])
        f.close() 
        
        f = open(path_input, 'r')
        for line in f:
            if "TEMPERATURE OPTION" in line:
                temp_option = bool(line.split()[-1])
        f.close() 
        
        f = open(path_input, 'r')
        for line in f:
            if "REQUIRED TEMPERATURE" in line:
                temperature = float(line.split()[-1])
        f.close()    
        
        f = open(path_input, 'r')
        for line in f:
            if "POTENTIAL CUTOFF" in line:
                potential_cutoff = float(line.split()[-1])
        f.close()           
        
    else:
        pass
    
    dic = {}
    
    dic['conf0_file'] = conf0_file
    dic['output_conf'] = output_conf
    dic['output_gr'] = output_gr
    dic['num_blocks'] = num_blocks
    dic['steps_per_block'] = steps_per_block
    dic['time_step'] = time_step
    dic['temp_option'] = temp_option
    dic['temperature'] = temperature
    dic['potential_cutoff'] = potential_cutoff
    
    dens = np.divide(n, np.prod(L))
    ro = L[0]/200.0
    nrad = int(proper_round(L[0]/(2.0*ro)))
    
    dic['dens'] = dens
    dic['ro'] = ro
    dic['nrad'] = nrad
        
    path_mdf_inp = "md.inp"
    path_mdf_inp = kwargs.get("path_mdf_inp", path_mdf_inp)


    keys_start = ['conf0_file', 'output_conf', 'output_gr', 'num_blocks',\
            'steps_per_block', 'time_step', 'temp_option', 'temperature', 'potential_cutoff']

    f = open(path_mdf_inp, 'w')
    for key in keys_start:
        value = dic[key]
        f.write("{}\n".format(value))
    f.close()
    
    return dic
    
    
    
    

    
    
    
