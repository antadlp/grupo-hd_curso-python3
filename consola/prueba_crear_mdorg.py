import sys
sys.path.append("/home/antadlp/Documents/grupo-hd_curso-python3/consola")
from conf_inicial_md import *
import subprocess

PATH_COORDS = "coords.dat"
PATH_ADDVEL_BIN = "./addvel.exe"

ADDVEL_FILE = "mdorg.dat"
INPUT_VEL_NAME = "addvelX.inp"

NUM_PARTICULAS = 33
NON_FROZEN = NUM_PARTICULAS
TEMPERATURA = 0.8
DENS = 0.3 

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
                                    ADDVEL_FILE,\
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








