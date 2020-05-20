import sys
sys.path.append("/home/antadlp/Documents/grupo-hd_curso-python3/consola")
from conf_inicial_md import *

PATH_COORDS = "coords.dat"

N = 20 
dens = 0.6

dic = conf_inicial_md(N=N, dens=dens, iter_check=15)
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








