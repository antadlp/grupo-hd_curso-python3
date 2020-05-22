import sys
sys.path.append("/home/antadlp/Documents/grupo-hd_curso-python3/consola4")
from metodos_md import *
import subprocess

n = 9
mdorg_file = "mdorg.dat" # CAMBIAR AQUI

crear_mdorg(num_particulas=n,
            mdorg_file = mdorg_file)
            
dic_conf = config_md(mdorg_file=mdorg_file)
            
n = dic_conf['n']
Lx = dic_conf['Lx']
Ly = dic_conf['Ly']
Lz = dic_conf['Lz']
L = [Lx, Ly, Lz]

dic_start = start_md(n=n,
                     L=L)

rcut = dic_start['potential_cutoff']

rcutsq = rcut*rcut
r2ij = np.divide(1.0, rcutsq)
r6ij = r2ij*r2ij*r2ij
r12ij = r6ij*r6ij
vcut = r12ij - r6ij

pot = 0.0
vir = 0.0

fx = np.zeros(n)
fy = np.zeros(n)
fz = np.zeros(n)

rx = dic_conf['rx']
ry = dic_conf['ry']
rz = dic_conf['rz']

for i in range(n-1):
    rxi = rx[i]
    ryi = ry[i]
    rzi = rz[i]
    
    fxi = fx[i]
    fyi = fy[i]
    fzi = fz[i]
    
    jlist = np.arange(i+1, n, 1)
    for j in jlist:
        rxij = rxi - rx[j]
        ryij = ryi - ry[j]
        rzij = rzi - rz[j]
        
        rxij = rxij - proper_round(rxij/Lx)*Lx
        ryij = ryij - proper_round(ryij/Ly)*Ly
        rzij = rzij - proper_round(rzij/Lz)*Lz
        
        rijsq = rxij*rxij + ryij*ryij + rzij*rzij
        
        if (rijsq < rcutsq):
            r2ij = np.divide(1.0, rijsq)
            r6ij = r2ij*r2ij*r2ij
            r12ij = r6ij*r6ij
            vij = r12ij - r6ij
            fij = (vij + r12ij)*r2ij
            vij = vij - vcut
            pot = pot + vij
            fxij = fij*rxij
            fyij = fij*ryij
            fzij = fij*rzij
            fxi = fxi + fxij
            fyi = fyi + fyij
            fzi = fzi + fzij
            vir = vir + rxij*fxij + ryij*fyij + rzij*fzij
            fxi   = fxi + fxij
            fyi   = fyi + fyij
            fzi   = fzi + fzij
            fx[j] = fx[j] - fxij
            fy[j] = fy[j] - fyij
            fz[j] = fz[j] - fzij 

            print("{}{}{}".format(fx[j], fy[j], fz[j]))










            


        
        
        
