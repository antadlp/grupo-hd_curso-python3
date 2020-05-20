import sys
sys.path.append("/home/antadlp/Documents/grupo-hd_curso-python3/consola")
from conf_inicial_2d import *

N = 5
dens = 1.0
iter_check = 40
elements = ['He']

dic = conf_inicial_2d(N=N, \
        dens=dens, elements=elements, \
        iter_check=iter_check)

df = dic['df']

df.to_excel("conf.xlsx")
df.to_csv("conf.csv")

