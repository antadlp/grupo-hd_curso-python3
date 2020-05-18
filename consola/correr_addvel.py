import subprocess

PATH_ADDVEL_BIN = "./addvel.exe"

coor_file = "azar.dat" # no necesariamente este nombre
addvel_file = "salida.dat"
NUM_PARTICULAS = 10
NON_FROZEN = 10
TEMPERATURA = 1.0

INPUT_VEL_NAME = "addvelX.inp"

f = open(INPUT_VEL_NAME, 'w')
f.write("{}\n{}\n{}\n{}\n{}".format(coor_file,\
                                    addvel_file,\
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


