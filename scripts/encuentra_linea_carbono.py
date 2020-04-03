path = '../inputs/isoctano.xyz'
file = open(path, 'r')

# 'w' escribir o crear archivo
# 'a' anexa sobre el ultimo renglon del archivo
# 'r' de read leer archivo

counter = 1

xc5 = 0
yc5 = 0
zc5 = 0
for line in file:
    
    if 'C' in line:
        
        if counter == 5:
            
            # en linea de interes
            # realizar instrucciones
            
            xc5 = float(line.split()[1])
            yc5 = float(line.split()[2])
            zc5 = float(line.split()[3])
            
            
            
        
        counter+=1
        
    


    
    
print("x: {}, y: {}, z: {}".format(xc5, yc5, zc5))   
    
    

