import pandas as pd
import numpy as np

def vv_posicion(x, v, F, m, dt):
    x_next = x + v*dt + np.divide(F, 2*m)*np.power(dt, 2)
    return x_next

def vv_velocidad(v, F_anterior, F_actual, m, dt):
    v_next = v + np.divide(F_anterior + F_actual, 2*m)*dt
    return v_next
    
def vv_posicion2d(**kwargs):
    
    df = kwargs.get('df', 'NULL')
    dt = kwargs.get('dt', 0.01)
    Lx = kwargs.get('Lx', 20.0)
    Ly = kwargs.get('Ly', 20.0)
    
    m = df['m'].values
    r = df[['x', 'y']].values
    v = df[['vx', 'vy']].values
    F = df[['fx', 'fy']].values
    
    r_next = r + v*dt + np.power(dt, 2)*np.divide(F.T, 2*m).T
    
    dic = {}
    dic['x'] = {}
    dic['y'] = {}
    for idx in df.index:
        dic['x'][idx] = np.mod(r_next[idx][0], Lx)
        dic['y'][idx] = np.mod(r_next[idx][1], Ly)
        
    _ = pd.DataFrame(dic)
    
    df.update(_)
    
    return 

def vv_velocidad2d(**kwargs):
    
    df = kwargs.get('df', 'NULL')
    dt = kwargs.get('dt', 0.01)
    df_fb = kwargs.get('df_fb', 'NULL')
    
    m = df['m'].values
    r = df[['x', 'y']].values
    v = df[['vx', 'vy']].values
    fa = df[['fx', 'fy']].values
    fb = df_fb[['fx', 'fy']].values
    
    v_next = v + dt*np.divide((fa+fb).T, 2*m).T
    
    dic = {}
    dic['vx'] = {}
    dic['vy'] = {}
    dic['fx'] = {}
    dic['fy'] = {}
    for idx in df.index:
        _ = v_next[idx][0]
        dic['vx'][idx] = np.sign(_)*np.mod(_, 5.0)
        _ = v_next[idx][1]
        dic['vy'][idx] = np.sign(_)*np.mod(_, 5.0)
        dic['fx'][idx] = df_fb['fx'].loc[idx]
        dic['fy'][idx] = df_fb['fy'].loc[idx]
        
    _ = pd.DataFrame(dic)
    df.update(_)    
    
    return    
