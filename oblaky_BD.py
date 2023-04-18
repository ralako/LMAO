import ROS  # Rednutie Oblaku Simulacia
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time
import math
import os


start = time.perf_counter()


#---CONFIG---

# Polomer (cm)
r_start = 0.1
r_stop = 50
r_num = 50
r_moznosti = np.logspace(math.log10(r_start), math.log10(r_stop), num=r_num)


# Hmotnost jednej castice (kg)
mC_start = 1e-18
mC_stop = 1e-9
mC_num = 50
mC_moznosti = np.logspace(math.log10(mC_start), math.log10(mC_stop), num=mC_num)


# Teplota (K)
T = 200

# Cas simulacii
t_fin = 30

# Rozmery matice simulacii
rozmery = str(r_num)+"x"+str(mC_num)

#------------



#THREADING

def simulacia(r_index, mC_index):
    oblak = ROS.Oblak_BD(ROS.round_sig(r_moznosti[r_index]),ROS.round_sig(mC_moznosti[mC_index]),
                         T,t_fin=t_fin,rozmery=rozmery)
    return (r_index, mC_index, t_fin, oblak.animacia())
    

def main():
    
    x_Hs = [[0 for i in range(r_num)] for j in range(mC_num)]
    t_Hs = [[0 for i in range(r_num)] for j in range(mC_num)]
    Ks = [[0 for i in range(r_num)] for j in range(mC_num)]
    
    r_index_pole = []
    mC_index_pole = []
    for i in range(mC_num):
        for j in range(r_num):
            r_index_pole.append(j)
            mC_index_pole.append(i)
            
    with ProcessPoolExecutor() as executor:
        for result in executor.map(simulacia, r_index_pole, mC_index_pole):
            
            r_index, mC_index, t_fin, (x_Hs[mC_index][r_index], t_Hs[mC_index][r_index], Ks[mC_index][r_index]) = result
            
            print("Done: "+str(ROS.round_sig(r_moznosti[r_index]))+" cm, "+
                  str(ROS.round_sig(mC_moznosti[mC_index]))+" kg")


    ROS.uloz("data_BD_"+rozmery+"_"+str(T)+"K/x_Hs.npy",x_Hs)
    ROS.uloz("data_BD_"+rozmery+"_"+str(T)+"K/t_Hs.npy",t_Hs)
    ROS.uloz("data_BD_"+rozmery+"_"+str(T)+"K/Ks.npy",Ks)
    ROS.uloz("data_BD_"+rozmery+"_"+str(T)+"K/config.npy",np.array([r_moznosti,mC_moznosti,t_fin],dtype=object))
            


if __name__ == '__main__':
    main()



    
stop = time.perf_counter()
print("------\nCelkový čas: "+str(int(stop-start))+" s")

