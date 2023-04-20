import ROS  # Rednutie Oblaku Simulacia
import matplotlib.pyplot as plt

r = 5       #cm
mC = 1e-17  #kg
T = 200     #K

oblak = ROS.Oblak(r,mC,T,uloz=True,video=True,fituj=True)
oblak.spracovanie()

#oblak = ROS.Oblak_BD(r,mC,T)
#x_H,t_H = oblak.animacia()
