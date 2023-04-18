import ROS  # Rednutie Oblaku Simulacia
import matplotlib.pyplot as plt

r = 26.5       #cm
mC = 4.5e-17  #kg
T = 200     #K

oblak = ROS.Oblak(r,mC,T,uloz=False,video=True,fituj=False)
oblak.spracovanie()

#oblak = ROS.Oblak_BD(r,mC,T)
#x_H,t_H = oblak.animacia()
