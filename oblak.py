import ROS  # Rednutie Oblaku Simulacia
import matplotlib.pyplot as plt

#r = 5       #cm
mC = 1.33e-18  #kg
T = 200     #K

def prech_obl(mC):
    return mC**(-0.491) * 10**(-7.8)   # r (cm)

r = prech_obl(mC)

oblak = ROS.Oblak(ROS.round_sig(r),ROS.round_sig(mC),T,
                  uloz=True,video=True,fituj=True)
oblak.spracovanie()

#oblak = ROS.Oblak_BD(r,mC,T)
#x_H,t_H = oblak.animacia()
