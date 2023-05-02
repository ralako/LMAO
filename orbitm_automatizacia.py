import pyautogui as pag
from time import sleep
import numpy as np
import ROS

# PUSTAT NA EXPANDNUTOM OKNE


# konstanty
GM = 3.986e14   # zakladne SI
R = 6378e3  # m


# h (m)
h_start = 500e3
h_stop = 1000e3
h_num = 11
h_moznosti = np.linspace(h_start, h_stop, num=h_num)

# Dv (m/s)
Dv_start = 0
Dv_stop = 200
Dv_num = 11
Dv_moznosti = np.linspace(Dv_start, Dv_stop, num=Dv_num)


# Safety stuff
pag.FAILSAFE = True
sleep(5)


def a_calc(h,Dv):
    vD = np.sqrt(GM/(R+h))
    menovatel = vD*vD + 2*vD*Dv - Dv*Dv
    return GM / menovatel


def e_calc(Q,a):
    return np.abs(Q/a - 1)


def klikacka(a,e,t):
    sleep(1)
    # velka poloos (km)
    pag.click(x=475, y=511)
    pag.typewrite(10*['backspace'], interval=0.02)
    pag.typewrite(str(ROS.round_sig(a/1000,sig=5)))

    # excentricita
    pag.click(x=478, y=548)
    pag.typewrite(15*['backspace'], interval=0.02)
    pag.typewrite(str(ROS.round_sig(e,sig=5)))

    # dlzka simulacie
    pag.click(x=445, y=408)
    pag.typewrite(['backspace','delete','delete'], interval=0.05)
    if t < 10:
        pag.typewrite("00"+str(t))
    elif t < 100:
        pag.typewrite("0"+str(t))
    else:
        pag.typewrite(str(t))

    # Tolerance Band (km)
    pag.click(x=444, y=713)
    pag.typewrite(10*['backspace'], interval=0.02)
    pag.typewrite(str(ROS.round_sig(h/1000,sig=5)))

    # Run
    pag.click(x=1386, y=240)
    sleep(0.5)
    pag.typewrite(['enter','enter'], interval=0.5)
    sleep(t)
    pag.click(x=670, y=662)
    sleep(0.5)
    pag.typewrite(['enter','enter'], interval=0.5)
    sleep(0.5)

    # Save
    pag.click(x=1064, y=835)
    pag.typewrite("h-"+str(ROS.round_sig(h/1000))+"_Dv-"+str(Dv)+"_t-"+str(t))
    sleep(1)
    pag.press('enter') ; sleep(1)
    pag.press('enter') ; sleep(1)

    # Clear plots
    pag.click(x=1206, y=236) ; sleep(1)
    pag.press('enter') ; sleep(1)


#--------------------------


for h in h_moznosti:
   for Dv in Dv_moznosti:

    h = h_moznosti[i]
    Dv = Dv_moznosti[i]
    
    a = a_calc(h,Dv)
    e = e_calc(R+h,a)
    q = a*(1-e)
    t = int((q-R)/1000-420)
    if t > 200:
        t = 200
    elif t < 2:
        t = 2
    klikacka(a, e, t)
    print("done",h,Dv,t)
    #print(h,Dv,a,e)
    #input()
