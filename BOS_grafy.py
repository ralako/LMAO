import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import ROS
import BOS

grid = False
uloz = False

# Graf prelet stredom farebny
rozmery_Dy0 = "50x50_prelet"
m_moznosti, mC_moznosti, mD,C,S,v = np.load("config_"+rozmery_Dy0+"_Dy0.npy", allow_pickle=True)
X_Dy0,Y_Dy0 = np.meshgrid(m_moznosti,mC_moznosti)
Dvs_Dy0 = np.load("Dvs_"+rozmery_Dy0+"_Dy0.npy")

fig, ax = plt.subplots(dpi=300)
levels = np.logspace(0, 3, 10)
locator = ticker.LogLocator(base=10)
Z_Dy0 = np.array(Dvs_Dy0)
cs = ax.contourf(X_Dy0,Y_Dy0,Z_Dy0,locator=locator,levels=levels,cmap='inferno')
cbar = fig.colorbar(cs,ticks=np.logspace(0,3,4))

plt.title("$\Delta v$ (m/s)")
plt.text(0.15,0.77,"$\Delta y=0$, $T=200$ K, $t_P=10$ s, $v=$"+str(v)+" km/s\n$C_D=$"+str(C)+", $m_D=$"+str(mD)+" g, $S_D=$"+str(ROS.round_sig(S*10000))+" cm$^2$",
         transform=fig.transFigure)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$m$ (g)")
plt.ylabel("$m_C$ (kg)")
if grid == True:
    plt.grid(linestyle = '--', linewidth = 0.5)

if uloz == True:
    if grid == True:
        plt.savefig("Dvs_"+rozmery_Dy0+"_Dy0_"+str(mD)+"_grid.png", dpi=300)
    else:
        plt.savefig("Dvs_"+rozmery_Dy0+"_Dy0_"+str(mD)+".png", dpi=300)
else:
    plt.show()



# Graf sigmy
fig, ax = plt.subplots(dpi=300)
sigmy = np.load("sigmy_"+rozmery_Dy0+"_Dy0.npy")
plt.plot(mC_moznosti,sigmy,lw=1,c='black',label="$\sigma$ (cm)")
plt.xscale("log")
plt.title("$\sigma(m_C)$")
#plt.show()
plt.savefig("sigmy_"+rozmery_Dy0+"_Dy0.png", dpi=300)




# Graf prelet bokom
rozmery_mCfix = "50x50_prelet_bokom"
m_moznosti, Dy_moznosti, mC = np.load("config_"+rozmery_mCfix+"_mCfix.npy",allow_pickle=True)
X_mCfix,Y_mCfix = np.meshgrid(m_moznosti,Dy_moznosti)
Dvs_mCfix = np.load("Dvs_"+rozmery_mCfix+"_mCfix.npy")

fig, ax = plt.subplots(dpi=300)
levels = np.logspace(0, 3, 10)
locator = ticker.LogLocator(base=10)
Z_mCfix = np.array(Dvs_mCfix)
cs = ax.contourf(X_mCfix,Y_mCfix,Z_mCfix,locator=locator,levels=levels,cmap='inferno')
cbar = fig.colorbar(cs,ticks=np.logspace(0,3,4))

plt.title("$\Delta v$ (m/s)")
plt.text(0.15,0.78,"$m_C=$"+str(mC)+" kg\n$r=$ "+str(ROS.round_sig(BOS.prech_obl(mC)))+" cm",
         transform=fig.transFigure)
plt.text(0.15,0.57,"$T=200$ K\n$t_P=10$ s\n$v=$"+str(v)+" km/s\n$C_D=$"+str(C)+",\n$m_D=$"+str(mD)+" g,\n$S_D=$"+str(ROS.round_sig(S*10000))+" cm$^2$",
         size='small',transform=fig.transFigure)
plt.xscale("log")
plt.xlabel("$m$ (g)")
plt.ylabel("$\Delta y$ (cm)")
if grid == True:
    plt.grid(linestyle = '--', linewidth = 0.5)
#plt.ylim(0,50)
#plt.show()

if uloz == True:
    if grid == True:
        plt.savefig("Dvs_"+rozmery_mCfix+"_mCfix_"+str(mD)+"_"+str(mC)+"_grid.png", dpi=300)
    else:
        plt.savefig("Dvs_"+rozmery_mCfix+"_mCfix_"+str(mD)+"_"+str(mC)+".png", dpi=300)
else:
    plt.show()
