import numpy as np
import matplotlib.pyplot as plt
import ROS

doby_zivota = [
  [  2.41,	1.16,	0.43,	0.12,	0,	0,	0,	0,	0,	0,	0   ],
  [  5.42,	2.71,	1.08,	0.34,	0.08,	0,	0,	0,	0,	0,	0   ],
  [  11.95,	6.10,	2.80,	0.90,	0.25,	0.05,	0,	0,	0,	0,	0   ],
  [  25.03,	13.34,	5.86,	2.16,	0.67,	0.17,	0,	0,	0,	0,	0   ],
  [  50.41,	27.76,	12.71,	4.94,	1.64,	0.46,	0.10,	0,	0,	0,	0   ],
  [  95.06,	55.09,	26.45,	10.77,	3.81,	1.18,	0.30,	0.06,	0,	0,	0   ],
  [  168.52,	103.51,	52.50,	22.65,	8.45,	2.76,	0.79,	0.18,	0,	0,	0   ],
  [  200,	179.33,	98.31,	45.19,	17.80,	6.15,	1.89,	0.50,	0.10,	0,	0   ],
  [  200,	200,	173.23,	86.12,	36.18,	13.26,	4.30,	1.24,	0.30,	0.05,	0   ],
  [  200,	200,	200,	153.40,	69.75,	27.11,	9.34,	2.87,	0.77,	0.17,	0   ],
  [  200,	200,	200,	200,	126.75,	53.15,	19.41,	6.29,	1.84,	0.46,	0.08]
]

Dv_moznosti = np.linspace(0, 200, num=11)
h_moznosti = np.linspace(500, 1000, num=11)  # km

X,Y = np.meshgrid(Dv_moznosti,h_moznosti)
Z = np.array(doby_zivota)

fig, ax = plt.subplots(dpi=300)
levels = np.linspace(0, 200, 9)

cs = ax.contourf(X,Y,Z, levels=levels,cmap='YlGn_r')
cbar = fig.colorbar(cs,ticks=np.linspace(0,200,5))
cbar.ax.set_yticklabels(['0','50','100','150','>200'])

plt.title("doba života $t_Z$ (roky)")

mD = 1.36 ; C = 0.5 ; S = 0.0000785
plt.text(0.55,0.75,"$C_D=$"+str(C)+",\n$m_D=$"+str(mD)+" g,\n$S_D=$"+str(ROS.round_sig(S*10000))+" cm$^2$",
         transform=fig.transFigure,c='white')

plt.xlabel("$\Delta v$ (m/s)")
plt.ylabel("$h$ (km)")

#plt.savefig("tZs.png", dpi=300)
#plt.show()


hranica_25r_Dv = [0,  11, 20, 36, 48, 63, 74, 91, 102,116]
hranica_25r_h  = [649,672,691,733,768,813,847,905,949,1000]

fig, ax = plt.subplots(dpi=300)
plt.grid(ls='--')
plt.plot(hranica_25r_h,hranica_25r_Dv, c='darkgreen')

plt.ylabel("$\Delta v$ (m/s)")
plt.xlabel("$h$ (km)")
plt.title("hranica $t_Z=25$ rokov")

plt.ylim(0,120)
plt.xlim(650,1000)

plt.savefig("hranica_25r.png", dpi=300)
