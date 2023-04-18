import ROS  # Rednutie Oblaku Simulacia
import matplotlib.pyplot as plt

r = 5       #cm
mC = [1e-17,1e-16,1e-15,1e-14,1e-13]  #kg
T = 200     #K

for mC_ in mC:
    oblak = ROS.Oblak(r,mC_,T,uloz=True,video=False,fituj=False)
    oblak.spracovanie()
    plt.plot(oblak.casy,oblak.x_H_od_t,label="$m_C$ = "+str(mC_))


plt.title("$r$ = "+str(r)+" cm, $T$ = "+str(T)+"K")
plt.xlabel("$t$ (s)")
plt.ylabel("$x_H$ (cm)")

plt.legend()

plt.savefig("x_H_od_t.png",dpi=300)
