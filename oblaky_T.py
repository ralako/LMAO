import ROS  # Rednutie Oblaku Simulacia
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time
import math
import os


start = time.perf_counter()


#---CONFIG---

# Polomer (cm)
r = 1


# Hmotnost jednej castice (kg)
mC_start = 1e-18
mC_stop = 1e-9
mC_num = 50
mC_moznosti = np.logspace(math.log10(mC_start), math.log10(mC_stop), num=mC_num)


# Teplota (K)
T_start = 50
T_stop = 350
T_num = 50
T_moznosti = np.linspace(T_start, T_stop, num=T_num)


# Rozmery matice simulacii
rozmery = str(T_num)+"x"+str(mC_num)

#------------



#Zaokruhlovacia funkcia
def round_sig(x, sig=3):
    return round(x, sig-int(math.floor(math.log10(abs(x))))-1)



#THREADING

def simulacia(T_index, mC_index):
    oblak = ROS.Oblak(r, round_sig(mC_moznosti[mC_index]), round_sig(T_moznosti[T_index]),
                                uloz=False, video=False, fituj=False, model="1D", rozmery=rozmery)
    oblak.spracovanie()
    return (T_index, mC_index, oblak.rhos_stred_vysl, oblak.sigmy_vysl, oblak.zac_fitu_A_sigma_cas)
    

def main():
    
    vysledky_rhos = [[0 for i in range(mC_num)] for j in range(T_num)]
    vysledky_sigmy = [[0 for i in range(mC_num)] for j in range(T_num)]
    vysledky_rhos_po10s = [[0 for i in range(mC_num)] for j in range(T_num)]
    vysledky_sigmy_po10s = [[0 for i in range(mC_num)] for j in range(T_num)]
    vysledky_zac_fitu_A_sigma_cas = [[0 for i in range(mC_num)] for j in range(T_num)]

    T_index_pole = []
    mC_index_pole = []
    for i in range(T_num):
        for j in range(mC_num):
            T_index_pole.append(i)
            mC_index_pole.append(j)
            
    with ProcessPoolExecutor() as executor:
        for result in executor.map(simulacia, T_index_pole, mC_index_pole):
            
            T_index, mC_index, rhos_stred_vysl, sigmy_vysl, zac_fitu_A_sigma_cas = result

            #print(rhos_stred_vysl)
            #print(rhos_stred_vysl[11])
            
            vysledky_rhos[T_index][mC_index] = rhos_stred_vysl
            vysledky_sigmy[T_index][mC_index] = sigmy_vysl
            vysledky_rhos_po10s[T_index][mC_index] = rhos_stred_vysl[11]
            vysledky_sigmy_po10s[T_index][mC_index] = sigmy_vysl[11]
            vysledky_zac_fitu_A_sigma_cas[T_index][mC_index] = zac_fitu_A_sigma_cas
            
            #print(vysledky_rhos_po10s)
            
            print("Done: "+str(round_sig(T_moznosti[T_index]))+" K, "+
                  str(round_sig(mC_moznosti[mC_index]))+" kg")


    # ulozenie
    #print("Zacinam ukladat")
    #print(vysledky_rhos_po10s)
    
    vysledky_nazvy = ["vysledky_rhos_"+rozmery+"_"+str(r)+"cm.npy",
                      "vysledky_sigmy_"+rozmery+"_"+str(r)+"cm.npy",
                      "vysledky_rhos_po10s_"+rozmery+"_"+str(r)+"cm.npy",
                      "vysledky_sigmy_po10s_"+rozmery+"_"+str(r)+"cm.npy",
                      "vysledky_zac_fitu_A_sigma_cas_"+rozmery+"_"+str(r)+"cm.npy",
                      "config_"+rozmery+"_"+str(r)+"cm.npy"]
    
    vysledky = [vysledky_rhos,vysledky_sigmy,vysledky_rhos_po10s,vysledky_sigmy_po10s,
                vysledky_zac_fitu_A_sigma_cas,np.array([T_moznosti,mC_moznosti],dtype=object) ]
    
    if not os.path.exists("vysledky_"+rozmery+"_"+str(r)+"cm"):
        os.makedirs("vysledky_"+rozmery+"_"+str(r)+"cm")
                
    for i in range(len(vysledky_nazvy)):
        if os.path.exists("vysledky_"+rozmery+"_"+str(r)+"cm/"+vysledky_nazvy[i]):
            os.remove("vysledky_"+rozmery+"_"+str(r)+"cm/"+vysledky_nazvy[i])
        np.save("vysledky_"+rozmery+"_"+str(r)+"cm/"+vysledky_nazvy[i], np.array(vysledky[i]))
            





if __name__ == '__main__':
    main()



    
stop = time.perf_counter()
print("------\nCelkový čas: "+str(int(stop-start))+" s")

