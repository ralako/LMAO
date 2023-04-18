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
r_num = 40
r_moznosti = np.logspace(math.log10(r_start), math.log10(r_stop), num=r_num)


# Hmotnost jednej castice (kg)
mC_start = 1e-18
mC_stop = 1e-9
mC_num = 40
mC_moznosti = np.logspace(math.log10(mC_start), math.log10(mC_stop), num=mC_num)


# Teplota (K)
T = 200


# Rozmery matice simulacii
rozmery = str(r_num)+"x"+str(mC_num)

#------------



#Zaokruhlovacia funkcia
def round_sig(x, sig=3):
    return round(x, sig-int(math.floor(math.log10(abs(x))))-1)



#THREADING

def simulacia(r_index, mC_index):
    oblak = ROS.Oblak(round_sig(r_moznosti[r_index]),round_sig(mC_moznosti[mC_index]),
                      T,uloz=True,video=False,fituj=True,model="1D",rozmery=rozmery)
    oblak.spracovanie()
    return (r_index, mC_index, oblak.rhos_stred_vysl, oblak.sigmy_vysl)
    

def main():
    
    vysledky_rhos = [[0 for i in range(mC_num)] for j in range(r_num)]
    vysledky_sigmy = [[0 for i in range(mC_num)] for j in range(r_num)]
    vysledky_rhos_po10s = [[0 for i in range(mC_num)] for j in range(r_num)]
    vysledky_sigmy_po10s = [[0 for i in range(mC_num)] for j in range(r_num)]

    r_index_pole = []
    mC_index_pole = []
    for i in range(r_num):
        for j in range(mC_num):
            r_index_pole.append(i)
            mC_index_pole.append(j)
            
    with ProcessPoolExecutor() as executor:
        for result in executor.map(simulacia, r_index_pole, mC_index_pole):
            
            r_index, mC_index, rhos_stred_vysl, sigmy_vysl = result

            #print(rhos_stred_vysl)
            #print(rhos_stred_vysl[11])
            
            vysledky_rhos[r_index][mC_index] = rhos_stred_vysl
            vysledky_sigmy[r_index][mC_index] = sigmy_vysl
            vysledky_rhos_po10s[r_index][mC_index] = rhos_stred_vysl[11]
            vysledky_sigmy_po10s[r_index][mC_index] = sigmy_vysl[11]
            
            #print(vysledky_rhos_po10s)
            
            print("Done: "+str(round_sig(r_moznosti[r_index]))+" cm, "+
                  str(round_sig(mC_moznosti[mC_index]))+" kg")


    # ulozenie
    #print("Zacinam ukladat")
    #print(vysledky_rhos_po10s)
    
    vysledky_nazvy = ["vysledky_rhos_"+rozmery+"_"+str(T)+"K.npy",
                      "vysledky_sigmy_"+rozmery+"_"+str(T)+"K.npy",
                      "vysledky_rhos_po10s_"+rozmery+"_"+str(T)+"K.npy",
                      "vysledky_sigmy_po10s_"+rozmery+"_"+str(T)+"K.npy",
                      "config_"+rozmery+"_"+str(T)+"K.npy"]
    
    vysledky = [vysledky_rhos,vysledky_sigmy,vysledky_rhos_po10s,vysledky_sigmy_po10s,
                np.array(r_moznosti,mC_moznosti)]
    
    if not os.path.exists("vysledky_"+rozmery+"_"+str(T)+"K"):
        os.makedirs("vysledky_"+rozmery+"_"+str(T)+"K")
                
    for i in range(len(vysledky_nazvy)):
        if os.path.exists("vysledky_"+rozmery+"_"+str(T)+"K/"+vysledky_nazvy[i]):
            os.remove("vysledky_"+rozmery+"_"+str(T)+"K/"+vysledky_nazvy[i])
        np.save("vysledky_"+rozmery+"_"+str(T)+"K/"+vysledky_nazvy[i], np.array(vysledky[i]))
            





if __name__ == '__main__':
    main()



    
stop = time.perf_counter()
print("------\nCelkový čas: "+str(int(stop-start))+" s")

