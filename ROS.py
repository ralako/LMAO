import matplotlib.pyplot as plt
from matplotlib import ticker
from celluloid import Camera
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import math
import os
import warnings

warnings.filterwarnings("ignore", message="Covariance of the parameters could not be estimated")
warnings.filterwarnings("ignore", message="Log scale: values of z <= 0 have been masked")


class Oblak:
    
    def __init__(self,r,mC,T,uloz=True,video=True,fituj=True,model="1D",opakuj_na_frame=200,n=100,dt=1e-4,t=0,
                 k=1.38065e-19,kraj=1e6,tolerancia_chyby_A=1e-3,t_posun = 1,vysl_poc=31,rozmery="1x1"):
        
        self.r = r              # Pološírka oblaku (cm)
        self.mC = mC            # Hmotnost jednej castice
        self.T = T              # Teplota (K)
        self.uloz = uloz        # Ak True, tak uklada datove subory
        self.video = video      # Ak True, tak generuje video
        self.fituj = fituj      # Ak True, tak fituje relativnu hustotu Gaussianom
        self.model = model      # Typ modelu podla ktoreho pocita (1D, polarny, sfericky)
        
        self.opakuj_na_frame = opakuj_na_frame  # Kolko simulacnych cyklov urobi pocas jedneho frame
        
        self.n = n              # Skalovací faktor elementu resp. pocet elementov
        self.dr = r/n           # Velkost elementu (cm)
        self.dt = dt            # Element casu (s)
        self.t = t              # Aktualny cas v simulacii (s)
        
        self.k = k              # Boltzmannova konstanta (cm)
        self.kTlmC = k*T/mC     # Pomocna konstanta (cm)
        self.kraj = kraj        # Hypoteticky element strasne daleko (cm)

        self.tolerancia_chyby_A = tolerancia_chyby_A    # urcuje zaciatok relevantnej oblasti
        
        self.t_posun = t_posun              # ako casto sa zapisuju vysledky (s)
        self.vysl_poc = vysl_poc            # kolkokrat sa zapisuju vysledky
        self.poc_vysl = 0                   # pocitadlo zapisanych vysledkov
        self.t_fin = t_posun*(vysl_poc-1)   # finalny cas po ktory sa simuluje (s)
        
        self.frame_posun = int(t_posun/(opakuj_na_frame*dt))        # po kolko frames sa zapisuju vysledky
        self.frames_poc = 1 + int(self.t_fin/(opakuj_na_frame*dt))  # celkovy pocet frames
        
        
        self.x = np.linspace(self.dr,r-self.dr,n-1) # Stredy elementov (cm)
        self.v = np.zeros(n-1)                      # Rychlosti stredov elementov (cm/s)
        self.a = np.zeros(n-1)                      # Zrychlenia stredov elementov (cm/s^2)
                                                    # element na pozicii x=0 sa nehybe 
        self.dxL = np.array((n-1)*[self.dr])        # vzidalenost ku stredu elementu nalavo (cm)
        self.dxR = np.array((n-2)*[self.dr]+[kraj]) # vzidalenost ku stredu elementu napravo (cm)
        self.rho = 2*self.dr/(self.dxL+self.dxR)    # relativna hustota
            # v elemente dr*(dxL+dxR) je 1+0.5+0.5 gulicky => faktor 2*
        self.rhos_stred = np.zeros(self.frames_poc) # stredove hustoty
        self.rhos_stred_vysl = np.zeros(vysl_poc)

        self.zac_rel_oblasti_frame = 0

        self.zac_fitu_A_sigma_cas = np.inf


        self.casy = np.linspace(0,self.t_fin,num=self.frames_poc)   # Casy od zaciatku rozpinania
        
        if fituj:
            
            self.Acka = np.zeros(self.frames_poc)                       # Konstanta A z Gaussianu
            self.Acka_se = np.zeros(self.frames_poc)
            self.sigmy = np.zeros(self.frames_poc)                      # Sigma Gaussianu
            self.sigmy_se = np.zeros(self.frames_poc)

            self.Acka_vysl = np.zeros(vysl_poc)
            self.Acka_se_vysl = np.zeros(vysl_poc)
            self.sigmy_vysl = np.zeros(vysl_poc)
            self.sigmy_se_vysl = np.zeros(vysl_poc)
            self.casy_vysl = np.linspace(0,self.t_fin,num=vysl_poc)

            self.fitol_som_hustotu = True

        self.rozmery = rozmery
        self.nazov = "oblak_"+str(self.r)+"cm_"+str(self.mC)+"kg_"+str(self.T)+"K"

        self.x_H_od_t = np.zeros(self.frames_poc)
        self.x_H_od_t[0] = self.r

        if self.mC < 1e-18:
            self.n = 50

        

        #VIDEO
        if self.video:
            self.fig = plt.figure()
            self.y = np.zeros(n-1)
            plt.xlim(0,50) ; plt.ylim(0,115)
            plt.xlabel("$x$ (cm)")
            plt.title("$r$ = "+str(self.r)+" cm, $m_C$ = "+str(self.mC)+" kg, $T$ = "+str(self.T)+" K")
            if uloz: self.fig.set_dpi(300)
            self.camera = Camera(self.fig)

        

    def simuluj_1D(self):
        self.a = self.kTlmC*(1/self.dxL - 1/self.dxR)
        dx = 0.5*self.a*self.dt*self.dt + self.v*self.dt
        self.x += dx ; self.dxL += dx ; self.dxR -= dx
        dx_ = np.append(0,dx)
        dxL_ = np.append(self.dxL,0)
        dxR_ = np.append([0,0],self.dxR)
        for i in range(self.n):
            dxL_[i] -= dx_[i]
            dxR_[i] += dx_[i]
        self.dxL = dxL_[0:-1]
        self.dxR = dxR_[2:]
        self.rho = 2*self.dr/(self.dxL+self.dxR)
        self.v += self.a*self.dt


    def ukonci_predcasne(self,poc):
        self.t_fin = (poc-1)*self.opakuj_na_frame*self.dt
        self.Acka = self.Acka[:poc]
        self.Acka_se = self.Acka_se[:poc]
        self.sigmy = self.sigmy[:poc]
        self.sigmy_se = self.sigmy_se[:poc]
        self.rhos_stred = self.rhos_stred[:poc]
        self.casy = self.casy[:poc]
        self.frames_poc = poc


    def fituj_hustotu(self,poc):     #gaussianom
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            try:
                popt,pcov = curve_fit(gauss,self.x,self.rho)
                popt[1] = abs(popt[1])
            except:
                print(poc, "neviem nafitovat")
                self.fitol_som_hustotu = False
                if self.t > 10:
                    self.ukonci_predcasne(poc)
                    return True
                else:
                    return False

        self.fitol_som_hustotu = True

        # zastavi simulaciu ak je nieco nafitovane ako zaporne a zaroven t > 10s
        if (popt[0] <= 0) or (pcov[0,0] <= 0) or (pcov[1,1] <= 0):
            if self.t > 10:
                self.ukonci_predcasne(poc)
                print(poc, "zly fit", popt[0], popt[1], pcov[0,0], pcov[1,1])
                return True
            
            else:
                self.zac_rel_oblasti_frame = poc+1
            
        else:
            perr = 2*np.sqrt(np.diag(pcov))
            
            self.Acka[poc] = popt[0]
            self.Acka_se[poc] = perr[0]
            self.sigmy[poc] = popt[1]
            self.sigmy_se[poc] = perr[1]

            if (poc % self.frame_posun) == 0:
                self.Acka_vysl[self.poc_vysl] = popt[0]
                self.Acka_se_vysl[self.poc_vysl] = perr[0]
                self.sigmy_vysl[self.poc_vysl] = popt[1]
                self.sigmy_se_vysl[self.poc_vysl] = perr[1]

        return False
    


    def animacia(self):
        for poc in range(self.frames_poc):

            if self.fituj:
                if self.fituj_hustotu(poc):
                    print("stop po "+str(self.t_fin)+" s")
                    break
                
            for j in range(self.opakuj_na_frame):
                if self.model == "1D":
                    self.simuluj_1D()
                self.t += self.dt

            # x_H (hranicne x kde dochadza ku zlomu rho)
            self.rho_rozsirene = np.append([self.dr/self.dxL[0]],self.rho)
            self.x_rozsirene = np.append([0],self.x)
            self.x_H_od_t[poc] = self.x_rozsirene[np.argmax(self.rho_rozsirene < 0.99)]
                
            self.rhos_stred[poc] = self.rho[0]
            #print(self.rhos_stred)
            if (poc % self.frame_posun) == 0:
                self.rhos_stred_vysl[self.poc_vysl] = self.rho[0]
                self.poc_vysl += 1
    
            if self.video:
                plt.plot(self.x,self.y,'|',color='black')
                plt.plot(self.x,100*self.rho,color='black',label="$\\rho_r(x,t)$ (%)")

                if self.fituj:
                    if (self.zac_rel_oblasti_frame < poc) and (self.fitol_som_hustotu == True):
                        plt.plot(self.x,100*gauss(self.x,self.Acka[poc],self.sigmy[poc]),color='red',alpha=0.5,
                                 label="fit $\\rho_r$ Gaussianom")
                
                plt.text(0.7,0.70,"$t$ = "+str('{:.2f}'.format(round(self.casy[poc],2)))+" s",
                         transform=self.fig.transFigure)
                
                self.camera.snap()

        if self.video:
            plt.text(0.6,0.83,"───  \\rho_r(x,t)$ (%)$",transform=self.fig.transFigure)
            plt.text(0.6,0.78,"───",c='red',transform=self.fig.transFigure)
            plt.text(0.654,0.78,"fit $\\rho_r$ Gaussianom",transform=self.fig.transFigure)
            self.ani = self.camera.animate(interval=int(1000*self.opakuj_na_frame*self.dt))




    def ulozenie(self):
        if self.fituj:
            self.zdrojove_cesty = ["data_"+self.rozmery+"/"+self.nazov+"/"+self.nazov+".mp4",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_casy.npy",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_Acka.npy",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_Acka_se.npy",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_sigmy.npy",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_sigmy_se.npy",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_rhos_stred.npy",
                            "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_parametre.npy"]
        else:
            self.zdrojove_cesty = ["data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+".mp4",
                                   "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_casy.npy",
                                   "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_rhos_stred.npy",
                                   "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_parametre.npy"]
        
        # vytvorenie priecinkovej struktury
        if not os.path.exists("data_"+self.rozmery+"_"+str(self.T)+"K"):
            os.makedirs("data_"+self.rozmery+"_"+str(self.T)+"K")
        if not os.path.exists("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov):
            os.makedirs("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov)
            
        #ulozenie videa
        if self.video:
            if os.path.exists(self.zdrojove_cesty[0]):
                os.remove(self.zdrojove_cesty[0])
            self.ani.save(self.zdrojove_cesty[0])
            plt.close(fig=self.fig)
        
        # ulozenie datovych suborov
        self.parametre = np.array([self.r,self.mC,self.T,self.opakuj_na_frame,self.n,self.dr,self.dt,self.t,
                               self.k,self.kTlmC,self.kraj,self.tolerancia_chyby_A,
                               self.t_posun,self.vysl_poc,self.poc_vysl,self.t_fin,
                               self.frame_posun,self.frames_poc])

        if self.fituj:
            data_ = [self.casy,self.Acka,self.Acka_se,self.sigmy,self.sigmy_se,self.parametre,self.rhos_stred]
        else:
            data_ = [self.casy,self.parametre,self.rhos_stred]
            
        for i in range(len(data_)):
            if os.path.exists(self.zdrojove_cesty[i+1]):
                os.remove(self.zdrojove_cesty[i+1])
            np.save(self.zdrojove_cesty[i+1],data_[i])



    def x_H_plot(self):
        
        self.fig3 = plt.figure()
        plt.plot(self.casy,self.x_H_od_t)
        
        plt.title("$r$ = "+str(self.r)+" cm, $m_C$ = "+str(self.mC)+", $T$ = "+str(self.T)+"K")
        plt.xlabel("$t$ (s)")
        plt.ylabel("$x_H$ (cm)")
        
        if self.uloz:
            plt.savefig("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_x_H_od_t.png",dpi=300)
            plt.close(fig=self.fig3)
        else:
            plt.show()
            
    

    def fitovanie(self):
        
        self.fig2 = plt.figure()

        try:
            frame_poc = np.argmax(self.Acka_se[self.zac_rel_oblasti_frame:] < self.tolerancia_chyby_A)
        except:
            print("nie je relevantna oblast Acka")
            frame_poc = 0

        hranica_relevantnosti = 0.95
        if frame_poc == 0 or frame_poc > hranica_relevantnosti * self.frames_poc:    # nebola cely cas resp. dlho (>95%)
            print("fit failed")                                     # dosiahnuta pozadovana presnost na fit
            if self.poc_vysl == self.vysl_poc:
                print("ale je to ok :)")
                
            plt.plot(self.casy,100*self.rhos_stred,c='black',label=r'$\rho_r$ v strede (%)')

            if self.uloz:
                if os.path.exists("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_rhos_stred_vysl.npy"):
                    os.remove("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_rhos_stred_vysl.npy")
                np.save("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_rhos_stred_vysl.npy",
                        self.rhos_stred_vysl)
            
            return
        
        frame_kon = self.frames_poc      # koniec relevantnej oblasti
        
        # hyperbolicky fit A
        poptH,pcovH = curve_fit(hyperbola,self.casy[frame_poc:frame_kon],self.Acka[frame_poc:frame_kon])
        perrH = 2*np.sqrt(np.diag(pcovH))
        zac_kreslenia_Acka = (poptH[1] < 0) * math.ceil(-poptH[1]/(self.opakuj_na_frame*self.dt))
        
        # lineárny fit sigmy
        poptL,pcovL = curve_fit(priamka,self.casy[frame_poc:frame_kon],self.sigmy[frame_poc:frame_kon])
        perrL = 2*np.sqrt(np.diag(pcovL))

        
        # A
        plt.plot(self.casy,100*self.Acka,c='red',label="$A$ (%)")
        plt.plot(self.casy,100*self.Acka_se*10,c='brown',alpha=.5,label="$\sigma_A \cdot 10$ (%)")
        plt.plot(self.casy[zac_kreslenia_Acka:],100*hyperbola(self.casy[zac_kreslenia_Acka:],
                                        poptH[0],poptH[1]),c='orange',label="hyperbolicky fit $A$")
        # sigma
        plt.plot(self.casy,self.sigmy,c='blue',label="$\sigma$ (cm)")
        plt.plot(self.casy,self.sigmy_se*100,c='navy',alpha=.5,label="$\sigma_\sigma$ (cm$\cdot 10^{-2}$)")
        plt.plot(self.casy,poptL[0]*self.casy+poptL[1],c='green',label="lineárny fit $\sigma$")

        # rho relativne
        plt.plot(self.casy,100*self.rhos_stred,c='black',label=r'$\rho_r$ v strede (%)')
        
        # relevantna oblast
        self.zac_fitu_A_sigma_cas = frame_poc*self.opakuj_na_frame*self.dt
        plt.axvline(x=self.zac_fitu_A_sigma_cas,c='black',linestyle='--',lw=1,alpha=.5,
                    label="zaciatok fitu $A$ a $\sigma$")

        #ZAPISANIE VYSLEDKOV
        self.Acka_vysl[self.poc_vysl:] = hyperbola(self.casy_vysl[self.poc_vysl:],poptH[0],poptH[1])
        self.sigmy_vysl[self.poc_vysl:] = priamka(self.casy_vysl[self.poc_vysl:],poptL[0],poptL[1])
        self.rhos_stred_vysl[self.poc_vysl:] = self.Acka_vysl[self.poc_vysl:]
        #print(self.poc_vysl)
        # metoda prenosu chyb
        clen1H = 1/(poptH[1]+self.casy_vysl[self.poc_vysl:])
        self.Acka_se_vysl[self.poc_vysl:] = clen1H * np.sqrt(
            perrH[0]*perrH[0] + (poptH[0]*poptH[0]*perrH[1]*perrH[1])/(clen1H*clen1H))
        self.sigmy_se_vysl[self.poc_vysl:] = np.sqrt(
            perrL[0]*perrL[0]*np.square(self.casy_vysl[self.poc_vysl:]) + perrL[1]*perrL[1])
        #ulozenie
        if self.uloz:
            self.zdrojove_cesty_2 = [
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_casy_vysl.npy",
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_Acka_vysl.npy",
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_Acka_se_vysl.npy",
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_sigmy_vysl.npy",
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_sigmy_se_vysl.npy",
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_rhos_stred_vysl.npy",
                    "data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_konstanty_fitu.npy" ]
            data_ = [self.casy_vysl,self.Acka_vysl,self.Acka_se_vysl,self.sigmy_vysl,self.sigmy_se_vysl,
                     self.rhos_stred_vysl,np.array([[poptH[0],poptH[1],poptL[0],poptL[1]],
                                                    [perrH[0],perrH[1],perrL[0],perrL[1]], self.zac_fitu_A_sigma_cas],dtype=object)]
            
            for i in range(len(data_)):
                if os.path.exists(self.zdrojove_cesty_2[i]):
                    os.remove(self.zdrojove_cesty_2[i])
                np.save(self.zdrojove_cesty_2[i],data_[i])
        else:
            print("B = "+str(poptH[0])+" s ; chyba = "+str(perrH[0])+" s")
            print("C = "+str(poptH[1])+" s ; chyba = "+str(perrH[1])+" s")
            print("D = "+str(poptL[0])+" cm/s ; chyba = "+str(perrL[0])+" cm/s")
            print("E = "+str(poptL[1])+" cm ; chyba = "+str(perrL[1])+" cm")





    def ulozenie_fitu(self):
        
        plt.xlim(0,self.t_fin) ; plt.ylim((0,115))
        plt.title("$r = $"+str(self.r)+" cm, $m_C = $"+str(self.mC)+" kg, $T = $"+str(self.T))
        plt.legend(loc=1)
        plt.xlabel("$t$ (s)")
        
        if self.uloz:
            plt.savefig("data_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"/"+self.nazov+"_graf.png",dpi=300)
            plt.close(fig=self.fig2)
        else:
            plt.show()


        
    def spracovanie(self):
        self.animacia()
        if self.uloz:
            self.ulozenie()
        else:
            plt.show()
        if self.fituj:
            self.fitovanie()
            self.ulozenie_fitu()





#ULOZENIE SUBORU

def uloz(cesta,data):
    if os.path.exists(cesta):
        os.remove(cesta)
    np.save(cesta,data)




#FITOVACIE FUNKCIE
    
def gauss(r, A, sigma):
        return A * np.exp(-0.5*(r*r)/(sigma*sigma))
    
def hyperbola(t, B, C):
    return B/(t+C)

def priamka(t,slope,intercept):
    return slope*t + intercept

def loglogfit(r,P,Q):
    return r**P * 10**Q




#ZAOKRUHLOVACIA FUNKCIA
def round_sig(x, sig=3):
    if x == 0:
        return 0
    else:
        return round(x, sig-int(math.floor(math.log10(abs(x))))-1)




#COLORMESH GRAFY a INÉ GRAFY

def colormesh(rozmery,t,T=200,fps=50,smooth=False,uloz=True):
    
    poc = t*fps
    r_moznosti, mC_moznosti = np.load("vysledky_"+rozmery+"_"+str(T)+"K/config_"+rozmery+"_"+str(T)+"K.npy", allow_pickle=True)
    r_num = r_moznosti.size
    mC_num = mC_moznosti.size
    
    Acka = [[0 for j in range(mC_num)] for i in range(r_num)]
    rhos_stred = [[0 for j in range(mC_num)] for i in range(r_num)]
    zac_fitu_A_sigma_cas = [[0 for j in range(mC_num)] for i in range(r_num)]
    prechodova_oblast = [[0 for j in range(mC_num)] for i in range(r_num)]
    prechodova_oblast_r = [] ; prechodova_oblast_r_log = []
    prechodova_oblast_mC = [] ; prechodova_oblast_mC_log = []
    
    for i in range(mC_num):
        for j in range(r_num):
            nazov = "oblak_"+str(round_sig(r_moznosti[j]))+"cm_"+str(round_sig(mC_moznosti[i]))+"kg_"+str(T)+"K"
            nazov_rhos_stred = "data_"+rozmery+"_"+str(T)+"K/"+nazov+"/"+nazov+"_rhos_stred_vysl.npy"
            nazov_Acka = "data_"+rozmery+"_"+str(T)+"K/"+nazov+"/"+nazov+"_Acka.npy"
            nazov_konstanty_fitu =  "data_"+rozmery+"_"+str(T)+"K/"+nazov+"/"+nazov+"_konstanty_fitu.npy"
            
            Acka[i][j] = np.load(nazov_Acka)[poc]
            rhos_stred[i][j] = np.load(nazov_rhos_stred)[t]
            try:
                zac_fitu_A_sigma_cas[i][j] = np.load(nazov_konstanty_fitu,allow_pickle=True)[2]
                #print(zac_fitu_A_sigma_cas[i][j])
            except:
                zac_fitu_A_sigma_cas[i][j] = np.inf

            if (zac_fitu_A_sigma_cas[i][j] < 10) and (rhos_stred[i][j] > 0.05):
                prechodova_oblast[i][j] = 1
                prechodova_oblast_r.append(r_moznosti[j])
                prechodova_oblast_r_log.append(np.log10(r_moznosti[j]))
                prechodova_oblast_mC.append(mC_moznosti[i])
                prechodova_oblast_mC_log.append(np.log10(mC_moznosti[i]))

                
    
    X,Y = np.meshgrid(np.array(r_moznosti),np.array(mC_moznosti))
    
    Z_Acka = np.array(Acka)
    Z_rhos_stred = np.array(rhos_stred)
    Z_zac_fitu_A_sigma_cas = np.array(zac_fitu_A_sigma_cas)
    Z_prechodova_oblast = np.array(prechodova_oblast)
    
    if smooth:
        levels_Acka = np.linspace(0, 1, 256)
        levels_rhos_stred = np.linspace(0, 1, 256)
        levels_zac_fitu_A_sigma_cas = np.linspace(0, 1, 256)
        levels_prechodova_oblast = [0,.5,1]
    else:
        levels_Acka = np.linspace(0, 1, 11)
        levels_rhos_stred = np.linspace(0, 1, 11)
        levels_zac_fitu_A_sigma_cas = np.linspace(0, 20, 11)
        levels_prechodova_oblast = [0,.5,1]


    # Acka
    fig, ax = plt.subplots(dpi=300)
    cs = ax.contourf(X,Y,Z_Acka,levels=levels_Acka)

    cbar = fig.colorbar(cs,ticks=[0,0.2,0.4,0.6,0.8,1])

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$A(t = $"+str('{:.1f}'.format(round(poc/fps,2)))+" s$)$, $T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/Acka_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/Acka_"+rozmery+"_"+str(T)+"K"+".png")

    else:
        plt.show()


    # Rhos stred
    fig, ax = plt.subplots(dpi=300)
    cs = ax.contourf(X,Y,Z_rhos_stred,levels=levels_rhos_stred)

    #print(rhos_stred)

    cbar = fig.colorbar(cs,ticks=[0,0.2,0.4,0.6,0.8,1])

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$\\rho_r(x = 0$ cm $t = $"+str('{:.1f}'.format(round(poc/fps,2)))+" s$)$, $T$ = "+str(T)+"K")
    
    if uloz:
        if smooth:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/Rhos_stred_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/Rhos_stred_"+rozmery+"_"+str(T)+"K"+".png")

    else:
        plt.show()


    # zac_fitu_A_sigma_cas
    fig,ax = plt.subplots(dpi=300)
    cs = ax.contourf(X,Y,Z_zac_fitu_A_sigma_cas,levels=levels_zac_fitu_A_sigma_cas,cmap='bone')

    cbar = fig.colorbar(cs,ticks=np.linspace(0, 20, 11))

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("zaciatok fitu $A, \sigma$ (s), $T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/zac_fitu_A_sigma_cas_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/zac_fitu_A_sigma_cas_"+rozmery+"_"+str(T)+"K"+".png")

    else:
        plt.show()



    # Prechodova oblast
    popt,pcov = curve_fit(priamka, prechodova_oblast_r_log, prechodova_oblast_mC_log)
    print("Prechodova oblast", popt,pcov)
    
    fig,ax = plt.subplots(dpi=300)
    cs = ax.contourf(X,Y,Z_prechodova_oblast,levels=levels_prechodova_oblast,colors=['white','black'])

    #cbar = fig.colorbar(cs,ticks=[0,1])
    plt.plot(r_moznosti, loglogfit(r_moznosti,popt[0],popt[1]), c='red', label="lineárny fit prechodovej oblasti v log log škále")

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("Prechodova oblast, $T$ = "+str(T)+"K")
    plt.legend()

    if uloz:
        if smooth:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/prechodova_oblast_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("vysledky_"+rozmery+"_"+str(T)+"K"+"/prechodova_oblast_"+rozmery+"_"+str(T)+"K"+".png")

    else:
        plt.show()




def colormesh_anim(rozmery,co,T,poc_moznosti,fps,smooth=False,uloz=True):

    fig,ax = plt.subplots()
    if uloz: fig.set_dpi(300)
    camera = Camera(fig)

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")

    plt.title("$\\rho_r$ v strede, $T$ = "+str(T)+"K")
    
    if smooth:
        levels = np.linspace(0, 1, 256)
    else:
        levels = np.linspace(0, 1, 11)
        
    
    for poc in poc_moznosti:
        r_moznosti, mC_moznosti = np.load("vysledky_"+rozmery+"_"+str(T)+"K/config_"+rozmery+"_"+str(T)+"K.npy")
        r_num = len(r_moznosti)
        mC_num = len(mC_moznosti)
        
        rhos_stred_po_t = [[0 for j in range(mC_num)] for i in range(r_num)]
        
        for i in range(mC_num):
            for j in range(r_num):
                nazov = "oblak_"+str(round_sig(r_moznosti[j]))+"cm_"+str(round_sig(mC_moznosti[i]))+"kg_"+str(T)+"K"
                nazov_rhos_stred_vysl = "data_"+rozmery+"_"+str(T)+"K/"+nazov+"/"+nazov+"_"+co+".npy"
                
                rhos_stred_po_t[i][j] = np.load(nazov_rhos_stred_vysl)[poc]
        
        X,Y = np.meshgrid(np.array(r_moznosti),np.array(mC_moznosti))
        Z = np.array(rhos_stred_po_t)

        #fig.clear()
        cs = plt.contourf(X,Y,Z,levels=levels)
        #cbar = fig.colorbar(cs,ticks=[0,0.2,0.4,0.6,0.8,1])

        plt.text(0.6,0.8,"$t$ = "+str('{:.2f}'.format(round(poc/fps,2)))+" s", transform=fig.transFigure)

        camera.snap()


    ani = camera.animate(interval=int(1000/fps))
    
    if uloz:
        if smooth:
            ani.save("vysledky_"+rozmery+"_"+str(T)+"K"+"/colormesh_anim_"+rozmery+"_"+str(T)+"K"+"_smooth.mp4")
        else:
            ani.save("vysledky_"+rozmery+"_"+str(T)+"K"+"/colormesh_anim_"+rozmery+"_"+str(T)+"K"+".mp4")

    else:
        plt.show()









class Oblak_BD:
    
    def __init__(self,r,mC,T,opakuj_na_frame=200,n=100,dt=1e-4,t=0,
                 k=1.38065e-19,kraj=1e6,t_fin=30,rozmery="1x1"):
        
        self.r = r              # Pološírka oblaku (cm)
        self.mC = mC            # Hmotnost jednej castice
        self.T = T              # Teplota (K)
        
        self.opakuj_na_frame = opakuj_na_frame  # Kolko simulacnych cyklov urobi pocas jedneho frame
        
        self.n = n              # Skalovací faktor elementu resp. pocet elementov
        self.dr = r/n           # Velkost elementu (cm)
        self.dt = dt            # Element casu (s)
        self.t = t              # Aktualny cas v simulacii (s)
        
        self.k = k              # Boltzmannova konstanta (cm)
        self.kTlmC = k*T/mC     # Pomocna konstanta (cm)
        self.kraj = kraj        # Hypoteticky element strasne daleko (cm)
        
        self.t_fin = t_fin      # finalny cas po ktory sa simuluje (s)
        self.frames_poc = 1 + int(self.t_fin/(opakuj_na_frame*dt))  # celkovy pocet frames       
        
        self.x = np.linspace(self.dr,r-self.dr,n-1) # Stredy elementov (cm)
        self.v = np.zeros(n-1)                      # Rychlosti stredov elementov (cm/s)
        self.a = np.zeros(n-1)                      # Zrychlenia stredov elementov (cm/s^2)
                                                    # element na pozicii x=0 sa nehybe 
        self.dxL = np.array((n-1)*[self.dr])        # vzidalenost ku stredu elementu nalavo (cm)
        self.dxR = np.array((n-2)*[self.dr]+[kraj]) # vzidalenost ku stredu elementu napravo (cm)
        self.rho = 2*self.dr/(self.dxL+self.dxR)    # relativna hustota
            # v elemente dr*(dxL+dxR) je 1+0.5+0.5 gulicky => faktor 2*
        self.rhos_stred = np.zeros(self.frames_poc) # stredove hustoty

        self.casy = np.linspace(0,self.t_fin,num=self.frames_poc)   # Casy od zaciatku rozpinania

        self.x_H_od_t = np.zeros(self.frames_poc)   # polohy zlomu v jednotlivych casoch

        self.rozmery = rozmery
        self.nazov = "oblak_"+str(self.r)+"cm_"+str(self.mC)+"kg_"+str(self.T)+"K"

        if self.mC < 1e-18:
            self.n = 50

        

    def simuluj_1D(self):
        self.a = self.kTlmC*(1/self.dxL - 1/self.dxR)
        dx = 0.5*self.a*self.dt*self.dt + self.v*self.dt
        self.x += dx ; self.dxL += dx ; self.dxR -= dx
        dx_ = np.append(0,dx)
        dxL_ = np.append(self.dxL,0)
        dxR_ = np.append([0,0],self.dxR)
        for i in range(self.n):
            dxL_[i] -= dx_[i]
            dxR_[i] += dx_[i]
        self.dxL = dxL_[0:-1]
        self.dxR = dxR_[2:]
        self.rho = 2*self.dr/(self.dxL+self.dxR)
        self.v += self.a*self.dt

    

    def animacia(self):
        
        for poc in range(self.frames_poc):
            for j in range(self.opakuj_na_frame):
                self.simuluj_1D()
                self.t += self.dt

            #rho_stred = self.rho[0]
            rho_stred = self.dr/self.dxL[0]
            self.rhos_stred[poc] = rho_stred

            if (self.dr/self.dxL[0]) < 0.99:
                self.x_H_od_t[poc] = 0
            else:
                self.x_H_od_t[poc] = self.x[np.argmax(self.rho < 0.99)]


        koniec_fitu = self.frames_poc
        for i in range(self.frames_poc):
            if self.x_H_od_t[i] == 0:
                koniec_fitu = i
                break
            
        popt,pcov = curve_fit(lambda t, slope: priamka(t,slope,self.r),
                              self.casy[:koniec_fitu],self.x_H_od_t[:koniec_fitu])
        perr = 2*np.sqrt(np.diag(pcov))
        self.K = popt[0]    # parameter fitu x_H


        self.t_H = np.inf
        for i in range(self.frames_poc):
            if self.rhos_stred[i] < 0.99:
                self.t_H = self.casy[i]
                break
        #self.t_H = self.casy[np.argmax(self.rhos_stred < 0.99)]
        
        self.ulozenie()
        
        return (self.x_H_od_t[-1], self.t_H, self.K)



    def ulozenie(self):
        if not os.path.exists("data_BD_"+self.rozmery+"_"+str(self.T)+"K"):
            os.makedirs("data_BD_"+self.rozmery+"_"+str(self.T)+"K")

        uloz("data_BD_"+self.rozmery+"_"+str(self.T)+"K/"+self.nazov+"_rhos_stred.npy", self.rhos_stred)






def colormesh_BD(rozmery,t,T=200,fps=50,smooth=False,uloz=True):

    r_moznosti, mC_moznosti, t_fin = np.load("data_BD_"+rozmery+"_"+str(T)+"K/config.npy",allow_pickle=True)
    
    #r_moznosti, mC_moznosti = np.load("data_BD_"+rozmery+"_"+str(T)+"K/config.npy") #docasne
    #t_fin = 30
    
    r_num = len(r_moznosti)
    mC_num = len(mC_moznosti)
    poc = int(t*fps)
    
    rhos_stred_po_t = [[0 for j in range(r_num)] for i in range(mC_num)]
    
    for i in range(r_num):
        for j in range(mC_num):
            nazov = "oblak_"+str(round_sig(r_moznosti[i]))+"cm_"+str(round_sig(mC_moznosti[j]))+"kg_"+str(T)+"K"     
            rhos_stred_po_t[j][i] = np.load("data_BD_"+rozmery+"_"+str(T)+"K/"+nazov+"_rhos_stred.npy")[poc]


    X,Y = np.meshgrid(np.array(r_moznosti),np.array(mC_moznosti))
    Z = np.array(rhos_stred_po_t)
    
    if smooth:
        levels = np.linspace(0, 1, 256)
    else:
        levels = np.linspace(0, 1, 11)

    fig, ax = plt.subplots(dpi=300)

    cs = ax.contourf(X,Y,Z, levels=levels)
    cbar = fig.colorbar(cs,ticks=np.linspace(0,1,6))

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$\\rho_r$ v strede, $t$ = "+str('{:.1f}'.format(round(poc/fps,2)))+" s, $T$ = "+str(T)+"K")

    #plt.text(0.6,0.8,"$t$ = "+str('{:.2f}'.format(round(poc/fps,2)))+" s", transform=fig.transFigure)

    if uloz:
        if smooth:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/colormesh_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/colormesh_"+rozmery+"_"+str(T)+"K"+".png")

    else:
        plt.show()



    # x_H pre t=max
    #Z = np.transpose(np.load("data_BD_"+rozmery+"_"+str(T)+"K/x_Hs.npy"))
    Z_x_H = np.load("data_BD_"+rozmery+"_"+str(T)+"K/x_Hs.npy")

    if smooth:
        levels_x_H = np.logspace(-2, 2, 256)
    else:
        levels_x_H = np.logspace(-2, 2, 13)
        
    fig, ax = plt.subplots(dpi=300)
    
    locator = ticker.LogLocator(base=10)
    cs = ax.contourf(X,Y,Z_x_H, locator=locator, levels=levels_x_H, cmap='inferno')
    cbar = fig.colorbar(cs,format='%.2f',ticks=np.logspace(-2, 2, 5))

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$x_H$ (cm) pre $t$ = "+str(t_fin)+" s, $T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/x_Hs_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/x_Hs_"+rozmery+"_"+str(T)+"K"+".png")
    else:
        plt.show()



    # x_H_relativne pre t=max
    #Z = np.transpose(np.load("data_BD_"+rozmery+"_"+str(T)+"K/x_Hs.npy"))/X
    Z_x_H_r = np.log10(1-np.load("data_BD_"+rozmery+"_"+str(T)+"K/x_Hs.npy")/X)
 
    fig, ax = plt.subplots(dpi=300)

    #locator = ticker.LogLocator(base=10)
    if smooth:
        cs = ax.contourf(X,Y,Z_x_H_r, levels=np.linspace(-2.25, 0, 256), cmap='inferno_r')
    else:
        cs = ax.contourf(X,Y,Z_x_H_r, levels=np.linspace(-2.25, 0, 10), cmap='inferno_r')
    
    cbar = fig.colorbar(cs,ticks=np.linspace(-2, 0, 5))
    cbar.ax.invert_yaxis()

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$\log(1-x_H/r)$ pre $t$ = "+str(t_fin)+" s, $T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/x_Hs_relativne_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/x_Hs_relativne_"+rozmery+"_"+str(T)+"K"+".png")
    else:
        plt.show()
        


    # t_H pre x=0
    #Z = np.transpose(np.load("data_BD_"+rozmery+"_"+str(T)+"K/t_Hs.npy"))
    Z_t_H = np.load("data_BD_"+rozmery+"_"+str(T)+"K/t_Hs.npy")

    if smooth:
        levels_t_H = np.linspace(0, t_fin, 256)
    else:
        levels_t_H = np.linspace(0, t_fin, 13)
        
    fig, ax = plt.subplots(dpi=300)
    
    cs = ax.contourf(X,Y,Z_t_H,levels=levels_t_H, cmap='copper')
    cbar = fig.colorbar(cs,ticks=np.linspace(0,t_fin, 7))

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$t_H$ (s) pre $x = 0$, $T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/t_Hs_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/t_Hs_"+rozmery+"_"+str(T)+"K"+".png")
    else:
        plt.show()



    # t_H pre x=0 dofitovane logaritmicky
    Ks = np.load("data_BD_"+rozmery+"_"+str(T)+"K/Ks.npy")
    
    Z_t_H_dofitovane = [[0 for j in range(r_num)] for i in range(mC_num)]
    for i in range(r_num):
        for j in range(mC_num):
            prvok = Z_t_H[j][i]
            if prvok == np.inf:
                Z_t_H_dofitovane[j][i] = - X[j][i]/Ks[j][i]
            else:
                Z_t_H_dofitovane[j][i] = Z_t_H[j][i]

    if smooth:
        levels_t_H_dofit = np.logspace(-2, 3.5, 256)
    else:
        levels_t_H_dofit = np.logspace(-2, 3.5, 12)
    
    fig, ax = plt.subplots(dpi=300)

    locator = ticker.LogLocator(base=10)
    cs = ax.contourf(X,Y,Z_t_H_dofitovane, locator=locator, levels=levels_t_H_dofit, cmap='copper')

    cbar = fig.colorbar(cs,format='%.2f',ticks=np.logspace(-2,3, 6))

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")
    plt.title("$t_H$ (s) pre $x = 0$, $T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/t_Hs_dofitovane_"+rozmery+"_"+str(T)+"K"+"_smooth.png")
        else:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/t_Hs_dofitovane_"+rozmery+"_"+str(T)+"K"+".png")
    else:
        plt.show()



    # x_H + t_H kombinovany
    fig,ax = plt.subplots(dpi=300)
    
    cs1 = ax.contourf(X,Y,Z_x_H, locator=locator, levels=levels_x_H, cmap='inferno')
    cbar1 = fig.colorbar(cs1,format='%.2f',ticks=np.logspace(-2, 2, 5),location='right',
                         shrink=0.72,anchor=(0,1),label = "$x_H$ (cm) pre $t$ = "+str(t_fin)+" s")

    cs2 = ax.contourf(X,Y,Z_t_H,levels=levels_t_H, cmap='copper')
    cbar2 = fig.colorbar(cs2,ticks=np.linspace(0,t_fin, 7), location='bottom', pad=0.13,
                         label = "$t_H$ (s) pre $x = 0$")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")

    fig.set_figwidth(8)
    fig.set_figheight(8)
    
    plt.title("$T$ = "+str(T)+"K")

    if uloz:
        if smooth:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/x_H_t_H_"+rozmery+"_"+str(T)+"K_smooth.png")
        else:
            plt.savefig("data_BD_"+rozmery+"_"+str(T)+"K"+"/x_H_t_H_"+rozmery+"_"+str(T)+"K.png")
    else:
        plt.show()



def colormesh_BD_anim(rozmery,T=200,fps=50,smooth=False,uloz=True):

    fig,ax = plt.subplots()
    if uloz: fig.set_dpi(300)
    camera = Camera(fig)

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("$r$ (cm)")
    plt.ylabel("$m_C$ (kg)")

    plt.title("$\\rho_r$ v strede, $T$ = "+str(T)+"K")
    
    if smooth:
        levels = np.linspace(0, 1, 256)
    else:
        levels = np.linspace(0, 1, 11)
        
    r_moznosti, mC_moznosti, t_fin = np.load("data_BD_"+rozmery+"_"+str(T)+"K/config.npy",allow_pickle=True)
    r_num = len(r_moznosti)
    mC_num = len(mC_moznosti)
    
    for poc in range(t_fin*fps):
        
        rhos_stred_po_t = [[0 for i in range(r_num)] for j in range(mC_num)]
        
        for i in range(r_num):
            for j in range(mC_num):
                nazov = "oblak_"+str(round_sig(r_moznosti[i]))+"cm_"+str(round_sig(mC_moznosti[j]))+"kg_"+str(T)+"K"     
                rhos_stred_po_t[j][i] = np.load("data_BD_"+rozmery+"_"+str(T)+"K/"+nazov+"_rhos_stred.npy")[poc]
        
        X,Y = np.meshgrid(np.array(r_moznosti),np.array(mC_moznosti))
        Z = np.array(rhos_stred_po_t)

        #fig.clear()
        cs = plt.contourf(X,Y,Z,levels=levels)
        #cbar = fig.colorbar(cs,ticks=[0,0.2,0.4,0.6,0.8,1])

        plt.text(0.7,0.8,"$t$ = "+str('{:.2f}'.format(round(poc/fps,2)))+" s", transform=fig.transFigure)

        camera.snap()


    ani = camera.animate(interval=int(1000/fps))
    
    if uloz:
        if smooth:
            ani.save("data_BD_"+rozmery+"_"+str(T)+"K"+"/colormesh_"+rozmery+"_"+str(T)+"K"+"_smooth.mp4")
        else:
            ani.save("data_BD_"+rozmery+"_"+str(T)+"K"+"/colormesh_anim_"+rozmery+"_"+str(T)+"K"+".mp4")

    else:
        plt.show()

