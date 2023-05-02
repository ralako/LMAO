import numpy as np
from math import pi, sqrt, log10
import ROS
import warnings

# BOS = Brzdenie Oblakom Simulácia

warnings.filterwarnings("ignore", message="invalid value encountered in double_scalars")


def Dv(mD,Dy,m,C,S,v,r,A,sigma):
    menovatel = 16 * mD * r*r*r
    citatel = C * S * v * m
    exponenciala = np.exp(-Dy*Dy / (sigma*sigma))
    return (citatel/menovatel) * A*A*A * exponenciala * sqrt(2*pi) * sigma


def prech_obl(mC):
    return mC**(-0.491) * 10**(-7.8)   # r (cm)



class Odpad:
    def __init__(self,mD,C,S,v):
        self.mD = mD
        self.C = C
        self.S = S
        self.v = v

    def prelet_stredom(self, mC_start,mC_stop,mC_num, m_start,m_stop,m_num, Dy=0,T=200,uloz=True):

        #init
        self.mC_start = mC_start
        self.mC_stop = mC_stop
        self.mC_num = mC_num
        self.m_start = m_start
        self.m_stop = m_stop
        self.m_num = m_num
        self.Dy = Dy
        self.T = T
        self.uloz = uloz
        
        # grid
        self.mC_moznosti = np.logspace(log10(self.mC_start), log10(self.mC_stop), num=self.mC_num)
        self.m_moznosti = np.logspace(log10(self.m_start), log10(self.m_stop), num=self.m_num)
        self.rozmery = str(self.m_num)+"x"+str(self.mC_num)+"_prelet"
        np.save("config_"+self.rozmery+"_Dy0.npy", np.array([self.m_moznosti, self.mC_moznosti,
                                                             self.mD, self.C, self.S, self.v],
                                                            dtype = object))
        self.X,self.Y = np.meshgrid(self.m_moznosti,self.mC_moznosti)

        # Delta v
        self.Dvs_Dy0 = [[0 for j in range(self.m_num)] for i in range(self.mC_num)]
        self.sigmy = [0 for j in range(self.m_num)]

        for i in range(self.mC_num):
            self.mC = self.mC_moznosti[i]
            self.r = prech_obl(self.mC)    # cm
            #print(mC,r)
            oblak = ROS.Oblak(self.r,self.mC,self.T, vysl_poc=11,
                              uloz=self.uloz,video=False,fituj=True,model="1D",rozmery=self.rozmery)
            oblak.spracovanie()
            
            self.A = oblak.Acka_vysl[10]
            if self.A == 0:
                self.A = ROS.hyperbola(10,oblak.B,oblak.C)
            self.sigma = oblak.sigmy_vysl[10]  # cm
            if self.sigma == 0:
                self.sigma = ROS.priamka(10,oblak.D,oblak.E)
                
            self.sigmy[i] = self.sigma
            
            for j in range(m_num):
                self.m = self.m_moznosti[j]
                #print(self.sigma,self.A,self.r,self.m)
                self.Dvs_Dy0[i][j] = Dv(self.mD,0.01*self.Dy,self.m,self.C,self.S,self.v,0.01*self.r,self.A,0.01*self.sigma)
            print("Done",ROS.round_sig(self.mC))

        # save
        np.save("Dvs_"+self.rozmery+"_Dy0.npy",self.Dvs_Dy0)
        np.save("sigmy_"+self.rozmery+"_Dy0.npy",self.sigmy)



    def prelet_bokom(self, Dy_start,Dy_stop,Dy_num, m_start,m_stop,m_num, mC=1e-16,T=200,uloz=True):

        # init
        self.Dy_start = Dy_start
        self.Dy_stop = Dy_stop
        self.Dy_num = Dy_num
        self.m_start = m_start
        self.m_stop = m_stop
        self.m_num = m_num
        self.mC = mC
        self.T = T
        self.uloz = uloz
        
        # grid
        self.Dy_moznosti = np.linspace(self.Dy_start, self.Dy_stop, num=self.Dy_num)
        self.m_moznosti = np.logspace(log10(self.m_start), log10(self.m_stop), num=self.m_num)
        self.rozmery = str(self.m_num)+"x"+str(self.Dy_num)+"_prelet_bokom"
        np.save("config_"+self.rozmery+"_mCfix.npy", np.array([self.m_moznosti, self.Dy_moznosti, self.mC],dtype=object))
        self.X,self.Y = np.meshgrid(self.m_moznosti,self.Dy_moznosti)

        # Delta v
        self.Dvs_mCfix = [[0 for j in range(self.m_num)] for i in range(self.Dy_num)]

        self.r = prech_obl(self.mC)    # cm
        #print(mC,r)
        oblak = ROS.Oblak(self.r,self.mC,self.T, vysl_poc=11,
                          uloz=self.uloz,video=False,fituj=True,model="1D",rozmery=self.rozmery)
        oblak.spracovanie()
        
        self.A = oblak.Acka_vysl[10]
        self.sigma = oblak.sigmy_vysl[10]  # cm
        
        for i in range(self.Dy_num):
            self.Dy = self.Dy_moznosti[i]
            for j in range(m_num):
                self.m = self.m_moznosti[j]
                self.Dvs_mCfix[i][j] = Dv(self.mD,0.01*self.Dy,self.m,self.C,self.S,self.v,0.01*self.r,self.A,0.01*self.sigma)


        # save
        np.save("Dvs_"+self.rozmery+"_mCfix.npy",self.Dvs_mCfix)

