# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 03:27:50 2018

the goal is to check the two order p


@author: wangj

"""

from hamiltonian import *
from PPFtfim import *
import matplotlib.pyplot as plt

def run(L):
    # L=8,10,12,14,16
    steps=10
    
    
    
    disOrder=np.zeros(steps)
    RR=np.zeros(steps)
    
    SG_mean=np.zeros(steps)
    SG_std=np.zeros(steps)
    
    
    for i in range(steps):
        RR[i]=10
        disOrder[i]=0.25*1.8**i
    
    for i in range(steps):
        print(i)
        output=np.zeros(int(RR[i]))
        for j in range(int(RR[i])):
            M=DisorderedHamiltonian(L,disOrder[i])
            t=tfim(M)
            output[j]=t.orderParameter_2()/L**2
        SG_mean[i]=np.mean(output)
        SG_std[i]=np.std(output)
        
    plt.semilogx(disOrder,SG_mean,'o-',disOrder,SG_std,'.-')
    plt.show()
    
    
    
    RR2=np.zeros(steps)
    
    SG_mean2=np.zeros(steps)
    
    
    for i in range(steps):
        RR2[i]=(SG_std[i]/1)**2+1
    print(RR2)
    
    for i in range(steps):
        print(i)
        output=np.zeros(int(RR2[i]))
        for j in range(int(RR2[i])):
            M=DisorderedHamiltonian(L,disOrder[i])
            t=tfim(M)
            output[j]=t.orderParameter_2()/L**2
        SG_mean2[i]=np.mean(output)
        SG_mean2[i]=(SG_mean2[i]*RR2[i]+SG_mean[i]*RR[i])/(RR2[i]+RR[i])
        
    plt.semilogx(disOrder,SG_mean2,'o-',disOrder,SG_std,'.-')
    plt.show()
    
    np.save("data2Log"+str(L), SG_mean2)


for i in [32,64,128]:
    print("*********************** Now the size is  *********************")
    print("                              "+str(i))
    print("*********************** Now the size is  *********************")

    run(i)