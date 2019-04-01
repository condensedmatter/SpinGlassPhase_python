# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 03:27:50 2018

the goal is to check the two order p


@author: wangj

"""

from hamiltonian import *
from PPFtfim import *
import matplotlib.pyplot as plt

steps=40
dlam2=0.25
L=10

lam2=np.zeros(steps)
OP1=np.zeros(steps)
OP2=np.zeros(steps)

for i in range(steps):
    lam2[i]=(-5)+dlam2*i

for i in range(steps):
    print(i)
    M=simplePureHamiltonian(L,1,lam2[i],1)
    t=tfim(M)
    OP1[i]=t.orderParameter_1()/L**2
    OP2[i]=t.orderParameter_2()/L**2
    
plt.plot(lam2,OP2,'o-',lam2,OP1,'.-')
plt.show()