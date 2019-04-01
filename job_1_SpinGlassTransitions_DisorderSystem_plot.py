# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 01:03:40 2018

@author: jian
"""
import matplotlib.pyplot as plt

for i in range(10):
    disOrder[i]=0.25*1.8**i
print(disOrder)


for i in [8,10,12,14,16]:
    s=np.load("data2Log"+str(i)+".npy")*i**1
    plt.semilogx(disOrder,s,'.-')

plt.xlabel("disorder strength $\delta J$",fontsize=20)
plt.ylabel("$ \chi_{SG} /L^{  } $",fontsize=20)
plt.legend(["8","10","12","14","16"],loc=4)

plt.show()



for i in [8,10,12,14,16]:
    s=np.load("data2Log"+str(i)+".npy")*i**0.02
    plt.plot(disOrder[4:6],s[4:6],'.-')

plt.xlabel("disorder strength $\delta J$",fontsize=20)
plt.ylabel("$ \chi_{SG} /L^{2} $",fontsize=20)
plt.legend(["8","10","12","14","16"],loc=4)

plt.show()