# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 21:15:00 2018

@author: wangj
"""
L1=3
L2=5
h1=np.hanning(2*L1-1)
h2=np.hanning(2*L2-1)
hh=np.zeros((2*L1-1,2*L2-1))
for i in range(2*L1-1):
    for j in range(2*L2-1):
        hh[i,j]=h1[i]*h2[j]
        
        
hh=np.roll(hh,L1,axis=0)        
hh=np.roll(hh,L2,axis=1)        

plt.imshow(hh)



