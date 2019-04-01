# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 14:35:30 2018

@author: wangj
"""

from hamiltonian import simplePureHamiltonian
from PPFtfim import tfim
import matplotlib.pyplot as plt
import numpy as np

L=16
STEPS=80
dt=0.125

M=simplePureHamiltonian(4*L,1,1.6, 0)
#plt.matshow(M)
#plt.colorbar()
#plt.show()

t=tfim(M)

AABB=t.AABBt0
#plt.matshow(AABB)
#plt.colorbar()
#plt.show()

swk=t.Swk(L,3*L,dt,STEPS)
cnt=t.correlator_dynamics_sector(L,3*L,dt,STEPS)

plt.imshow(np.abs(cnt))
plt.show()

plt.matshow(swk.real)
plt.colorbar()
plt.matshow(np.log(np.abs(swk.real)+np.e*0.02))
plt.colorbar
plt.show()
print(np.min(swk))
'''
i=L
j=3*L
tSteps=STEPS
sec1=t.correlator_dynamics_sector_AABB(i,j,dt,tSteps)

sec2=t.correlator_dynamics_sector(i,j,dt,tSteps)
print(np.max(np.abs(sec1-sec2)))
'''