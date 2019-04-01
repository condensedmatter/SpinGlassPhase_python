# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 02:07:53 2018

@author: wangj
"""

from hamiltonian import *
from PPFtfim import *
import matplotlib.pyplot as plt


M=simplePureHamiltonian(20,1,1,-0.3)
plt.matshow(M)
plt.colorbar()
plt.show()

t=tfim(M)

print(t.correlator_equal_time(20,1))

mat=t.correlator_equal_time_Matrix()
plt.matshow(mat)
plt.colorbar()
plt.show()