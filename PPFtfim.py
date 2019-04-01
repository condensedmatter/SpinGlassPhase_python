# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 01:43:56 2018

Python Package For transverse field ising model (PPFtfim)

M represents the Hamiltonian


@author: Wang, Jian
"""

import numpy as np
from pf import pf   # the Python Package of Assessing Pfaffian (PPAP)
import time # for test 



class tfim(object):
    def __init__(self,M): 
        # each M gives quenched disorder Hamiltonian configuration, 
        # in case of disorder average, this object need to be constructed multiple times
        # with different Ms satisfying some distribution.
        self.M=M # M is an L by L matrix, contains all parameters of a chain
        self.L=M.shape[0]  # L is the length of 1 dimensional Ising chain
        self.u, self.s, self.vh = np.linalg.svd(M)
        
        self.AABBt0=self.aux_matrix_dynamics(0).real
        

   
    def aux_matrix_dynamics_block_for_test(self,t):
        # this is for test purpose
        A1=np.zeros((self.L, self.L), dtype=complex)
        A2=A1+0
        A3=A1+0
        A4=A1+0
        for i in range(self.L):
            for j in range(self.L):
                A1[i,j]=(i+1)+(j+1)*1j
                A2[i,j]=(i+1)-(j+1)*1j
                A3[i,j]=-(i+1)+(j+1)*1j
                A4[i,j]=-(i+1)-(j+1)*1j
                
        for i in range(self.L):
            for j in range(self.L):
                A1[i,j]=1+t
                A2[i,j]=2+t
                A3[i,j]=3+t
                A4[i,j]=4+t
        A12=np.concatenate((A1, A2), axis=1)
        A34=np.concatenate((A3, A4), axis=1)
        return np.concatenate((A12,A34),axis=0)   # 这里之所以没有乘以2，参见“推导3”
        
    def aux_matrix_dynamics(self,t):
        W=np.diag(np.exp(-2.0*1j*t*self.s))
        A1=np.matmul(np.matmul(self.vh.transpose(), W), self.vh)
        A2=np.matmul(np.matmul(self.vh.transpose(), W), self.u.transpose())
        A3=np.matmul(np.matmul(self.u, W), self.vh)
        A4=np.matmul(np.matmul(self.u, W), self.u.transpose())
        A12=np.concatenate((A1, A2), axis=1)
        A34=np.concatenate((A3, A4), axis=1)
        return np.concatenate((A12,-A34),axis=0)   # 这里之所以没有乘以2，参见“推导3”
    
    def aux_pfaffian_constructor_t(self,i,j,t):
        #  i,j are integer
        #  0<i<j<=L
        # 这里要特别注意 i,j和 数组的index相差1
        #  S is a 2(i+j-1) dimensional matrix
        S=np.zeros((2*(i+j-1),2*(i+j-1)), dtype=complex)  # DEBUG, dtype=complex can not miss, otherwise only cast to real part
        Mt=self.aux_matrix_dynamics(t)
        #M0=self.aux_matrix_dynamics(0.0)  # NEED TO BE IMPROVED
        M0=self.AABBt0        
        AtA=Mt[0:self.L,0:self.L]
        AtB=Mt[0:self.L,self.L:]
        BtA=Mt[self.L:,0:self.L]
        BtB=Mt[self.L:,self.L:]        
        AA=M0[0:self.L,0:self.L]
        AB=M0[0:self.L,self.L:]
        BA=M0[self.L:,0:self.L]
        BB=M0[self.L:,self.L:]
        #print(BB.shape)
        a=2*i-1
        b=2*(i+j-1)
     
        S[0:a:2,0:a:2]=AA[0:i,0:i]
        S[0:a:2,1:a:2]=AB[0:i,0:(i-1)]
        S[1:a:2,0:a:2]=BA[0:(i-1),0:i]
        S[1:a:2,1:a:2]=BB[0:(i-1),0:(i-1)]

        S[a:b:2,a:b:2]=AA[0:j,0:j]
        S[a:b:2,(a+1):b:2]=AB[0:j,0:(j-1)]
        S[(a+1):b:2,a:b:2]=BA[0:(j-1),0:j]
        S[(a+1):b:2,(a+1):b:2]=BB[0:(j-1),0:(j-1)]
        
        S[0:a:2,a:b:2]=AtA[0:i,0:j]
        S[0:a:2,(a+1):b:2]=AtB[0:i,0:(j-1)]
        S[1:a:2,a:b:2]=BtA[0:(i-1),0:j]
        S[1:a:2,(a+1):b:2]=BtB[0:(i-1),0:(j-1)]
        
        for inew in range(b):
            for jnew in range(b):
                if inew>=jnew:
                    S[inew,jnew]=0
        S=S-S.transpose()
        return S


    def aux_pfaffian_constructor_t_AABB(self,i,j,t,AABB):
        S=np.zeros((2*(i+j-1),2*(i+j-1)), dtype=complex)  # DEBUG, dtype=complex can not miss, otherwise only cast to real part

        Mt=AABB[t,:,:]
        #M0=self.aux_matrix_dynamics(0.0)  # NEED TO BE IMPROVED
        M0=self.AABBt0
        
        AtA=Mt[0:self.L,0:self.L]
        AtB=Mt[0:self.L,self.L:]
        BtA=Mt[self.L:,0:self.L]
        BtB=Mt[self.L:,self.L:]
        
        AA=M0[0:self.L,0:self.L]
        AB=M0[0:self.L,self.L:]
        BA=M0[self.L:,0:self.L]
        BB=M0[self.L:,self.L:]
        #print(BB.shape)

        a=2*i-1
        b=2*(i+j-1)
     
        S[0:a:2,0:a:2]=AA[0:i,0:i]
        S[0:a:2,1:a:2]=AB[0:i,0:(i-1)]
        S[1:a:2,0:a:2]=BA[0:(i-1),0:i]
        S[1:a:2,1:a:2]=BB[0:(i-1),0:(i-1)]

        S[a:b:2,a:b:2]=AA[0:j,0:j]
        S[a:b:2,(a+1):b:2]=AB[0:j,0:(j-1)]
        S[(a+1):b:2,a:b:2]=BA[0:(j-1),0:j]
        S[(a+1):b:2,(a+1):b:2]=BB[0:(j-1),0:(j-1)]
        
        S[0:a:2,a:b:2]=AtA[0:i,0:j]
        S[0:a:2,(a+1):b:2]=AtB[0:i,0:(j-1)]
        S[1:a:2,a:b:2]=BtA[0:(i-1),0:j]
        S[1:a:2,(a+1):b:2]=BtB[0:(i-1),0:(j-1)]
        
        for inew in range(b):
            for jnew in range(b):
                if inew>=jnew:
                    S[inew,jnew]=0
        S=S-S.transpose()
        return S


    def aux_pfaffian_constructor(self,i,j):
        S=np.zeros((2*(i+j-1),2*(i+j-1)))  
        M0=self.AABBt0
        AA=M0[0:self.L,0:self.L]
        AB=M0[0:self.L,self.L:]
        BA=M0[self.L:,0:self.L]
        BB=M0[self.L:,self.L:]
        a=2*i-1
        b=2*(i+j-1)     
        S[0:a:2,0:a:2]=AA[0:i,0:i]
        S[0:a:2,1:a:2]=AB[0:i,0:(i-1)]
        S[1:a:2,0:a:2]=BA[0:(i-1),0:i]
        S[1:a:2,1:a:2]=BB[0:(i-1),0:(i-1)]
        S[a:b:2,a:b:2]=AA[0:j,0:j]
        S[a:b:2,(a+1):b:2]=AB[0:j,0:(j-1)]
        S[(a+1):b:2,a:b:2]=BA[0:(j-1),0:j]
        S[(a+1):b:2,(a+1):b:2]=BB[0:(j-1),0:(j-1)]        
        S[0:a:2,a:b:2]=AA[0:i,0:j]
        S[0:a:2,(a+1):b:2]=AB[0:i,0:(j-1)]
        S[1:a:2,a:b:2]=BA[0:(i-1),0:j]
        S[1:a:2,(a+1):b:2]=BB[0:(i-1),0:(j-1)]
        
        for inew in range(b):
            for jnew in range(b):
                if inew>=jnew:
                    S[inew,jnew]=0
        S=S-S.transpose()
        return S
        
    def correlator_equal_time(self,i,j):
        S=self.aux_pfaffian_constructor(i,j)
        return np.real(pf(S))
    
    def correlator_dynamics(self,i,t,j):
        S=self.aux_pfaffian_constructor_t(i,j,t)
        return pf(S)

    def correlator_dynamics_AABB(self,i,t,j,AABB):
        S=self.aux_pfaffian_constructor_t_AABB(i,j,t,AABB)
        return pf(S)
    
    def correlator_equal_time_Matrix(self):
        output=np.zeros((self.L, self.L))
        for i in range(self.L):
            for j in range(i,self.L):
                output[i,j]=self.correlator_equal_time(i+1,j+1)
                # attension (1) i,j runs over all sites, i,j is not the ith or jth site
                # attension (2) equal time correlator suppose to be real
                
        for i in range(self.L):
            for j in range(i):
                output[i,j]=output[j,i]
        return output
    
    def correlator_dynamics_sector(self,i,j,dt,tSteps):
        # return a complex matrix of dimension (tSteps,(i-j+1))
        # 1<=i<j<=self.L
        start = time.time()
        sector=np.zeros((tSteps,(j-i+1)), dtype=complex)  # DEBUG, dtype=complex can not miss, otherwise only cast to real part
        for r in range(0,(j-i+1)):
            for t in range(0,tSteps):
                sector[t,r]=self.correlator_dynamics(i,t*dt,i+r)
        end = time.time()
        print("time spent in sector")
        print(end - start)

        return sector
        
    
    def correlator_dynamics_sector_AABB(self,i,j,dt,tSteps):
        # return a complex matrix of dimension (tSteps,(i-j+1))
        # 1<=i<j<=self.L
        start = time.time()
        AABB=np.zeros((tSteps,2*self.L,2*self.L), dtype=complex)
        for t in range(0,tSteps):
            AABB[t,:,:]=self.aux_matrix_dynamics(t*dt)
        sector=np.zeros((tSteps,(j-i+1)), dtype=complex)  # DEBUG, dtype=complex can not miss, otherwise only cast to real part
        for r in range(0,(j-i+1)):
            for t in range(0,tSteps):
                sector[t,r]=self.correlator_dynamics_AABB(i,t,i+r,AABB)
        end = time.time()
        print("time spent in sector")
        print(end - start)

        return sector
    
    def hanningWindow2D(self,L1,L2):
        # the center is rolled to hh(0,0)
        h1=np.hanning(2*L1-1)
        h2=np.hanning(2*L2-1)
        hh=np.zeros((2*L1-1,2*L2-1))
        for i in range(2*L1-1):
            for j in range(2*L2-1):
                hh[i,j]=h1[i]*h2[j]
        hh=np.roll(hh,L1,axis=0)        
        hh=np.roll(hh,L2,axis=1) 
        return hh
    
    def Swk(self,i,j,dt,tSteps):
        Cnt=self.correlator_dynamics_sector(i,j,dt,tSteps)
        R=(j-i+1)
        CCCCnt=np.zeros((tSteps*2-1,R*2-1), dtype=complex)
        CCCCnt[0:tSteps,0:R]=Cnt
        CCCCnt[0:tSteps,R:2*R-1]=Cnt[:,1:R][:,::-1]
        CCCCnt[tSteps:2*tSteps-1,0:R]=np.conjugate(Cnt[1:tSteps,:][::-1,:])
        CCCCnt[tSteps:2*tSteps-1,R:2*R-1]=np.conjugate(Cnt[1:tSteps,1:R][::-1,::-1])

        hh=self.hanningWindow2D(tSteps,R)
        output=np.fft.fft2(np.multiply(CCCCnt,hh))
        
        return output.real

        
    
    
    
    def orderParameter_1(self):
        m=self.correlator_equal_time_Matrix()
        output=0
        for i in range(self.L):
            for j in range(self.L):
                output=output+m[i,j]
        return output
    
    def orderParameter_2(self):
        # this is the Edward-Anderson order parametr for a single quenced configuration
        # disorder average need to be done, for better convergence
        m=self.correlator_equal_time_Matrix()
        output=0
        for i in range(self.L):
            for j in range(self.L):
                output=output+abs(m[i,j])**2
        return output
                
        
