# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 23:56:57 2018

@author: wangj
"""


import numpy as np
import h5py

class Bank:
    def __init__(self,filename):
        #如果文件不存在 就创建
        #如果文件存在 就打开
        self.f=h5py.File(filename, "a")
        
        self.CurrentID=0
        
        if(self.f.__contains__("SizeTable")):
            self.SizeTable=self.f.get("SizeTable")
            print("exist")
        else:
            self.SizeTable=self.f.create_dataset("SizeTable",(0,2),maxshape=(None,2))
            print("new")
            
            
    def add_new_size(self,L):
        # 警告，只有在没有L的时候追加
        IDname=str(self.CurrentID)
        self.SizeTable.resize(self.SizeTable.shape[0]+1,axis=0)
        print(self.SizeTable.shape)

        self.SizeTable[-1,0]=L
        self.SizeTable[-1,1]=IDname

        temp=self.f.create_group(IDname)
        temp.create_dataset("ParameterTable",(0,5),maxshape=(None,5))
        self.CurrentID=self.CurrentID+1
        
        
            
    #def store(self,L,dJ,chi):
        

        
        
    def close(self):
        self.f.close()
        
        
        

T=Bank("bankdata.h5")
#T.add_new_size(3)
T.close()