#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:37:46 2016

@author: sdebuyl
"""

import matplotlib.pyplot as plt
import odetest
import os
import numpy as np

problem = odetest.Problem1()

u, t =odetest.solver(problem, np.ones(2, dtype=np.float64),np.arange(0,10,.01, dtype=np.float64), odetest.RK2)

tableTOSAVE=[[2,1,2],[3,5.5,2.22]]
cwd = os.getcwd()
filepath=''.join([cwd,'/tt.txt'])
np.savetxt(filepath,u)

print(len(u))
print(len(t))
print(u[0:10,1])



fig=plt.figure()
plt.plot(t,u[:,1],label=r"output/input1")#,label="computed points"
plt.show()
fig.savefig('full_figure.pdf')      
