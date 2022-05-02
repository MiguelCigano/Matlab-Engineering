#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 18:04:59 2022

@author: jmiguel
"""

import numpy as np 
from control.matlab import*
import matplotlib.pyplot as plt

A = np.array([[0, 1],[-0.16, -1 ]]);
B = np.array([[0],[1]]);
K = np.array([[0.34, -2]]);
Q = np.eye(2, dtype= int);
R = np.array([[2]]);
T = 0;
N = 19;

x = [None]*(N+1)
X = np.array([[10],[-10]])
x[0] = X;


for i in range(N):
      x[i+1] = A@x[i]-B@K@x[i]
      
      

vec = np.zeros((len(x),1));
for i in range(20):
    vec_x1 = x[i];
    vec[i] = vec_x1[0];
    

vec2 = np.zeros((len(x), 1));
for i in range(20):
    vec_x2 = x[i];
    vec2[i] = vec_x2[1]
    
    
plt.plot(vec)
plt.plot(vec2)
plt.grid()
plt.show()
