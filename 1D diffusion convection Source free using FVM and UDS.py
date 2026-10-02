# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 10:21:48 2026

@author: stan
"""
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

#python code solves the 1D convection diffusion eqn on a 1D domain 
#note that T represents a fluid property phi here to go faster we write T 
#parameters 
L=1
N=int(50) #number of control volumes in the domain 
h=L/N

Pe=5 #grid peclet number Pe = 5
Pe_L = (Pe/h)*L #global peclet number
print(f"Global peclet number = {Pe_L}")
#Initialise temperatures
T=np.zeros(N)
#initialise iterated temperatures 
T_new=T.copy()

#error threshold
nummerical_error = 1 
epsi=1.E-8
iteration= 0

#Coefficients 
aW_val=1+Pe
aE_val=1
aP_val=aE_val+aW_val

#convegence criteria

while nummerical_error > epsi:
    for i in range(0,N):
      if i == 0 :  #Give special treatment to first node
       aP = aE_val + (2 + Pe)
       T_new[i] = (aE_val*T[i+1]+ 2 + Pe)/aP
      elif i == N-1:
        aP = aW_val + 2 
        T_new[i] = (aW_val * T[i-1])/aP 
      else :
        T_new[i] = (aE_val * T[i+1] + aW_val * T[i-1])/aP_val
     #Reset numerical error 
    nummerical_error=0
    for i in range(0,N):
          nummerical_error = nummerical_error + np.abs(T_new[i]-T[i])
    iteration = iteration + 1
    T=T_new.copy()

print(f"Temperature distribution on 1-D domain takes {iteration} iterations to get there") 

#computing the analytical solution 
x = np.linspace(0, 1, 100)
x_dom = (np.arange(N) + 0.5) * h 
T_Analytical=(np.exp(Pe_L * x)-1)/(np.exp(Pe_L)-1)*-1  +1 
plt.plot(x_dom, T, "ro-", label="Numerical")
plt.plot(x, T_Analytical, "b-", label="Analytical")
plt.xlabel("x")
plt.ylabel("phi")
plt.title("Numerical vs Analytical Comparison of phi for Pe = 5")
plt.legend()
plt.grid()
plt.show()

#interpolation
T_Analytical_interp = np.interp(x_dom, x, T_Analytical)

#table 
df5=pd.DataFrame({
    "Node": np.arange(1,N+1),
    "Distance" :x_dom,
    "Numerical solution of phi " :T,
    "Analytical solution of phi ":T_Analytical_interp   
    })
df5_display=df5.round(6)
print(df5.to_string(index=False))
df5.to_csv("table for phi for Pe = 5 ", index=False)


        
