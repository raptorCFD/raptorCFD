#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def getThickParamIGP(xt,T,rho_bar,beta_bar):

    rho = rho_bar * (T/xt)**2
    t1 = np.sqrt(2 * rho)
    beta = beta_bar * np.arctan(T/(1-xt))

    invA =  np.array([[-(2*(2*xt - 1))/(- xt**4 + 3*xt**3 - 3*xt**2 + xt),           
     -1/(xt**2 - 2*xt + 1),        
      -(2*xt**2)/(xt**2 - 2*xt + 1),  
       -(2*(- xt**3 + 2*xt**2))/(xt**3 - 3*xt**2 + 3*xt - 1)],\
       
                  [-(- 8*xt**2 + xt + 1)/(xt*(- xt**4 + 3*xt**3 - 3*xt**2 + xt)), 
                  (2*xt + 1)/(xt**3 - 2*xt**2 + xt), 
                  (2*(xt**2 + 2*xt))/(xt**2 - 2*xt + 1),    
                  -(xt**3 + xt**2 - 8*xt)/(xt**3 - 3*xt**2 + 3*xt - 1)],\
                  
              [-(2*(2*xt**2 + 2*xt - 1))/(xt*(- xt**4 + 3*xt**3 - 3*xt**2 + xt)), 
               -(xt + 2)/(xt**3 - 2*xt**2 + xt), 
                 -(2*(2*xt + 1))/(xt**2 - 2*xt + 1), 
                 -(2*(- xt**2 + 2*xt + 2))/(xt**3 - 3*xt**2 + 3*xt - 1)],\
                 
                           [ (3*xt - 1)/(xt*(- xt**4 + 3*xt**3 - 3*xt**2 + xt)),     
                                1/(xt**3 - 2*xt**2 + xt),           
                                      2/(xt**2 - 2*xt + 1),           
                                           -(xt - 3)/(xt**3 - 3*xt**2 + 3*xt - 1)]])

    
    b = np.array([T - xt ** 0.5 * t1, -0.5 * t1 * xt ** -0.5, -np.tan(beta/2) - 0.25 * t1, -t1])
    
    t = invA.dot(b)

    t2 = t[0]
    t3 = t[1]
    t4 = t[2]
    t5 = t[3]

    return t1,t2,t3,t4,t5
