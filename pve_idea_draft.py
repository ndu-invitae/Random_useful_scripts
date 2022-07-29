#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:07:14 2022

@author: niudu


Fitting Platt equation 
"""

import numpy as np
from scipy import stats
from scipy.optimize import least_squares
from typing import List, Union
Vector = list[Union[int, float]]


def flat_list(bulk_list):
    """
    Flattern lists within a list.
    """
    return [i for sublist in bulk_list for i in sublist]

def platt_equation_w_inhibition(coef,x,y):

    return coef[0]*(1-np.exp(-(coef[1]*x)/coef[0])) - y 

def platt_equation_wo_inhibition(coef,x,y):

    return coef[0]*(1-np.exp(-(coef[1]*x)/coef[0]))*np.exp((-1*coef[2]*x)/coef[0]) - y


def r_squared(yhat:np.array, y:np.array) -> float:
    SSE = sum(np.array(yhat - y)**2)
    SST = sum(np.array(yhat - np.mean(y))**2)
    return 1-(SSE/SST)


def fit_model(X:List[int,float], Y:List[int,float])-> np.array:
    """
    
    :param X: Independent variable of selected sample dataset
    :param Y: Dependent variable of selected sample dataset
    :return: fitting coefficients
    """
    def div0(a, b):
        if b == 0:
            return 1.
        return a/b
    
    X = np.array(X, dtype = float)
    Y = np.array(Y, dtype = float)
    assert len(X) == len(Y) >= 3, "X and Y length must be equal with at least 3 data points"
    
    init_ymax = max(Y)
    init_slope = div0(Y[1]-Y[0]),(X[1]-X[0])
    res_lsq = least_squares(platt_equation_w_inhibition, [init_ymax, init_slope], loss='soft_l1', bounds=(0,np.inf),  args=(X, Y)) 






def PvE_func(Light_level,P_data_ori):
    # E_prim and E_prim_culture refers to the E' of the sample in pH-MIMS system and culture. The values are not used in this function
    # Here they are simply transfferred to the data framework that will be referened in other functions
    try:
        scale_up = 10**int(-np.log10(np.median(P_data_ori)/100))
    except ValueError:
        scale_up = 10**10

    # Check if the first data point is negative
    offset = 0
    if min(P_data_ori)<0:
        offset = -1*min(P_data_ori)*scale_up
        P_data = P_data_ori*scale_up + offset
    else: P_data = P_data_ori*scale_up
    G_alpha = stats.linregress(Light_level[0:4],P_data[0:4]).slope
    G_alpha = 10**(-10) if (G_alpha < 0.0) else G_alpha
    # G_beta = 0
    G_Ps = max(P_data)
    # The non-linear regression funciton in scipy does not work well on extremely small numbers, therefore must be scaled up
    # Scale = -round(np.log10(G_alpha))
    
    coef_0 = np.array([G_Ps,G_alpha], dtype=float)
    
    # Start non-linear fitting
    res_lsq = least_squares(Platt, coef_0, loss='soft_l1', bounds=(0,np.inf), args=(Light_level, P_data)) 
    coef = res_lsq.x/scale_up
    #Scale back

    # The calculation of P shall be done outside this function
    # E0 = pd.DataFrame(list(np.arange(0,max(Light_level)*1.01,max(Light_level)/100)))  
    # Prod = Platt(coef[0:2],E0,offset/scale_up)
    Ek = coef[0]/coef[1]
   
    # Pm = coef[0]*(coef[1]/(coef[1]+coef[2]))*((coef[2]/(coef[1]+coef[2]))**(coef[2]/coef[1]))
    
    return list(coef)+[Ek]+[offset/scale_up]

def R_squared(X, Y):

    SSE = sum(np.array(X - Y)**2)
    average = np.mean(Y)
    SST = sum(np.array(X - average)**2)
    return 1-(SSE/SST)