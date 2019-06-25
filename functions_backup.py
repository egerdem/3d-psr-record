# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 21:27:37 2019

@author: ABRA A5
"""

"""
@author: Ege Erdem

3D PSR Extrapolation - Recording - FUNCTIONS

"""
import time
from scipy.io.wavfile import read
import numpy as np
from scipy.special import sph_harm, spherical_jn
from scipy import signal
import pickle
from convertangle import cart2sph, cart2sph_single, cart2sphr, sph2cart
from sphericals import spherical_bn

pi = np.pi

def interpolation(N, k, ti, fi):
    print("eski kod, backup, yavaş")

    pressure_grad = []
    velocity = []
    pressure = []
    I_active = []
    I_n = []
#    if np.isscalar(rc):
#        rc = np.full((len(teta_ls),1)  , rc)
        
    for i in range(len(teta_ls)):
        
        p_grad = grad(teta_ls[i], phi_ls[i], k, N, ti, fi,i) # grad of the pressure        
        u = calc_u(p_grad, ro, k, f[fi])  # local velocity calculation using euler momentum equation with pressure grad        
        p = pressure_field(teta_ls[i], phi_ls[i],k, fi, ti,i,N)  # extrapolated pressure at r (due to monochromatic plane wave coming from mt,mp)

        I_act, In = ia(p, u)  # local active intensity
#        print("I_act =", I_act)
#        print("In =", In)

        pressure_grad.append(p_grad)
        velocity.append(u)
        pressure.append(p)
        I_active.append(I_act)
        I_n.append(In)
    
    pressure_grad = np.array(pressure_grad)
    velocity = np.array(pressure_grad)
    pressure = np.array(pressure)
    I_active = np.array(I_active)
    
    I_n = np.array(I_n)
#    If,Ig,Ih = list(zip(*I_active))
#    Iu = np.array(If)    
#    Iv = np.array(Ig)  
#    Iw = np.array(Ih)
    return(pressure,I_active)


def grad(teta_ls, phi_ls, k, N, ti, fi,i):
    teta = teta_ls
    phi = phi_ls
    """ pressure grad in spherical coordinates, to be used for obtaining velocity and then intensity """
    r_unit = np.array([np.sin(teta)*np.cos(phi), np.sin(teta)*np.sin(phi) , np.cos(teta)])

    p_unit = np.array([-np.sin(phi)*np.sin(teta), np.sin(teta)*np.cos(phi) , 0])

    t_unit = np.array([np.cos(teta)*np.cos(phi), np.cos(teta)*np.sin(phi) , -np.sin(teta)])
    
    counter = 0
    gradp = 0.
    gradpX = 0.
    pr_total = 0.
    pteta_der_total = 0.
    pphi_der_total = 0
    A= np.ones((10,2))*1j
    for n in range(N+1):
        for m in range(-n,n+1):
                    
            pr_der = anm_offline[counter, fi,ti]*4*pi*((1j)**n)*k*spherical_jn_der_offline[fi,n]*sph_harm(m, n, phi, teta)
            pr_total = pr_total + pr_der
            Y = sph_harm(m+1, n, phi, teta)
            if np.isnan(Y):               
                Y = 0.
            pteta_der = anm_offline[counter, fi,ti]*4*pi*((1j)**n)*spherical_jn_offline[fi,n]*(m*(1/np.tan(teta))*sph_harm(m, n, phi, teta)
            + np.sqrt((n-m)*(n+m+1))*np.exp(-1j*phi)*Y)    
            
            pteta_der_total = pteta_der_total + pteta_der
            
            pphi_der = anm_offline[counter, fi,ti]*4*pi*((1j)**n)*spherical_jn_offline[fi,n]*1j*m*sph_harm(m, n, phi, teta)
            pphi_der_total = pphi_der_total + pphi_der 


            grad_sum = pr_der*r_unit+(1/rc)*pteta_der*t_unit + (1/(rc*np.sin(teta)))*pphi_der*p_unit
            gradp = gradp + grad_sum
            counter = counter +1
            
    A[i] = pr_total
        
#    print("pr_total backup",pr_total)
#    print("pteta_der_total backup",pteta_der_total)
#    print("pphi_der_total backup",pphi_der_total)

    
    grad_sumX = pr_total*r_unit+(1/rc)*pteta_der_total*t_unit + (1/(rc*np.sin(teta)))*pphi_der_total*p_unit
    
#    print("gradp backup = ", gradp)
#    print("grad_sumXX = ", grad_sumX,"\n")
#    print("gradp=",gradp)

    return(gradp)
    
def a_nm(n, fi, ti, counter):
    p_nm = pq_Ymn[counter, fi, ti]                
    anm[counter, fi, ti] = p_nm / spherical_bn_offline[fi,n]  
    
    return(anm[counter, fi, ti])

    
def calc_u(pressure_grad, ro, k,f):
    w = 2*np.pi*f
    u = -pressure_grad/(1j*w*ro)
    
    return(u)


def pressure_field(teta_ls, phi_ls, k, fi, ti,i,N): 
    """ extrapolated pressure at distance r due to a monochromatic plane wave coming from mt,mp (teta,azimuth) """

    pr = 0
    counter = 0
    
#        pressure_ls, I_ls = interpolation(rc, teta_ls, phi_ls, N, k, f[fi],counter)
    
    for n in range(N+1):
        for m in range(-n,n+1):  

            term = anm_offline[counter, fi,ti]*4*np.pi*((1j)**n)*spherical_jn_offline[fi,n]*sph_harm_offline[counter,i]
            
            pr = pr + term
            counter = counter + 1
    
#    print("pr backup =", pr)

#    print(counter, "coefflik toplandı")
    
    return(pr)    
    
def ia(p, u):
    """ local active intensity with known pressure and velocity """
    I_a = 0.5*(p*np.conj(u)).real
#    print(p,u)
    In = I_a/np.linalg.norm(I_a)
    return(I_a, In)
    
def g_tid(teta):
    
    g = 0.00137175+0.457839*np.cos(teta)+0.535779*(np.cos(teta)**2) +0.0401825*(np.cos(teta)**3) - 0.126105*(np.cos(teta)**4)

    return g

def ls_signals(I_n_ls, pressure_ls):


    angle_psi_sh = []
    angle_psi_degree = []
    gains = []
    
    
    for j in range(len(I_n_ls)):
        """ t: angle between local active intensity and ls vectors"""
        t = np.arccos(np.sum(I_n_ls[j]*-L_ra_unit[j])/(np.linalg.norm(I_n_ls[j])*np.linalg.norm(L_ra_unit[j])))
        angle_psi_sh.append(t)
        angle_psi_degree.append(np.degrees(t))      
        
        gains.append(g_tid(t))
        
    prsum = gains*pressure_ls.T

    return(prsum.T)