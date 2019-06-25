# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 02:23:42 2018

@author: ABRA A5
"""
from scipy.special import sph_harm, lpmv, spherical_jn, spherical_yn
from scipy.special import sph_harm
import numpy as np

#jn, jnp, yn , ynp = sph_jnyn(n, kr)
#hn = jn - 1j*yn
#hnp = jnp - 1j*ynp

def spherical_hn2(n,kr):    
    
    hn2 = spherical_jn(n, kr, derivative=False) - 1j*spherical_yn(n, kr, derivative=False)
    
    return(hn2)
    
def spherical_hn2der(n,kr):
    
#    hn2der = 0.5*(spherical_hn2(n-1,kr)- (spherical_hn2(n,kr)+kr*spherical_hn2(n+1,kr)/kr))   
    hn2der = spherical_jn(n, kr, derivative=True) - 1j*spherical_yn(n, kr, derivative=True)
    
    
    return(hn2der)

def spherical_bn(n,k,r,ra):
    
    kr = k*r
    kra = k*ra
    bn = spherical_jn(n, kr, derivative=False) - spherical_jn(n, kra, derivative=True)*(spherical_hn2(n,kr))/spherical_hn2der(n,kra)
#    bn = spherical_jn(n, kr, derivative=False) - spherical_jn(n, kra, derivative=True)/spherical_hn2der(n,kra)

    return(bn)
  
#rc = 0.155
#c = 343.0 
#ra = 0.12


#N = 4
#cc = 0
#zz = 0
#dd = 0


"""
N=4
sperical_harm_j = np.zeros((25,10))*1j

for i in range(10):
    for n in range(N+1):
        for m in range(-n,n+1):  
            

            sperical_harm_j[zz,i] = sph_harm(m, n, phi_ls[i], teta_ls[i])*((1j)**n)
            zz += 1
    zz = 0

dosya = open("sph_harm_j","wb")
pickle.dump(sperical_harm_j, dosya)
dosya.close()

sperical_harm_der = np.zeros((25,10))*1j

for i in range(10):
    for n in range(N+1):
        for m in range(-n,n+1):  
            
            Y = sph_harm(m+1, n, phi_ls[i], teta_ls[i])
            if np.isnan(Y):               
                Y = 0.
            sperical_harm_der[cc,i] = (m*(1/np.tan(teta_ls[i]))*sph_harm(m, n, phi_ls[i], teta_ls[i])
            + np.sqrt((n-m)*(n+m+1))*np.exp(-1j*phi_ls[i])*Y) * ((1j)**n)
            cc += 1
    cc = 0

dosya = open("sph_harm_offline_der","wb")
pickle.dump(sperical_harm_der, dosya)
dosya.close()

    
  
        
sperical_harm_m = np.zeros((25,10))*1j

for i in range(10):
    for n in range(N+1):
        for m in range(-n,n+1):  
            

            sperical_harm_m[dd,i] = sph_harm(m, n, phi_ls[i], teta_ls[i])*m *((1j)**n)
            dd += 1
    dd = 0

dosya = open("sph_harm_offline_m","wb")
pickle.dump(sperical_harm_m, dosya)
dosya.close()


sperical_harm_ege = np.zeros((25,10))*1j
ee = 0
for i in range(10):
    for n in range(N+1):
        for m in range(-n,n+1):  
            
            sperical_harm_ege[ee,i] = sph_harm(m, n, phi_ls[i], teta_ls[i])*((1j)**n)
            ee += 1
    ee = 0
""" 

#dosya = open("sph_harm_offline","wb")
#pickle.dump(sperical_harm, dosya)
#dosya.close()

#spherical_jn_offline = np.zeros((513,N+1))*1j
#spherical_jn_der_offline = np.zeros((513,N+1))*1j

"""
#spherical_bn_offline = np.zeros((513,N+1))*1j

N=4   

for n in range(N+1):
    print(n)
    for fi in range(513):
        k = 2*np.pi*f[fi]/c
        spherical_bn_offline[fi, n] = spherical_bn(n,k,rc,ra)
#        spherical_jn_offline[fi, n] = spherical_jn(n, k*rc, derivative=False)
#        spherical_jn_der_offline[fi, n] = spherical_jn(n, k*rc, derivative=True)
#
#
spherical_bn_offline[0,:] = 0

dosya = open("spherical_bn_offline","wb")
pickle.dump(spherical_bn_offline, dosya)
dosya.close()   

"""

       
#
#dosya = open("spherical_jn_offline","wb")
#pickle.dump(spherical_jn_offline, dosya)
#dosya.close()  
#    
#dosya = open("spherical_jn_der_offline","wb")
#pickle.dump(spherical_jn_der_offline, dosya)
#dosya.close()   
#        
#        

"""
spherical_jn_offline_ext = np.zeros((513,25))*1j
spherical_jn_der_offline_ext = np.zeros((513,25))*1j

for fi in range(len(f)):
    
    for i in range(25):
        if i==0:
            spherical_jn_offline_ext[fi,i] = spherical_jn_offline[fi, 0]
        if i>0 and i<=3:
            spherical_jn_offline_ext[fi,i] = spherical_jn_offline[fi, 1]
        if i>3 and i<=8:
            spherical_jn_offline_ext[fi,i] = spherical_jn_offline[fi, 2]
        if i>8 and i<=15:
            spherical_jn_offline_ext[fi,i] = spherical_jn_offline[fi, 3]
        if i>15 and i<=24:
            spherical_jn_offline_ext[fi,i] = spherical_jn_offline[fi, 4]

for fi in range(len(f)):
    
    for i in range(25):
        if i==0:
            spherical_jn_der_offline_ext[fi,i] = spherical_jn_der_offline[fi, 0]
        if i>0 and i<=3:
            spherical_jn_der_offline_ext[fi,i] = spherical_jn_der_offline[fi, 1]
        if i>3 and i<=8:
            spherical_jn_der_offline_ext[fi,i] = spherical_jn_der_offline[fi, 2]
        if i>8 and i<=15:
            spherical_jn_der_offline_ext[fi,i] = spherical_jn_der_offline[fi, 3]
        if i>15 and i<=24:
            spherical_jn_der_offline_ext[fi,i] = spherical_jn_der_offline[fi, 4]

dosya = open("spherical_jn_offline_ext","wb")
pickle.dump(spherical_jn_offline_ext, dosya)
dosya.close()  
  
dosya = open("spherical_jn_der_offline_ext","wb")
pickle.dump(spherical_jn_der_offline_ext, dosya)
dosya.close()       

""" 