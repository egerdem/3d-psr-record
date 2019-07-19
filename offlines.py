# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 23:37:55 2019

@author: ABRA A5
"""
""" Nemeth recordings loaded into 32*1488000 vstacked numpy array """

import pickle

file = open("spherical_jn_der_offline","rb")
spherical_jn_der_offline = pickle.load(file)
file.close()

file = open("spherical_jn_offline","rb")
spherical_jn_offline = pickle.load(file)
file.close()

file = open("sph_harm_offline","rb")
sph_harm_offline = pickle.load(file)
file.close()

file = open("spherical_bn_offline","rb")
spherical_bn_offline = pickle.load(file)
file.close()

file = open("sph_harm_offline_der","rb")
sph_harm_offline_der = pickle.load(file)
file.close()

file = open("sph_harm_j","rb")
sph_harm_j = pickle.load(file)
file.close()

file = open("sph_harm_offline_m","rb")
sph_harm_offline_m = pickle.load(file)
file.close()

file = open("spherical_jn_offline_ext","rb")
spherical_jn_offline_ext = pickle.load(file)
file.close()

file = open("spherical_jn_der_offline_ext","rb")
spherical_jn_der_offline_ext = pickle.load(file)
file.close()

file = open("spherical_bn_offline","rb")
spherical_bn_offline = pickle.load(file)
file.close()

file = open("nemeth_vstacked","rb")
audio = pickle.load(file)
file.close()
print("kodu")
#file = open("audio_stft","rb")
#audio_stft = pickle.load(file)
#file.close()  
def off():
    print("off fonks çalıştı")