"""
@author: Ege Erdem

3D PSR Extrapolation - Recording - Run """

""" FIRST RUN pq_Ynm_anm_offline TO UPDATE "anm_offline_whatever" """

file = open('anm_offlineCOORorderconj_512hann448',"rb")
anm_offline = pickle.load(file)
file.close()

import time
from scipy.io.wavfile import read, write
import numpy as np
from scipy.special import sph_harm, spherical_jn
from scipy import signal
import pickle
from convertangle import cart2sph, cart2sph_single, cart2sphr, sph2cart
from sphericals import spherical_bn
#from functions import *

#file = open("nemeth_vstacked","rb")
#audio = pickle.load(file)
#file.close()
rate = 48000

""" *****************************************************************  """
size = 512
w = 'hann'

#overlap = size-1
overlap = 7/8
noverlap = overlap*size
#tt = 'anm_offline_%d%s%d' %(size,w,overlap)


print(" window = %s, fft size = %d, overlap = %f yani %d sample " %(w,size,overlap, overlap*size))
""" *****************************************************************  """


f, t, _ = signal.stft(audio[0], rate, window=w,nperseg=size, noverlap = noverlap)
#f, t, _ = signal.stft(audio[0], rate, window=w,nperseg=size, noverlap = overlap)

start = time. time()

n0 =  np.zeros((10,25))
signals_ = np.zeros([10, len(f), len(t)])*1j
#for fi in range(len(f)):


for ti in range(len(t)):

    print("ti = ", ti)

    for fi in range(len(f)):
        
        k = 2*np.pi*f[fi]/c
        kr = k*rc
        
        if f[fi]<652:
            N=2
            
        elif f[fi]>652 and f[fi]<1301:
            N=2
            
        elif f[fi]>1301 and f[fi]<2607:
            N=3
            
        elif f[fi]>2607:
            N=4
                
        pressure_ls, I_ls = interpolation(N, k, ti, fi, anm_offline)
        l_signals = ls_signals(I_ls, pressure_ls)
        signals_[:,fi,ti] = l_signals        

end = time. time()
print(end - start)   


ttx = 'l_signalsCORRorderconj_%d%s%d' %(size,w,100*overlap)

dosya = open(ttx,"wb")
pickle.dump(signals_, dosya)
dosya.close()

#file = open("l_signals_512hann50","rb")
#signals_ = pickle.load(file)
#file.close()
#overlap = 50

signals_[:,0,:] = 0

_, xrec = signal.istft(signals_[0,:,:], rate, window=w,nperseg=size,noverlap = noverlap)

inversed_signal = np.zeros((10,len(xrec)))

for y in range(10):
    
    _, xrec = signal.istft(signals_[y,:,:], rate, window=w,nperseg=size,noverlap = noverlap)

    inversed_signal[y,:] = xrec

    
inversed_signal = inversed_signal/np.max(np.abs(inversed_signal))


for i in range(10):
    j = i+1    
    file = r'D:\MASTER\_Tez_SPARG\CODES\3d-psr-record\psr%dCORRoc_%d%s%d.wav' %(j,size,w,overlap*100)
    write(file, rate, inversed_signal[i,:])

