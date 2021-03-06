"""
@author: Ege Erdem

3D PSR Extrapolation - Recording - Run """

""" FIRST RUN pq_Ynm_anm_offline TO UPDATE "anm_offline_whatever" """

import time
from scipy.io.wavfile import read, write
import numpy as np
from scipy.special import sph_harm, spherical_jn
from scipy import signal
import pickle
from convertangle import cart2sph, cart2sph_single, cart2sphr, sph2cart
from sphericals import spherical_bn
import functions3 as ff
#from functions import *

file = open("nemeth_vstacked","rb")
audio = pickle.load(file)
file.close()
    
# loudspekaer vectors 
#  2 pentagon layered 10 loudspeakers, look at the paper for the specific array 
file = open("L","rb")
L = pickle.load(file)
file.close()
L[abs(L) < 1e-8 ] = 0
L_sph = cart2sph(L)
L_ra = 0.155*L # loudspeaker coordinates (as row vectors) as if they are at a 15.5cm virtual spherical surface

L_ra_sph = cart2sph(L_ra)

#unit loudspeaker vectors
L_ra_unit = np.array([L_ra[i]/np.linalg.norm(L_ra[i]) for i in range(len(L_ra))])
L_ra_sph = L_ra_sph[0]

teta_ls = L_ra_sph[:,0]
phi_ls = L_ra_sph[:,1]

   
pi = np.pi
c = 343.0 
ro = 1.0 # ambient density
rc = 0.155
ra = 0.12


def wwrite(signals_, rate, size, w, noverlap, overlap):
    
    _, xrec = signal.istft(signals_[0,:,:], rate, window=w,nperseg=size,noverlap = noverlap)
    
    inversed_signal = np.zeros((10,len(xrec)))
    
    for y in range(10):
        
        _, xrec = signal.istft(signals_[y,:,:], rate, window=w,nperseg=size,noverlap = noverlap)
    
        inversed_signal[y,:] = xrec
    
        
    inversed_signal = inversed_signal/np.max(np.abs(inversed_signal))
    
    
    for i in range(10):
        j = i+1
    
        file = r'C:\Users\INFORMATICS\Desktop\3d-psr-record2\3d-psr-record\psr%d_%d%s%d_r15.wav' %(j,size,w,overlap*100)
#        file = r'/Users/egeerdem/Desktop/3d-psr-record\zpsr%d_%d%s%d.wav' %(j,size,w,overlap*100)
        
        write(file, rate, inversed_signal[i,:])


def main(audio, rate, size, w, overlap, noverlap):
    
#    import pickle
    files = ["spherical_jn_der_offline", "spherical_jn_offline", "sph_harm_offline", "spherical_bn_offline", "sph_harm_offline_der",
         "sph_harm_j", "sph_harm_offline_m","spherical_jn_offline_ext","spherical_jn_der_offline_ext","spherical_bn_offline"]
    
    offlines = []

    filety = open('anm_offline_512hann448_r12',"rb")
    anm_offline = pickle.load(filety)
    filety.close()
    
    for file in files:
        fi = open(file,"rb")
        offlines.append(pickle.load(fi))
        fi.close()
    
    o = {files[0]: offlines[0], files[1]: offlines[1], files[2]: offlines[2],
         files[3]: offlines[3], files[4]: offlines[4], files[5]: offlines[5],
         files[6]: offlines[6], files[7]: offlines[7], files[8]: offlines[8],
         files[9]: offlines[9]}
    
    f, t, _ = signal.stft(audio[0], rate, window=w,nperseg=size, noverlap = noverlap)

    start = time. time()
    signals_ = np.zeros([10, len(f), len(t)])*1j
    
    for ti in range(len(t)):
#    for ti in range(10,20):
        print("ti = ", ti)
    
        for fi in range(len(f)):
#        for fi in range(20,30):
            k = 2*np.pi*f[fi]/c
#            kr = k*rc
            
            if f[fi]<652:
                N=1
                
            elif f[fi]>652 and f[fi]<1301:
                N=2
                
            elif f[fi]>1301 and f[fi]<2607:
                N=3
                
            elif f[fi]>2607:
                N=4
                    
            pressure_ls, I_ls = ff.interpolation(N, k, f, ti, fi, anm_offline, teta_ls, phi_ls, rc, **o)
            l_signals = ff.ls_signals(I_ls, pressure_ls, L_ra_unit)
            signals_[:,fi,ti] = l_signals        
    
    end = time. time()
    sec = end - start
    print(end - start)  
    print(sec/60)
    
    ttx = 'l_signals_%d%s%d_r15' %(size,w,100*overlap)
    
    dosya = open(ttx,"wb")
    pickle.dump(signals_, dosya)
    dosya.close()
    
#    file = open("l_signals_512hann50","rb")
#    signals_ = pickle.load(file)
#    file.close()

    signals_[:,0,:] = 0

    
    print(" window = %s, fft size = %d, overlap = %f yani %d sample " %(w,size,overlap, overlap*size))
    """ *****************************************************************  """
    
    wwrite(signals_, rate, size, w, noverlap, overlap)

if __name__ == "__main__":

    size = 512
    overlap = 7/8

    main(audio, rate=48000, size=512, w ='hann', overlap = 7/8, noverlap = overlap*size)



