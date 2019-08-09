
#from scipy.io.wavfile import read
import numpy as np
import pickle
from scipy.special import sph_harm
from scipy import signal

""" Offline calculation of Ynm and Pq """

teta_q = np.array([69.0, 90.0, 111.0, 90.0, 32.0, 55.0,
                       90.0, 125.0, 148.0, 125.0, 90.0, 55.0, 21.0, 58.0,
                       121.0, 159.0, 69.0, 90.0, 111.0, 90.0, 32.0, 55.0,
                       90.0, 125.0, 148.0, 125.0, 90.0, 55.0, 21.0, 58.0,
                       122.0, 159.0]) / 180.0 * np.pi

phi_q = np.array([0.0, 32.0, 0.0, 328.0, 0.0, 45.0, 69.0, 45.0, 0.0, 315.0,
                     291.0, 315.0, 91.0, 90.0, 90.0, 89.0, 180.0, 212.0, 180.0, 148.0, 180.0,
                     225.0, 249.0, 225.0, 180.0, 135.0, 111.0, 135.0, 269.0, 270.0, 270.0,
                     271.0]) / 180.0 * np.pi

N = 4
counter = 0
Ymn = np.zeros((25,len(teta_q)))*1j

for q in range(len(teta_q)):
    
    for n in range(N+1):
        for m in range(-n,n+1): 
            
            Ymn[counter][q] =  sph_harm(m, n, phi_q[q], teta_q[q]) 
            counter += 1
    counter = 0

file = open("nemeth_vstacked","rb")
audio = pickle.load(file)
file.close()

rate = 48000

#signal.get_window(('kaiser', 4.0), 9)

f, t, _ = signal.stft(audio[0], rate, window='hann',nperseg=512)

audio_stft_32 = np.ones((32,len(f),len(t)))*1j

pq_Ymn = np.ones((25,len(f),len(t)))*1j


for q in range(32):
    
#    f, t, X = signal.stft(audio[q], rate, window=('kaiser', 16.0), nperseg=1024)
    f, t, X = signal.stft(audio[q], rate, window='hann', nperseg=512)

    audio_stft_32[q,:,:] = X

for x in range(len(f)):
    for y in range(len(t)):
        pq_Ymn[:,x,y] = np.matmul(Ymn,audio_stft_32[:,x,y])

dosya = open("pq_Ymn_512hann","wb")
pickle.dump(pq_Ymn, dosya)
dosya.close()
