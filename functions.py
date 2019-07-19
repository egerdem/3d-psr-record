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

if __name__ == "__main__":
    
    pi = np.pi
    """
    # change this later
    rate = 48000
    """
    c = 343.0 
    ro = 1.0 # ambient density
    rc = 0.155
    ra = 0.12
    
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
        
    teta_q = np.array([69.0, 90.0, 111.0, 90.0, 32.0, 55.0,
                           90.0, 125.0, 148.0, 125.0, 90.0, 55.0, 21.0, 58.0,
                           121.0, 159.0, 69.0, 90.0, 111.0, 90.0, 32.0, 55.0,
                           90.0, 125.0, 148.0, 125.0, 90.0, 55.0, 21.0, 58.0,
                           122.0, 159.0]) / 180.0 * np.pi
    
    phi_q = np.array([0.0, 32.0, 0.0, 328.0, 0.0, 45.0, 69.0, 45.0, 0.0, 315.0,
                         291.0, 315.0, 91.0, 90.0, 90.0, 89.0, 180.0, 212.0, 180.0, 148.0, 180.0,
                         225.0, 249.0, 225.0, 180.0, 135.0, 111.0, 135.0, 269.0, 270.0, 270.0,
                         271.0]) / 180.0 * np.pi
    print("functions loaded from functions.py")

#def psr():
    
    
def grad(teta_ls, phi_ls, k, N, ti, fi,anm_offline,**o):
    teta = teta_ls
    phi = phi_ls
    """ pressure grad in spherical coordinates, to be used for obtaining velocity and then intensity """
    r_unit = np.array([np.sin(teta)*np.cos(phi), np.sin(teta)*np.sin(phi) , np.cos(teta)]).T
    p_unit = np.array([-np.sin(phi)*np.sin(teta), np.sin(teta)*np.cos(phi) , np.array([0,0,0,0,0,0,0,0,0,0])]).T
    t_unit = np.array([np.cos(teta)*np.cos(phi), np.cos(teta)*np.sin(phi) , -np.sin(teta)]).T
#    kr = k*rc
    
    # pr_der
    
    matrix1 = o['sph_harm_offline'].T
    matrix2 = anm_offline[:, fi,ti][:, np.newaxis]*4*pi*k*o['spherical_jn_der_offline_ext'][fi,:][:, np.newaxis]
    pr_der = np.matmul(matrix1, matrix2).T

    # pteta_der
    
    matrix3 = o['sph_harm_offline_der'].T
    matrix4 = anm_offline[:, fi,ti][:, np.newaxis]*4*pi*o['spherical_jn_offline_ext'][fi,:][:, np.newaxis]
    pteta_der = np.matmul(matrix3, matrix4).T

    
    # pphi_der
    
    matrix5 = o['sph_harm_offline_m'].T
    matrix6 = anm_offline[:, fi,ti][:, np.newaxis]*4*pi*o['spherical_jn_offline_ext'][fi,:][:, np.newaxis]*1j
    pphi_der = np.matmul(matrix5, matrix6).T

    grad_sum = r_unit * pr_der.T + t_unit * pteta_der.T*(1/rc) + p_unit * pphi_der.T*(1/(rc*np.sin(teta)[:, np.newaxis]))

    return(grad_sum)
    
    
def calc_u(pressure_grad, ro, k,f):
    w = 2*np.pi*f
    u = -pressure_grad/(1j*w*ro)
    
    return(u)


def pressure_field(fi, ti, N,anm_offline,**o): 
    """ extrapolated pressure at distance r due to a monochromatic plane wave coming from *** (teta,azimuth) """
    
    matrixa1 = o['sph_harm_offline'].T
    matrixb1 = anm_offline[:, fi,ti][:, np.newaxis]*4*pi*o['spherical_jn_offline_ext'][fi,:][:, np.newaxis]
    
    pr1 = np.matmul(matrixa1, matrixb1).T

    return(np.squeeze(pr1))
    
def ia(p, u):
    """ local active intensity with known pressure and velocity """
    I_a = 0.5*(p.T[:, np.newaxis]*np.conj(u)).real
    In = np.array([I_a[i]/np.linalg.norm(I_a[i]) for i in range(len(I_a))])
    
    return(I_a, In)
    
def interpolation(N, k, f, ti, fi,anm_offline,**o):
      
    p_grad = grad(teta_ls, phi_ls, k, N, ti, fi,anm_offline, **o) # grad of the pressure
    u = calc_u(p_grad, ro, k, f[fi])  # local velocity calculation using euler momentum equation with pressure grad
    p = pressure_field(fi, ti, N,anm_offline,**o)  # extrapolated pressure at r 
    
    I_act, In = ia(p, u)  # local active intensity   

    return(p,I_act)

    
    
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
    
