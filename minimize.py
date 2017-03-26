import cosmo
import constants

import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize
from scipy import integrate


def rms_function(x, org_cosmo, tar_cosmo, zf, R1, R2):
	z = x[0]
	s = x[1]
	y, err = integrate.quad(intg_function_kernel, R1, R2, (z,s, org_cosmo, tar_cosmo, zf), epsabs=0.0, epsrel=1e-6)
	return y/np.log(R2/R1)

def intg_function_kernel(R, z, s, org_cosmo, tar_cosmo, zf):
	grow_org = org_cosmo.gr(z)                # growth factor of original cosmology
	grow_tar = tar_cosmo.gr(zf)               # growth factor of target cosmology
	logM_org = np.log10(org_cosmo.mass(R/s)) 
	logM_tar = np.log10(tar_cosmo.mass(R))
	sigma_org = org_cosmo.Sp_logsigma2_logm(logM_org)
	sigma_org = np.power(10., sigma_org)
	sigma_org = np.sqrt(sigma_org) * grow_org
	sigma_tar = tar_cosmo.Sp_logsigma2_logm(logM_tar)
	sigma_tar = np.power(10., sigma_tar)
	sigma_tar = np.sqrt(sigma_tar) * grow_tar
	y = 1.0 - sigma_org / sigma_tar
	return np.power(y, 2)/R

def mini_rms(org_cosmo, tar_cosmo, zf, R1, R2):
	x0 = [0.5, 2.]                           # initial guess
	res = optimize.minimize(rms_function, x0, (org_cosmo, tar_cosmo, zf, R1, R2), method='nelder-mead', options={'disp':False, 'xtol':1e-5, 'ftol':1e-8})
	return res.x[0], res.x[1]