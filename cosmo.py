import constants

import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate
from scipy import interpolate

# for interpolate.InterpolatedUnivariateSpline, the input x array must be increasing.

class cosmo:
	Omega_m = 0.0
	Omega_v = 0.0
	Omega_b = 0.0
	
	# Varibales for power spectrum
	Aps = 0.0                     # the normalization parameter of the PS at z=0
	Agr = 0.0                     # the normalization parameter of linear growth
	ns = 0.0                      #
	Sigma8 = 0.0                  #
	Ps_flag = 0.0                 # flags for different types of power spectrum
	                              # flag == 1,   scale free
	                              # flag == 2,   BBKS (SCDM, default \Gamma = 0.5)
	                              # flag == 3,   E&H  (power with BAO)
	                              # flag == 4,   Read transfer function from file
	                              
	# Variables for spline functions
	Tb_M = 0.0                    # list of the mass of interest
	Tb_sigma2 = 0.0               # list of the sigma^2 at z=0 as function of M
	Sp_logsigma2_logm = 0.0       # the 1d spline function of log10(sigma2) with log10(M) at z = 0
	Num_M = 0.0                   # number of list Tb_M
	M_min = 0.0                   # minimum value of Tb_M in unit of 10^10 M_sun/h
	M_max = 0.0                   # maximum value of Tb_M
	Tb_tk = 0.0                   # the list of transfer function in case of reading from file.
	Sp_tk = 0.0                   # the spline funcion of tk from file
	
#### initilization
	def __init__(self, Omega_m, Omega_v, Omega_b, ns, Sigma8, Ps_flag, M_min, M_max, Num_M):
		self.set_spline_parameters(M_min, M_max, Num_M)
		self.set_cosmological_parameters(Omega_m, Omega_v, Omega_b)
		self.set_power_parameters(ns, Sigma8, Ps_flag)
		self.init_M()
		self.init_growth_factor()
		self.init_power_spectrum()
		self.init_logsigma2_slpine()	

#### functions used to initialize the objects ###########
	def set_cosmological_parameters(self, Omega_m, Omega_v, Omega_b):
		self.Omega_m = Omega_m
		self.Omega_v = Omega_v
		self.Omega_b = Omega_b

	def set_power_parameters(self, ns, Sigma8, Ps_flag):
		self.ns = ns
		self.Sigma8 = Sigma8
		self.Ps_flag = Ps_flag		

	def set_spline_parameters(self, M_min, M_max, Num_M):
		self.M_min = M_min
		self.M_max = M_max
		self.Num_M = Num_M	

	def init_M(self):
		dis_M = np.log10(self.M_max / self.M_min) / (self.Num_M - 1)
		self.Tb_M = np.zeros((self.Num_M,))	
		for i in range(self.Num_M):
			self.Tb_M[i] = np.power(10, np.log10(self.M_min) + i * dis_M)

	def init_growth_factor(self):
		self.Agr = 1.0 / self.unnorm_growth(0.0)	

	def init_power_spectrum(self):
		if self.Ps_flag == 4:
			self.init_tk_from_file()
		self.norm_ps()

	def init_logsigma2_slpine(self):
		self.Tb_sigma2 = np.zeros((self.Num_M,))               # the first step is to prepare the data
		for i in range(self.Num_M):
			self.Tb_sigma2[i] = self.calculate_sigma2(self.Tb_M[i])
		self.Sp_logsigma2_logm = interpolate.InterpolatedUnivariateSpline(np.log10(self.Tb_M), np.log10(self.Tb_sigma2))		

#### functions used to calculate the growth factor ###
	def g_factor(self, z):
		x = self.Omega_m * np.power(1.0+z, 3)
		x += (1.0 - self.Omega_m - self.Omega_v) * np.power(1.0+z, 2) + self.Omega_v
		return x  

	def unnorm_growth(self, z):
		g = self.g_factor(z)
		y1 = self.Omega_m * np.power(1.0+z, 3)/g
		y2 = self.Omega_v /g
		x = 1.0 / (1.0 + z) * 2.5 * y1
		x /= np.power(y1, 4.0/7.0) - y2 + (1.0 + 0.5 * y1) * (1.0 + y2/70.0)
		return x

	def gr(self, z):
		return self.Agr * self.unnorm_growth(z)

#### functions used to normalize the power spectrum and calculate sigma2--M ###
	def calculate_sigma2(self, M):
		R = self.radius(M)
		return self.Aps * self.unnorm_sigma2(R)

	def win_th(self, y):
		x = 1
		if y <= 1e-8:
			x = 1.0 - 0.1 * y * y
		else:
			x = 3.0 / np.power(y, 3) * (np.sin(y) - y * np.cos(y))
		return x

	def sigma2_int_kernel(self, k, R):
		x = 0.5 / np.pi / np.pi
		x *= np.power(k, self.ns+2.0)
		x *= np.power(self.prim_tk(k), 2)
		x *= np.power(self.win_th(k*R), 2)
		return x

	def unnorm_sigma2(self, R):
		x, err = integrate.quad(self.sigma2_int_kernel, 0, np.inf, (R,), epsabs=0.0, epsrel=1e-6)
		return x

	def norm_ps(self):
		if self.Ps_flag == 1:
			self.Aps = 1.0
		else:	
			self.Aps = self.Sigma8 * self.Sigma8 / self.unnorm_sigma2(8.0)

#### functions used to calculate the transfer function ###
	def prim_tk(self, k):
		x = 1
		if self.Ps_flag == 1:
			x = 1
		if self.Ps_flag == 2:
			x = self.tk_bbks(k)
		if self.Ps_flag == 3:
			pass                          # return E&H transfer function
		if self.Ps_flag == 4:
			x = self.tk_file(k)           # return tk from file
		return x

	def tk_bbks(self, k):
		x = 1.0
		T_sig = 0.5                       # shape parameter for SCDM
		q = k / T_sig
		a0 = 2.34
		a1 = 3.89
		a2 = 16.19
		a3 = 5.46
		a4 = 6.71
		if q >= 1e-8:                
			b1 = 1.0 
			b1 += a1 * q
			b1 += a2 * a2 * q * q
			b1 += a3 * a3 * a3 * q * q * q
			b1 += a4 * a4 * a4 * a4 * q * q * q * q
			b1 = 1.0 / np.power(b1, 0.25)
			b2 = np.log(1.0 + a0 * q) / a0 / q
			x = b1 * b2
		return x

	def init_tk_from_file(self):
		file_name = './Input/transf0.dat'
		self.Tb_tk = np.loadtxt(file_name)
		fb = self.Omega_b / self.Omega_m
		fb = 0.166                       # for the WMAP shape of Jing's simulations
		self.Tb_tk[:,1] = (1-fb) * self.Tb_tk[:,1] + fb * self.Tb_tk[:,2]
		self.Sp_tk = interpolate.InterpolatedUnivariateSpline(np.log10(self.Tb_tk[:,0]), self.Tb_tk[:,1]/self.Tb_tk[0,1]) 

	def tk_file(self, k):    
		if k >= self.Tb_tk[0,0] and k<= self.Tb_tk[-1,0]:
			return self.Sp_tk(np.log10(k))
		if k < self.Tb_tk[0,0]:
			return 1.0
		if k > self.Tb_tk[-1,0]:
			return 0.0

#### functions used to calculate the linear power ########
	def lin_ps(self, k, z):
		x = self.Aps * np.power(k, self.ns)
		x *= np.power(self.prim_tk(k) , 2)
		x *= np.power(self.gr(z), 2)
		return x

	def lin_dt(self, k, z):
		x = self.lin_ps(k, z)
		x *= np.power(k, 3)/2/np.pi/np.pi
		return x

#### functions used to calculate the R-M relation
	def radius(self, M):
		x = M * 3.0 / 4.0 / np.pi / (constants.rho_crit * self.Omega_m)
		x = np.power(x, 1.0/3.0)
		return x

	def mass(self, R):
		x = 4.0 * np.pi / 3.0 * np.power(R, 3.0) * (constants.rho_crit * self.Omega_m)
		return x	
