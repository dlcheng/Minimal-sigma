import numpy as np
import matplotlib.pyplot as plt

import cosmo
import minimize

def show_legend(ax, location='best'):
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, loc=location)

def plot_variance(org_cosmo, tar_cosmo, zf, zs, s):
	fig = plt.figure(figsize=(5, 5))
	ax = fig.add_subplot(1, 1, 1)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(0.5, 20)
	ax.set_ylim(1e-1, 1e2)
	R = np.linspace(np.log10(0.5), np.log10(60), 15)
	R = np.power(10, R)
	logM_org = np.log10(org_cosmo.mass(R))
	sigma_org = org_cosmo.Sp_logsigma2_logm(logM_org)
	sigma_org = np.power(10., sigma_org)
	logM_tar = np.log10(tar_cosmo.mass(R))
	sigma_tar = tar_cosmo.Sp_logsigma2_logm(logM_tar)
	sigma_tar = np.power(10., sigma_tar)
	ax.plot(R, sigma_org, "ro",  label="SCDM, z=0")
	ax.plot(R, sigma_tar, "ko", label="LCDM, z=0")
	grow_org = org_cosmo.gr(zs)
	R1 = R/s
	logM_org1  = np.log10(org_cosmo.mass(R1))
	sigma_org1 = org_cosmo.Sp_logsigma2_logm(logM_org1)
	sigma_org1 = np.power(10., sigma_org1) * grow_org * grow_org
	ax.plot(R, sigma_org1, "k-", label="Scaled, SCDM")
	ax.set_xlabel("R [Mpc/h]")
	ax.set_ylabel(r"$\sigma^2(R)$")
	show_legend(ax)
	fig.set_tight_layout(True)		
	fig.savefig("Sigma-scaling.pdf")

def plot_fit_countor(org_cosmo, tar_cosmo, R1, R2, zf, zs, s):
	x_dis = 0.01
	y_dis = 0.01
	N_x = 50
	N_y = 50
	x_left  = s - (N_x + 0.5) * x_dis
	x_right = s + (N_x + 0.5) * x_dis
	y_left  = zs - (N_y + 0.5) * y_dis
	y_right = zs + (N_y + 0.5) * y_dis	
	x, y = np.mgrid[slice(x_left, x_right+x_dis, x_dis), slice(y_left, y_right+y_dis, y_dis)]
	rms = np.zeros((x.shape[0]-1, x.shape[1]-1))
	for i in range(rms.shape[0]):
		for j in range(rms.shape[1]):
			rms[i][j] = minimize.rms_function([y[i][j]+y_dis/2., x[i][j]+x_dis/2.], org_cosmo, tar_cosmo, zf, R1, R2)
	cmap = plt.get_cmap('PuRd')
	rms = np.log10(rms)
	fig = plt.figure(figsize=(5, 5))
	ax  = fig.add_subplot(1, 1, 1)
	ax.set_xlabel(r"$s$")
	ax.set_ylabel(r"$z_{*}$")
	cf = ax.contourf(x[:-1, :-1]+x_dis/2., y[:-1,:-1]+y_dis/2., rms, cmap=cmap)
	fig.colorbar(cf, ax=ax)
	fig.set_tight_layout(True)	
	fig.savefig("z-s-countor.pdf")
	plt.show()


	
