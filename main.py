# Minimize Eq.(4) of Angulo and White (2009) to obtain the best length scaling relation and redshift 
# mapping between Simulation 8411 and 8600.

import constants
import cosmo
import plot
import minimize
import plot

import numpy as np
import matplotlib.pyplot as plt

reload(constants)
reload(cosmo)
reload(plot)
reload(minimize)
reload(plot)


#set parameters for 8600, 8610
cos_8600 = cosmo.cosmo(1.0,   0.0,   0.166, 1.0,   0.83, 2, 1e5, 1e19, 1000)  #SCDM
cos_8610 = cosmo.cosmo(1.0,   0.0,   0.166, 0.968, 0.83, 4, 1e5, 1e19, 1000)  #WMAP ps
cos_8411 = cosmo.cosmo(0.268, 0.732, 0.166, 0.968, 0.83, 4, 1e5, 1e19, 1000)  #LCDM WMAP

#
R1 = 1.
R2 = 10.
zs, s = minimize.mini_rms(cos_8600, cos_8411, 0., R1, R2)                     #zs is the starting redshift of 8600
                                                                              #s is the length scale ratio

plot.plot_variance(cos_8600, cos_8411, 0, zs, s)
plot.plot_fit_countor(cos_8600, cos_8411, R1, R2, 0, zs, s)

print "From 8600 to 8411 at z=0, we need"
print "z* = %.6e, s=%.6e"%(zs, s)
