#compare results between "slow" FT and FFT

import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
l = float(sys.argv[2])
tau = float(sys.argv[3])

c_r = np.loadtxt(folder + '/spat_corr.txt')
c_t = np.loadtxt(folder + '/time_corr.txt')

c_r_fft = np.loadtxt(folder + '/fft/spat_corr.txt')
c_t_fft = np.loadtxt(folder + '/fft/time_corr.txt')

fig = plt.figure()
plt.plot(c_r[:,0],c_r[:,1]/c_r[0,1],'.-',label='simulation')
plt.plot(c_r[:,0],np.exp(-c_r[:,0]/l),linewidth=0.9,label='theory')
plt.plot(c_r_fft[:,0],c_r_fft[:,1]/c_r_fft[0,1],'.-',linewidth=0.85,label='simulation, FFT')
plt.axhline(y=0, color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
plt.xlim([-0.1,7.1])
plt.ylim([-0.05,1.05])
plt.legend()
plt.savefig(folder + '/spat_corr_compare.pdf')
#plt.show()
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/spat_corr_log_compare.pdf')

fig = plt.figure()
plt.plot(c_t[:,0],c_t[:,1]/c_t[0,1],label='simulation')
plt.plot(c_t[:,0],np.exp(-c_t[:,0]/tau),linewidth=0.9,label='theory')
plt.plot(c_t_fft[:,0],c_t_fft[:,1]/c_t_fft[0,1],linewidth=0.85,label='simulation, FFT')
plt.axhline(y=0, color='black')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
plt.xlim([-0.001,0.05])
plt.legend()
plt.savefig(folder + '/time_corr_compare.pdf')
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/time_corr_log_compare.pdf')

plt.show()