import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
l = float(sys.argv[2])
tau = float(sys.argv[3])

c_r = np.loadtxt(folder + '/spat_corr.txt')
c_t = np.loadtxt(folder + '/time_corr.txt')

fig = plt.figure()
plt.plot(c_r[:,0],c_r[:,1]/c_r[0,1],'.-',label='simulation')
plt.plot(c_r[:,0],np.exp(-c_r[:,0]/l),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
plt.xlim([-0.1,7.1])
plt.ylim([-0.05,1.05])
plt.legend()
plt.savefig(folder + '/spat_corr.pdf')
#plt.show()
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/spat_corr_log.pdf')

fig = plt.figure()
plt.plot(c_t[:,0],c_t[:,1]/c_t[0,1],label='simulation')
plt.plot(c_t[:,0],np.exp(-c_t[:,0]/tau),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
plt.xlim([-0.001,0.05])
plt.legend()
plt.savefig(folder + '/time_corr.pdf')
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/time_corr_log.pdf')

plt.show()