import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
nx = int(sys.argv[2])

lambdas = [1.0,2.0,4.0,8.0]

fig, ax = plt.subplots(4,1,sharey=True, figsize=(3.0,8.0))
for i in range(len(lambdas)):
    l = lambdas[i]
    c_r = np.loadtxt('../active-noise-results/hpcc/Cr_total_nx=%d_tau=0.010000_lambda=%f.txt' % (nx, l))
    dft = np.loadtxt('/home/lfrechette/research/fourier/inverse-exp_3d_onedir_nx=%d_l=%f.txt' % (nx, l))
    ax[i].plot(c_r[:,0],c_r[:,1],'.-',label='simulation')
    ax[i].plot(c_r[:,0],c_r[:,2],linewidth=0.9,label='theory')
    ax[i].plot(dft[:,0],dft[:,1]/dft[0,1],linewidth=0.8,color='red',label='DFT')
    ax[i].axhline(y=0, color='black', linestyle='--')
    ax[i].set_xlabel(r'$r$')
    ax[i].set_ylabel(r'$C(r)$')
    ax[i].text(7, 0.8, r'$\lambda=%.00f$' % l)
    ax[i].set_xlim([-0.1, nx/2.0+0.1])
#plt.xlim([-0.1,7.1])
#plt.ylim([-0.05,1.05])
    if i==0:
        ax[i].set_title(r'$\tau=0.01$')
    if i==3:
        ax[i].legend()
plt.savefig(folder + '/spat_corr_vs_lambda_nx=%d_tau=0.01.pdf' % nx)

taus = [0.01, 0.1, 1.0, 10.0]
fig, ax = plt.subplots(4,1,sharey=True, figsize=(3.0,8.0))
for i in range(len(taus)):
    tau = taus[i]
    c_t = np.loadtxt('../active-noise-results/hpcc/Ct_total_nx=%d_tau=%f_lambda=1.000000.txt' % (nx, tau))
    ax[i].plot(c_t[:,0],c_t[:,1],'.-',label='simulation')
    ax[i].plot(c_t[:,0],c_t[:,2],linewidth=0.9,label='theory')
    ax[i].axhline(y=0, color='black', linestyle='--')
    ax[i].set_xlabel(r'$t$')
    ax[i].set_ylabel(r'$C(t)$')
    ax[i].text(c_t[-1,0]*0.8, 0.8, r'$\tau=%.00f$' % tau)
#plt.xlim([-0.1,7.1])
#plt.ylim([-0.05,1.05])
    if i==0:
        ax[i].set_title(r'$\lambda=1$')
    if i==3:
        ax[i].legend()
plt.savefig(folder + '/time_corr_vs_tau_nx=%d_lambda=1.pdf' % nx)
#plt.show()
'''
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
'''
plt.show()