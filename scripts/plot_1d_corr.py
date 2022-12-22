import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
nx = int(sys.argv[2])

lambdas = [1.0,2.0,4.0,8.0]

fig, ax = plt.subplots(4,1,sharey=True, figsize=(3.0,8.0))
for i in range(len(lambdas)):
    l = lambdas[i]
    c_r = np.loadtxt('../active-noise-results/hpcc/1d/nx=%d/tau=0.01/lambda=%s/spat_corr.txt' % (nx, l))
    ax[i].plot(c_r[:,0],c_r[:,1]/c_r[0,1],'.-',label='simulation')
    ax[i].plot(c_r[:,0],np.exp(-c_r[:,0]/float(l)),linewidth=0.9,label='theory')
    ax[i].axhline(y=0, color='black', linestyle='--')
    ax[i].set_xlabel(r'$r$')
    ax[i].set_ylabel(r'$C(r)$')
    ax[i].text(c_r[-1,0]*0.7, 0.8, r'$\lambda=%.00f$' % l)
#plt.xlim([-0.1,7.1])
#plt.ylim([-0.05,1.05])
    if i==0:
        ax[i].set_title(r'$\tau=0.01$')
    if i==3:
        ax[i].legend()
plt.savefig(folder + '/spat_corr_1d_nx=%d_vs_lambda_tau=0.01.pdf' % nx)

taus = [0.01, 0.1, 1.0, 10.0]
fig, ax = plt.subplots(4,1,sharey=True, figsize=(3.0,8.0))
for i in range(len(taus)):
    tau = taus[i]
    c_t = np.loadtxt('../active-noise-results/hpcc/1d/nx=%d/tau=%s/lambda=1.0/time_corr.txt' % (nx, tau))
    ax[i].plot(c_t[:,0],c_t[:,1]/c_t[0,1],'.-',label='simulation')
    ax[i].plot(c_t[:,0],np.exp(-c_t[:,0]/float(tau)),linewidth=0.9,label='theory')
    ax[i].axhline(y=0, color='black', linestyle='--')
    ax[i].set_xlabel(r'$t$')
    ax[i].set_ylabel(r'$C(t)$')
    ax[i].text(c_t[-1,0]*0.7, 0.8, r'$\tau=%.02f$' % tau)
#plt.xlim([-0.1,7.1])
#plt.ylim([-0.05,1.05])
    if i==0:
        ax[i].set_title(r'$\lambda=1$')
    if i==3:
        ax[i].legend()
plt.savefig(folder + '/time_corr_1d_nx=%d_vs_tau_lambda=1.pdf' % nx)
#plt.show()

###################################
#Full tau vs lambda
###################################
#spatial correlation
fig, ax = plt.subplots(4,4,sharey=True, figsize=(8.0,8.0))
for i in range(len(lambdas)):
    for j in range(len(taus)):
        l = lambdas[i]
        tau = taus[j]
        c_r = np.loadtxt('../active-noise-results/hpcc/1d/nx=%d/tau=%s/lambda=%s/spat_corr.txt' % (nx, tau, l))
        ax[i,j].plot(c_r[:,0],c_r[:,1]/c_r[0,1],'.-',label='simulation')
        ax[i,j].plot(c_r[:,0],np.exp(-c_r[:,0]/float(l)),linewidth=0.9,label='theory')
        ax[i,j].axhline(y=0, color='black', linestyle='--')
        #ax[i,j].set_xlim([-0.1,0.1+nx/2])
        ax[i,j].set_xlim([-0.1,0.1+3*l])
        ax[i,j].set_xlabel(r'$r$')

        if j==len(taus)-1:
            ax2 = ax[i,j].twinx()
            ax2.set_yticks([],[])
            ax2.set_ylabel('$\lambda=%.00f$' % l)
#plt.xlim([-0.1,7.1])
#plt.ylim([-0.05,1.05])
        if j==0:
            ax[i,j].set_ylabel(r'$C(r)$')
        if i==0:
            ax[i,j].set_title(r'$\tau=%s$' % tau)
        if i==0 and j==0:
            ax[i,j].legend(fontsize=10)
plt.savefig(folder + '/spat_corr_1d_nx=%d.pdf' % nx)

#time correlation
fig, ax = plt.subplots(4,4,sharey=True, figsize=(8.0,8.0))
for i in range(len(lambdas)):
    for j in range(len(taus)):
        l = lambdas[i]
        tau = taus[j]
        c_t = np.loadtxt('../active-noise-results/hpcc/1d/nx=%d/tau=%s/lambda=%s/time_corr.txt' % (nx, tau, l))
        ax[i,j].plot(c_t[:,0],c_t[:,1]/c_t[0,1],'.-',label='simulation')
        ax[i,j].plot(c_t[:,0],np.exp(-c_t[:,0]/float(tau)),linewidth=0.9,label='theory')
        ax[i,j].axhline(y=0, color='black', linestyle='--')
        #ax[i,j].set_xlim([-0.1,0.1+nx/2])
        ax[i,j].set_xlabel(r'$t$')

        if j==len(taus)-1:
            ax2 = ax[i,j].twinx()
            ax2.set_yticks([],[])
            ax2.set_ylabel('$\lambda=%.00f$' % l)
#plt.xlim([-0.1,7.1])
#plt.ylim([-0.05,1.05])
        if j==0:
            ax[i,j].set_ylabel(r'$C(t)$')
        if i==0:
            ax[i,j].set_title(r'$\tau=%s$' % tau)
        if i==0 and j==0:
            ax[i,j].legend(fontsize=10)
plt.savefig(folder + '/time_corr_1d_nx=%d.pdf' % nx)

#plt.show()