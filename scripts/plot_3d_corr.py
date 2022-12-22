import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
nx = int(sys.argv[2])
lim_type = sys.argv[3]

lambdas = [1.0,2.0,4.0,8.0]
taus = [0.01, 0.1, 1.0, 10.0]

###################################
#Full tau vs lambda
###################################
#spatial correlation
fig, ax = plt.subplots(4,4,sharey=True, figsize=(8.0,8.0))
for i in range(len(lambdas)):
    for j in range(len(taus)):
        l = lambdas[i]
        tau = taus[j]
        c_r = np.loadtxt('../active-noise-results/hpcc/nx=%d_ny=%d_nz=%d/tau=%s/lambda=%s/Cx_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, nx, nx, tau, l, nx, tau, l))
        dft = np.loadtxt('/home/lfrechette/research/fourier/inverse-exp_3d_onedir_nx=%d_l=%f.txt' % (nx, l))
        ax[i,j].plot(c_r[:,0],c_r[:,1]/c_r[0,1],'.-',label='simulation')
        ax[i,j].plot(c_r[:,0],np.exp(-c_r[:,0]/float(l)),linewidth=0.9,label='theory')
        ax[i,j].plot(dft[:,0], dft[:,1]/dft[0,1], linewidth=0.8, label='DFT', color='red')
        ax[i,j].axhline(y=0, color='black', linestyle='--')
        if lim_type=='corr_len':
            ax[i,j].set_xlim([-0.1,0.1+3*l])
        elif lim_type=='system_size':
            ax[i,j].set_xlim([-0.1,0.1+nx/2-1])
        else:
            print('error: need limit type')
            exit()
        
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
plt.savefig(folder + '/spat_corr_3d_y=0_z=0_nx=%d.pdf' % nx)

#time correlation
fig, ax = plt.subplots(4,4,sharey=True, figsize=(8.0,8.0))
for i in range(len(lambdas)):
    for j in range(len(taus)):
        l = lambdas[i]
        tau = taus[j]
        c_t = np.loadtxt('../active-noise-results/hpcc/nx=%d_ny=%d_nz=%d/tau=%s/lambda=%s/Ct_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, nx, nx, tau, l, nx, tau, l))
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
plt.savefig(folder + '/time_corr_3d_y=0_z=0_nx=%d.pdf' % nx)

#plt.show()