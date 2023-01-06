import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
nx = int(sys.argv[2])
lim_type = sys.argv[3]

#lambdas = [1.0,2.0,4.0,8.0]
lambdas = [8.0,16.0,32.0,64.0]

###################################
#Full tau vs lambda
###################################
#spatial correlation
fig, ax = plt.subplots(len(lambdas),1,sharey=True, figsize=(3.0,8.0))
for i in range(len(lambdas)):
    l = lambdas[i]
    dft = np.loadtxt('/home/lfrechette/research/fourier/inverse-exp_3d_onedir_nx=%d_l=%f.txt' % (nx, l))
    ax[i].plot(dft[:,0],np.exp(-dft[:,0]/float(l)),linewidth=0.9,label='FT', color='orange')
    ax[i].plot(dft[:,0], dft[:,1]/dft[0,1], linewidth=0.8, label='DFT', color='red')
    ax[i].axhline(y=0, color='black', linestyle='--')
    if lim_type=='corr_len':
        ax[i].set_xlim([-0.1,0.1+3*l])
    elif lim_type=='system_size':
        ax[i].set_xlim([-0.1,0.1+nx/2-1])
    else:
        print('error: need limit type')
        exit()
    
    ax[i].set_xlabel(r'$r$')

    ax2 = ax[i].twinx()
    ax2.set_yticks([],[])
    ax2.set_ylabel('$\lambda=%.00f$' % l)
    ax[i].set_ylabel(r'$C(r)$')
    if i==0:
        ax[i].legend(fontsize=10)
plt.savefig(folder + '/spat_corr_3d_dft_vs_ft_%s_y=0_z=0_nx=%d.pdf' % (lim_type, nx))

#plt.show()