import numpy as np
import pylab as plt
import sys

folder = sys.argv[1]
l = float(sys.argv[2])
tau = float(sys.argv[3])
nchunks = int(sys.argv[4])

c_r_avg = np.loadtxt(folder + '/spat_corr_chunk=0.txt')
c_r_3d_avg = np.loadtxt(folder + '/spat_corr_3d_chunk=0.txt')
c_t_avg = np.loadtxt(folder + '/time_corr_chunk=0.txt')

for i in range(nchunks):
    c_r = np.loadtxt(folder + '/spat_corr_chunk=%d.txt' % i)
    c_r_3d = np.loadtxt(folder + '/spat_corr_3d_chunk=%d.txt')
    c_t = np.loadtxt(folder + '/time_corr_chunk=%d.txt' % i)

    if i>0:
        c_r_avg += c_r
        c_t_avg += c_t
        c_r_3d_avg += c_r_3d
    
    fig = plt.figure()
    plt.plot(c_r[:,0],c_r[:,1]/c_r[0,1],'.-',label='simulation')
    plt.plot(c_r[:,0],np.exp(-c_r[:,0]/l),linewidth=0.9,label='theory')
    plt.axhline(y=0, color='black')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$C(r)$')
    plt.xlim([-0.1,7.1])
    plt.ylim([-0.05,1.05])
    plt.legend()
    plt.savefig(folder + '/spat_corr_chunk=%d.pdf' % i)
    #plt.show()
    plt.yscale('log')
    plt.ylim([10**(-3),10**0])
    plt.savefig(folder + '/spat_corr_log_chunk=%d.pdf' % i)
    
    fig = plt.figure()
    plt.plot(c_t[:,0],c_t[:,1]/c_t[0,1],label='simulation')
    plt.plot(c_t[:,0],np.exp(-c_t[:,0]/tau),linewidth=0.9,label='theory')
    plt.axhline(y=0, color='black')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$C(t)$')
    plt.xlim([-0.001,0.05])
    plt.legend()
    plt.savefig(folder + '/time_corr_chunk=%d.pdf' % i)
    plt.yscale('log')
    plt.ylim([10**(-3),10**0])
    plt.savefig(folder + '/time_corr_log_chunk=%d.pdf' % i)
    
plt.close()
c_r_avg /= nchunks
c_t_avg /= nchunks

fig = plt.figure()
plt.plot(c_r_avg[:,0],c_r_avg[:,1]/c_r_avg[0,1],'.-',label='simulation')
plt.plot(c_r_avg[:,0],np.exp(-c_r_avg[:,0]/l),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
plt.xlim([-0.1,7.1])
plt.ylim([-0.05,1.05])
plt.legend()
plt.savefig(folder + '/spat_corr_chunk_avg.pdf')
#plt.show()
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/spat_corr_log_chunk_avg.pdf')

fig = plt.figure()
plt.plot(c_t_avg[:,0],c_t_avg[:,1]/c_t_avg[0,1],label='simulation')
plt.plot(c_t_avg[:,0],np.exp(-c_t_avg[:,0]/tau),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
plt.xlim([-0.001,0.05])
plt.legend()
plt.savefig(folder + '/time_corr_chunk_avg.pdf')
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/time_corr_log_chunk_avg.pdf')

#plt.show()