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

r_list = []
for i in range(c_r_3d_avg.shape[0]):
    r = np.sqrt(c_r_3d_avg[i,0]**2 + c_r_3d_avg[i,1]**2 + c_r_3d_avg[i,2]**2)
    if not(r in r_list):
        r_list.append(r)
r_list.sort()
print(r_list)
c_r_alt_avg = np.zeros(len(r_list))
r_hist_avg = np.zeros(len(r_list))

for i in range(nchunks):
    c_r = np.loadtxt(folder + '/spat_corr_chunk=%d.txt' % i)
    c_r_3d = np.loadtxt(folder + '/spat_corr_3d_chunk=%d.txt' % i)
    c_t = np.loadtxt(folder + '/time_corr_chunk=%d.txt' % i)
    
    c_r_alt = np.zeros(len(r_list))
    r_hist = np.zeros(len(r_list))

    if i>0:
        c_r_avg += c_r
        c_t_avg += c_t
        
    for j in range(c_r_3d.shape[0]):
        r = np.sqrt(c_r_3d[j,0]**2 + c_r_3d[j,1]**2 + c_r_3d[j,2]**2)
        index = r_list.index(r)
        c_r_alt[index] += c_r_3d[j,3]
        r_hist[index] += 1.0
        
    c_r_alt /= r_hist
    c_r_alt_avg += c_r_alt    
    
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
    plt.plot(np.array(r_list),c_r_alt/c_r_alt[0],'.-',label='simulation')
    plt.plot(np.array(r_list),np.exp(-np.array(r_list)/l),linewidth=0.9,label='theory')
    plt.axhline(y=0, color='black')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$C(r)$')
    plt.xlim([-0.1,7.1])
    plt.ylim([-0.05,1.05])
    plt.legend()
    plt.savefig(folder + '/spat_corr_alt_chunk=%d.pdf' % i)
    #plt.show()
    plt.yscale('log')
    plt.ylim([10**(-3),10**0])
    plt.savefig(folder + '/spat_corr_alt_log_chunk=%d.pdf' % i)
    
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
c_r_alt_avg /= nchunks
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
plt.plot(np.array(r_list),c_r_alt_avg/c_r_alt_avg[0],'.-',label='simulation')
plt.plot(np.array(r_list),np.exp(-np.array(r_list)/l),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
plt.xlim([-0.1,7.1])
plt.ylim([-0.05,1.05])
plt.legend()
plt.savefig(folder + '/spat_corr_alt_chunk_avg.pdf')
#plt.show()
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(folder + '/spat_corr_alt_log_chunk_avg.pdf')

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