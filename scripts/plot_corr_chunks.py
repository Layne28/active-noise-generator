import numpy as np
import pylab as plt
import sys

base_folder = sys.argv[1]
out_folder = sys.argv[2]
nx = int(sys.argv[3])
l = sys.argv[4]
tau = sys.argv[5]
nchunks = int(sys.argv[6])

folder = base_folder + '/nx=%d_ny=%d_nz=%d/tau=%s/lambda=%s/' % (nx, nx, nx, tau, l)

l = float(l)
tau = float(tau)

c_r_avg = np.loadtxt(folder + '/spat_corr_chunk=0.txt')
c_r_3d_avg = np.loadtxt(folder + '/spat_corr_3d_chunk=0.txt')
c_t_avg = np.loadtxt(folder + '/time_corr_chunk=0.txt')

r_list = []
for i in range(c_r_3d_avg.shape[0]):
    r = np.sqrt(c_r_3d_avg[i,0]**2 + c_r_3d_avg[i,1]**2 + c_r_3d_avg[i,2]**2)
    if not(r in r_list):
        r_list.append(r)
r_list.sort()
#print(r_list)
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
    plt.savefig(out_folder + '/chunks/spat_corr_nx=%d_tau=%f_lambda=%f_chunk=%d.pdf' % (nx, tau, l, i))
    #plt.show()
    plt.yscale('log')
    plt.ylim([10**(-3),10**0])
    plt.savefig(out_folder + '/chunks/spat_corr_log_nx=%d_tau=%f_lambda=%f_chunk=%d.pdf' % (nx, tau, l, i))
    
    fig = plt.figure()
    plt.plot(np.array(r_list),c_r_alt/c_r_alt[0],'.-',label='simulation')
    plt.plot(np.array(r_list),np.exp(-np.array(r_list)/l),linewidth=0.9,label='theory')
    plt.axhline(y=0, color='black')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$C(r)$')
    #plt.xlim([-0.1,7.1])
    plt.ylim([-0.05,1.05])
    plt.legend()
    plt.savefig(out_folder + '/chunks/spat_corr_alt_nx=%d_tau=%f_lambda=%f_chunk=%d.pdf' % (nx, tau, l, i))
    #plt.show()
    plt.yscale('log')
    plt.ylim([10**(-3),10**0])
    plt.savefig(out_folder + '/chunks/spat_corr_alt_log_nx=%d_tau=%f_lambda=%f_chunk=%d.pdf' % (nx, tau, l, i))
    
    fig = plt.figure()
    plt.plot(c_t[:,0],c_t[:,1]/c_t[0,1],label='simulation')
    plt.plot(c_t[:,0],np.exp(-c_t[:,0]/tau),linewidth=0.9,label='theory')
    plt.axhline(y=0, color='black')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$C(t)$')
    #plt.xlim([-0.001,0.05])
    plt.legend()
    plt.savefig(out_folder + '/chunks/time_corr_nx=%d_tau=%f_lambda=%f_chunk=%d.pdf' % (nx, tau, l, i))
    plt.yscale('log')
    plt.ylim([10**(-3),10**0])
    plt.savefig(out_folder + '/chunks/time_corr_log_nx=%d_tau=%f_lambda=%f_chunk=%d.pdf' % (nx, tau, l, i))
    
    plt.close()
    
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
plt.savefig(out_folder + '/spat_corr_nx=%d_tau=%f_lambda=%f_chunk_avg.pdf' % (nx, tau, l))
#plt.show()
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(out_folder + '/spat_corr_log_nx=%d_tau=%f_lambda=%f_chunk_avg.pdf' % (nx, tau, l))

fig = plt.figure()
plt.plot(np.array(r_list),c_r_alt_avg/c_r_alt_avg[0],'.-',label='simulation')
plt.plot(np.array(r_list),np.exp(-np.array(r_list)/l),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
#plt.xlim([-0.1,7.1])
plt.ylim([-0.05,1.05])
plt.legend()
plt.savefig(out_folder + '/spat_corr_alt_nx=%d_tau=%f_lambda=%f_chunk_avg.pdf' % (nx, tau, l))
#plt.show()
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(out_folder + '/spat_corr_alt_log_nx=%d_tau=%f_lambda=%f_chunk_avg.pdf' % (nx, tau, l))

fig = plt.figure()
plt.plot(c_t_avg[:,0],c_t_avg[:,1]/c_t_avg[0,1],label='simulation')
plt.plot(c_t_avg[:,0],np.exp(-c_t_avg[:,0]/tau),linewidth=0.9,label='theory')
plt.axhline(y=0, color='black')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
plt.xlim([-0.001,0.05])
plt.legend()
plt.savefig(out_folder + '/time_corr_nx=%d_tau=%f_lambda=%f_chunk_avg.pdf' % (nx, tau, l))
plt.yscale('log')
plt.ylim([10**(-3),10**0])
plt.savefig(out_folder + '/time_corr_log_nx=%d_tau=%f_lambda=%f_chunk_avg.pdf' % (nx, tau, l))

np.savetxt(folder + '/Cr_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, tau, l), np.c_[np.array(r_list), c_r_alt_avg/c_r_alt_avg[0], np.exp(-np.array(r_list)/l)])
np.savetxt(folder + '/Ct_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, tau, l), np.c_[c_t_avg[:,0], c_t_avg[:,1]/c_t_avg[0,1], np.exp(-c_t_avg[:,0]/tau)])

#plt.show()