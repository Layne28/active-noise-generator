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

c_x_avg = c_r_3d_avg[c_r_3d_avg[:,1]==0]
c_x_avg = c_x_avg[c_x_avg[:,2]==0]
c_x_avg = c_x_avg[c_x_avg[:,0]>=0]
c_x_avg = np.delete(c_x_avg, 1, 1)
c_x_avg = np.delete(c_x_avg, 1, 1)
c_x_avg = c_x_avg[c_x_avg[:,0].argsort()]
print(c_x_avg)
print(c_x_avg.shape)

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
    
    c_x = c_r_3d[c_r_3d[:,1]==0]
    c_x = c_x[c_x[:,2]==0]
    c_x = c_x[c_x[:,0]>=0]
    c_x = np.delete(c_x, 1, 1)
    c_x = np.delete(c_x, 1, 1)
    c_x = c_x[c_x[:,0].argsort()]

    if i>0:
        c_r_avg += c_r
        c_t_avg += c_t
        c_x_avg += c_x
        
    for j in range(c_r_3d.shape[0]):
        r = np.sqrt(c_r_3d[j,0]**2 + c_r_3d[j,1]**2 + c_r_3d[j,2]**2)
        index = r_list.index(r)
        c_r_alt[index] += c_r_3d[j,3]
        r_hist[index] += 1.0
        
    c_r_alt /= r_hist
    c_r_alt_avg += c_r_alt    
    
c_r_avg /= nchunks
c_r_alt_avg /= nchunks
c_t_avg /= nchunks
c_x_avg /= nchunks

np.savetxt(folder + '/Cr_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, tau, l), np.c_[np.array(r_list), c_r_alt_avg/c_r_alt_avg[0], np.exp(-np.array(r_list)/l)])
np.savetxt(folder + '/Cx_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, tau, l), c_x_avg)
np.savetxt(folder + '/Ct_total_nx=%d_tau=%f_lambda=%f.txt' % (nx, tau, l), np.c_[c_t_avg[:,0], c_t_avg[:,1]/c_t_avg[0,1], np.exp(-c_t_avg[:,0]/tau)])

#plt.show()