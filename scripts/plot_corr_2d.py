import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import h5py
import sys

myfile = sys.argv[1]

data = h5py.File(myfile,'r')

c_t = np.array(data['/corr_t/value'])
times = np.array(data['/corr_t/time'])
tau = np.array(data['/parameters/tau'])
Lambda = np.array(data['/parameters/lambda'])
D = np.array(data['/parameters/D'])
n = np.array(data['/grid/dimensions'])

c_r = np.array(data['/corr_r/value'])
positions = np.array(data['/corr_r/position'])

r = np.zeros(c_r.shape)
for i in range(c_r.shape[0]):
    for j in range(c_r.shape[1]):
        r[i,j] = la.norm(positions[i,j,:])

print(c_t[0])

fig = plt.figure()
plt.plot(times,c_t,label='simulation')
plt.plot(times,D*np.exp(-times/tau),linewidth=0.9,label='theory') #include tau in exponential
plt.axhline(y=0, color='black', linestyle='--')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
#plt.xlim([0,5])
plt.legend()
plt.savefig('./plots/2d/time_corr_nx=%d_ny=%d_tau=%f_lambda=%f.png' % (n[0], n[1], tau, Lambda))

fig, ax = plt.subplots(1,2)
sim = ax[0].imshow(c_r,vmin=0,vmax=1)#, extent=positions)
plt.colorbar(sim, ax=ax[0])

theory = ax[1].imshow(np.exp(-r/Lambda),vmin=0,vmax=1)
plt.colorbar(theory, ax=ax[1])
#plt.plot(times,D*np.exp(-times/tau),linewidth=0.9,label='theory') #include tau in exponential
ax[0].set_xlabel(r'$x$')
ax[1].set_xlabel(r'$x$')
ax[0].set_ylabel(r'$y$')

ax[0].set_title(r'simulation')
ax[1].set_title(r'theory')
#plt.xlim([0,5])
#plt.legend()
plt.tight_layout()
plt.savefig('./plots/2d/spat_corr_nx=%d_ny=%d_tau=%f_lambda=%f.png' % (n[0], n[1], tau, Lambda))


plt.show()