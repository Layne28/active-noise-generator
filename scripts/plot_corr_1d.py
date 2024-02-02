import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import h5py
import sys

myfile = sys.argv[1]

data = h5py.File(myfile)

c_t = np.array(data['/corr_t/value'])
times = np.array(data['/corr_t/time'])
tau = np.array(data['/parameters/tau'])
Lambda = np.array(data['/parameters/lambda'])
D = np.array(data['/parameters/D'])
nx = np.array(data['/grid/dimensions'])

c_r = np.array(data['/corr_r/value'])
positions = np.array(data['/corr_r/position'])

r = np.zeros(c_r.shape)
for i in range(c_r.shape[0]):
    r[i] = la.norm(positions[i])

print(c_t[0])

fig = plt.figure()
plt.plot(times,c_t,label='simulation')
plt.plot(times,D*np.exp(-times/tau),linewidth=0.9,label='theory') #include tau in exponential
plt.axhline(y=0, color='black', linestyle='--')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
plt.xlim([0,5*tau])
plt.legend()
plt.savefig('./plots/1d/time_corr_nx=%d_tau=%f_lambda=%f.png' % (nx, tau, Lambda))

fig = plt.figure()
sim = plt.plot(positions, c_r, label='simulation')
theory = plt.plot(positions, np.exp(-np.abs(positions)/Lambda), label='theory')
plt.xlabel(r'$x$')
plt.ylabel(r'$C(x)$')
#plt.xlim([0,5])
plt.legend()
plt.tight_layout()
plt.savefig('./plots/1d/spat_corr_nx=%d_tau=%f_lambda=%f.png' % (nx, tau, Lambda))

plt.show()