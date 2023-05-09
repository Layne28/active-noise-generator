import numpy as np
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

print(c_t[0])

fig = plt.figure()
plt.plot(times,c_t,label='simulation')
plt.plot(times,D*np.exp(-times/tau),linewidth=0.9,label='theory') #include tau in exponential
plt.axhline(y=0, color='black', linestyle='--')
plt.xlabel(r'$t$')
plt.ylabel(r'$C(t)$')
plt.xlim([0,5])
plt.legend()
plt.savefig('./time_corr.png')
#plt.yscale('log')
#plt.ylim([10**(-3),10**0])
#plt.savefig(folder + '/time_corr_log.pdf')

plt.show()