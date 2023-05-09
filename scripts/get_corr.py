#Compute space and time correlation functions of active noise

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la
import h5py
import numba

def main():

    in_folder = sys.argv[1]
    out_folder = sys.argv[2]

    data = h5py.File(in_folder + '/noise_traj.h5', 'r')

    tau_max = 5 #max time to which to compute correlation function

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    times = np.array(data['/noise/time'])
    noise_x = np.array(data['/noise/value/x'])
    noise_y = np.array(data['/noise/value/y'])
    dims = np.array(data['/grid/dimensions'])
    dx = np.array(data['/grid/spacing'])
    tau = np.array(data['/parameters/tau'])
    Lambda = np.array(data['/parameters/lambda'])
    D = np.array(data['/parameters/D'])

    data.close()

    noise = np.stack([noise_x, noise_y], axis=-1)
    print(noise.shape)

    delta_t = times[1]-times[0]
    frame_diff_max = int(tau_max/delta_t) #max number of data points to which to compute correlation function
    print(frame_diff_max)


    print('Loaded data.')

    #Create output file
    corr_file = h5py.File(out_folder + '/noise_corr.h5', 'w')
    corr_file.create_dataset('/corr_t/time', data=delta_t*np.arange(frame_diff_max))
    #corr_file.create_dataset('/corr_r/displacements', data=displacements)
    corr_file.create_dataset('/grid/dimensions', data=dims)
    corr_file.create_dataset('/grid/spacing', data=dx)
    corr_file.create_dataset('/parameters/tau', data=tau)
    corr_file.create_dataset('/parameters/lambda', data=Lambda)
    corr_file.create_dataset('/parameters/D', data=D)

    #Compute C(t)
    print('c0_test: ', get_corr_t0_test(noise, frame_diff_max))
    c_t = get_corr_t(noise, frame_diff_max)
    print('time zero corr: ', c_t[0])

    corr_file.create_dataset('/corr_t/value', data=c_t) 

    #Compute C(r)
    #corr_file.create_dataset('/corr_q/corr', data=c_q)

@numba.jit(nopython=True) 
def get_corr_t(xi_mat, tmax):
    nx = xi_mat.shape[1]
    ny = xi_mat.shape[2]
    ndim = xi_mat.shape[3]
    print(ndim)
    c_t = np.zeros(tmax)
    for t in range(tmax):
        #take inner product of xi_r(t) and xi_r(t+delta_t) and average over t and r
        for i in range(nx): #sum over x
            for j in range(ny): #sum over y
                for d in range(ndim): #sum over dimensions
                    c_t[t] += np.mean(xi_mat[:-tmax,i,j,d]*xi_mat[t:(-tmax+t),i,j,d])
        c_t[t] /= (nx*ny*ndim)
    return c_t

@numba.jit(nopython=True) 
def get_corr_t0_test(xi_mat, tmax):
    nx = xi_mat.shape[1]
    ny = xi_mat.shape[2]
    ndim = xi_mat.shape[3]
    print(ndim)
    c_0 = 0.0
    for i in range(nx): #sum over x
        for j in range(ny): #sum over y
            for d in range(ndim): #sum over dimensions
                c_0 += xi_mat[0,i,j,d]*xi_mat[0,i,j,d]
    c_0 /= (nx*ny*ndim)
    return c_0

main()