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

    in_file = sys.argv[1] #expects noise_traj.h5
    out_folder = "/".join(in_file.split('/')[:-1])
    print(out_folder)

    data = h5py.File(in_file, 'r')

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    times = np.array(data['/noise/time'])
    noise_x = np.array(data['/noise/value/x'])
    dims = np.array(data['/grid/dimensions'])
    dx = np.array(data['/grid/spacing'])
    tau = np.array(data['/parameters/tau'])
    Lambda = np.array(data['/parameters/lambda'])
    D = np.array(data['/parameters/D'])

    tau_max = 5*tau #max time to which to compute correlation function

    data.close()

    noise = noise_x

    if times.shape[0]>1:
        delta_t = times[1]-times[0]
    else:
        delta_t = 1
    frame_diff_max = int(tau_max/delta_t) #max number of data points to which to compute correlation function

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
    ct0 = get_corr_t0_test(noise)
    print('initial frame correlation: ', ct0)
    c_t = get_corr_t(noise, frame_diff_max)
    print('time zero corr: ', c_t[0])

    corr_file.create_dataset('/corr_t/value', data=c_t) 

    #Compute C(r)
    c_r, pos = get_corr_r(noise, dx)
    #print(pos)
    corr_file.create_dataset('/corr_r/value', data=c_r)
    corr_file.create_dataset('/corr_r/position', data=pos)

@numba.jit(nopython=True) 
def get_corr_t(xi_mat, tmax):
    nx = xi_mat.shape[1]
    c_t = np.zeros(tmax)
    for t in range(tmax):
        #take inner product of xi_r(t) and xi_r(t+delta_t) and average over t and r
        for i in range(nx): #sum over x
            c_t[t] += np.mean(xi_mat[:-tmax,i]*xi_mat[t:(-tmax+t),i])
        c_t[t] /= nx
    return c_t

@numba.jit(nopython=True) 
def get_corr_r(xi_mat, dx):
    tmax = xi_mat.shape[0]
    nx = xi_mat.shape[1]

    c_r = np.zeros(nx)
    pos = np.zeros(nx)
    counts = np.zeros(nx)

    for i1 in range(nx):
        for i2 in range(nx):

            #get "distance" on grid
            xind = i2-i1

            #apply pbc
            if xind<(-nx//2):
                xind += nx
            if xind>=(nx//2):
                xind -= nx

            pos[xind+nx//2] = dx[0]*xind
            counts[xind+nx//2] += 1

    print('r=0 count: ', counts[nx//2])

    #for t in range(tmax):
    for i1 in range(nx):
        for i2 in range(nx):
                    
            #get "distance" on grid
            xind = i2-i1

            #apply pbc
            if xind<(-nx//2):
                xind += nx
            if xind>=(nx//2):
                xind -= nx

            c_r[xind+nx//2] += np.mean(xi_mat[:,i1]*xi_mat[:,i2])/counts[xind+nx//2]

    print(c_r[nx//2])
    print(c_r[nx//2-1])
    print(np.max(c_r))
    
    return c_r, pos
                                                  
@numba.jit(nopython=True) 
def get_corr_t0_test(xi_mat):
    nx = xi_mat.shape[1]
    ndim = 1
    c_0 = 0.0
    for i in range(nx): #sum over x
        c_0 += xi_mat[0,i]*xi_mat[0,i]
    c_0 /= nx
    return c_0

main()