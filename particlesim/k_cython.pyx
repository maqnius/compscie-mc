#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen,
#                       Simon Treu, Hanna Wulkow, Thomas Hadler
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#


import numpy as np
import cython
cimport numpy as np
from libc.math cimport sqrt

DTYPE = np.int
ctypedef np.int_t DTYPE_t

def calc_k_vectors(int K):
    """
    Calculates the k vectors of our lattice until a cutoff Value K.

    Parameters
    ----------
    K : int
        Cutoff value for the absolute value of the k-vectors


    Returns
    -------
    k_vectors : ndarray
        Array of k vectors that have an absolute value below cutoff K
    """
    cdef int b_limit, c_limit

    k_vectors = []

    # Create all k-vectors with absolute value <= K
    for a in range(-K, K + 1):
        b_limit = int((K ** 2 - a ** 2)**.5)
        for b in range(-b_limit, b_limit + 1):
            c_limit = int((K ** 2 - a ** 2 - b ** 2)**.5)
            for c in range(-c_limit, c_limit + 1):
                k_vectors.append([a, b, c])

    # Remove k = [0, 0, 0]
    k_vectors.remove([0, 0, 0])

    return np.array(k_vectors)

def calc_k_vectors_test(int K):
    """
    Function for testing and speedup

    Parameters
    ----------
    K : int
        cutoff radius in reciprocal space
    """
    cdef int b_limit, c_limit, a, b, c

    # Array size can be estimated by the ratio of the Volume of a sphare compared to a cube which is ~ 0.53
    cdef int max_number = int(4.2 * K**3)
    cdef np.ndarray[DTYPE_t, ndim=2] k_vectors = np.empty((max_number, 3), dtype=DTYPE)
    # Create all k-vectors with absolute value <= K
    cdef int i = 0
    for a in range(-K, K + 1):
        b_limit = int(sqrt(K ** 2 - a ** 2))
        for b in range(-b_limit, b_limit + 1):
            c_limit = int(sqrt(K ** 2 - a ** 2 - b ** 2))
            for c in range(-c_limit, c_limit + 1):
                k_vectors[i] = [a, b, c]
                i += 1

    k_vectors = k_vectors[:i]
    k_vectors = np.delete(k_vectors, i//2, 0) # Remove [0,0,0]

    return k_vectors

@cython.boundscheck(False)
@cython.wraparound(False)
def fast_distances(double[:, :] xyz, double box_len, double[:,:] distances):
    """
    Calculates the distance of all positions to all other positions in the xyz array. It says doubles here but
    python floats are cython doubles

    Parameters
    ----------
    xyz :       double[:,:]
                position array
    box_len :    double
                 the length of the box for the periodic boundry
    distances :  double[:,:]
                 an array passed to the function to avoid the necessity of returning a new array

    Return
    ------
        void
    """

    cdef:
        Py_ssize_t nrow = xyz.shape[0]
        Py_ssize_t ncol = xyz.shape[1]
        int i, j
        double vec, dist_tmp, distance
        double box_half = box_len / 2.0

    for i in range(nrow):
        for j in range(i+1, nrow):
            distance = 0
            for k in range(ncol):
                vec = xyz[i,k]-xyz[j,k]
                dist_tmp  = box_half - ((vec+box_half) % box_len)
                distance += dist_tmp*dist_tmp
            distance = sqrt(distance)
            distances[i,j] = distance
            distances[j,i] = distance

    return

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fast_distances_cpdef(double[:, :] xyz, double box_len, double[:,:] distances):
    """
    Calculates the distance of all positions to all other positions in the xyz array. It says doubles here but
    python floats are cython doubles

    Parameters
    ----------
    xyz :       double[:,:]
                position array
    box_len :    double
                 the length of the box for the periodic boundry
    distances :  double[:,:]
                 an array passed to the function to avoid the necessity of returning a new array

    Return
    ------
        void
    """

    cdef:
        Py_ssize_t nrow = xyz.shape[0]
        Py_ssize_t ncol = xyz.shape[1]
        int i, j
        double vec, dist_tmp, distance
        double box_half = box_len / 2.0

    for i in range(nrow):
        for j in range(i+1, nrow):
            distance = 0
            for k in range(ncol):
                vec = xyz[i,k]-xyz[j,k]
                dist_tmp  = box_half - ((vec+box_half) % box_len)
                distance += dist_tmp*dist_tmp
            distance = sqrt(distance)
            distances[i,j] = distance
            distances[j,i] = distance

    return