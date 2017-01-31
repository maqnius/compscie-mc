import numpy as np
cimport numpy as np
from libc.math cimport sqrt

DTYPE = np.int
ctypedef np.int_t DTYPE_t

def calc_k_vectors(int K):
    cdef int b_limit, c_limit

    k_vectors = []

    # Create all k-vectors with absolute value <= K
    for a in range(-K, K + 1):
        b_limit = int(np.sqrt(K ** 2 - a ** 2))
        for b in range(-b_limit, b_limit + 1):
            c_limit = int(np.sqrt(K ** 2 - a ** 2 - b ** 2))
            for c in range(-c_limit, c_limit + 1):
                k_vectors.append([a, b, c])

    # Remove k = [0, 0, 0]
    k_vectors.remove([0, 0, 0])

    return np.array(k_vectors)

def calc_k_vectors_test(int K):
    cdef int b_limit, c_limit, a, b, c

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
    print(i/k_vectors.shape[0])
    k_vectors = k_vectors[:i]
    k_vectors = np.delete(k_vectors, i//2, 0)

    return k_vectors