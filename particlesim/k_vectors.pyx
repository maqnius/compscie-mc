def calc_k_vectors(int K):


    import numpy as np
    cdef int b_limit, c_limit


    k_vectors = []
    # Create all k-vectors with absolute value <= K
    for a in range(-K, K + 1):
        b_limit = int(np.sqrt(K ** 2 - a ** 2))
        for b in range(-b_limit, b_limit + 1):
            c_limit = int(np.sqrt(K ** 2 - a ** 2 - b ** 2))
            for c in range(-c_limit, c_limit + 1):
                k_vectors.append([a, b, c])

    return k_vectors