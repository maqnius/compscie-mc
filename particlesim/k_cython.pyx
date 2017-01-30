def calc_k_vectors(int K):


    import numpy as np
    cdef int b_limit, c_limit


    k_vectors_part = []
    # Create all k-vectors with absolute value <= K
    for a in range(0, K + 1):
        b_limit = int(np.sqrt(K ** 2 - a ** 2))
        for b in range(0, b_limit + 1):
            c_limit = int(np.sqrt(K ** 2 - a ** 2 - b ** 2))
            for c in range(0, c_limit + 1):
                k_vectors_part.append([a, b, c])
    k_vectors_part.pop(0)
    k_vectors_part = np.array(k_vectors_part)

    permutations = [[-1, -1, -1], [-1, -1, 1], [-1, 1, 1], [1, -1, -1], [1, 1, -1], [1, -1, 1], [-1, 1, -1]]

    k_vectors = k_vectors_part
    for permutation in permutations:
        k_vectors = np.concatenate((k_vectors_part, np.multiply(k_vectors_part, permutation)))

    return k_vectors