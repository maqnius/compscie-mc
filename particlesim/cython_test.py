import pyximport; pyximport.install()
import particlesim.k_cython as k_cython

print(k_cython.calc_k_vectors(2))


