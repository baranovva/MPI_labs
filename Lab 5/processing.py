import matplotlib.pyplot as plt
import numpy as np

indexes = ['N=1', 'N=2', 'N=3', 'N=4']
threads = np.arange(1, 5)

times_openmp_integral = np.array([0.624834, 0.489784, 0.391744, 0.30332])
times_mpi_integral = np.array([0.622505, 0.513214, 0.455214, 0.404095])

plt.scatter(indexes, times_openmp_integral, label='Integral, OpenMP', s=60, c='b')
plt.scatter(indexes, times_mpi_integral, label='Integral, MPI', s=60, c='r')
plt.ylabel("t, s")
plt.legend()
plt.grid(True)
plt.show()

s_openmp_integral = [times_openmp_integral[0] / time for time in times_openmp_integral]
s_mpi_integral = [times_mpi_integral[0] / time for time in times_mpi_integral]

print(s_openmp_integral)
print(s_mpi_integral)

plt.scatter(indexes, s_openmp_integral, label='Integral, OpenMP', s=60, c='b')
plt.scatter(indexes, s_mpi_integral, label='Integral, MPI', s=60, c='r')
plt.ylabel("Sp", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

e_openmp_integral = [s_openmp_integral[i] / threads[i] for i in range(len(s_openmp_integral))]
e_mpi_integral = [s_mpi_integral[i] / threads[i] for i in range(len(s_mpi_integral))]

print(e_openmp_integral)
print(e_mpi_integral)

plt.scatter(indexes, e_openmp_integral, label='Integral, OpenMP', s=60, c='b')
plt.scatter(indexes, e_mpi_integral, label='Integral, MPI', s=60, c='r')
plt.ylabel("Ep", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()
