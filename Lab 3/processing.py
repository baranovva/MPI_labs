import matplotlib.pyplot as plt
import numpy as np

part = 2
indexes = ['N=1', 'N=2', 'N=3', 'N=4']
threads = np.arange(1, 5)

if part == 1:
    times_openmp_integral = np.array([1.2946848, 0.68593885, 0.47992434, 0.42186248])
    times_mpi_integral = np.array([1.399542, 0.732657, 0.539793, 0.438021])

    plt.scatter(indexes, times_openmp_integral, label='Integral, OpenMP', s=60, c='b')
    plt.scatter(indexes, times_mpi_integral, label='Integral, MPI', s=60, c='r')
    plt.ylabel("t, sec")
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

elif part == 2:
    times_openmp_array = np.array([3.115, 1.6588, 1.2798, 1.0856])
    times_mpi_array = np.array([3.53399, 1.81284, 1.62305, 1.44017])

    plt.scatter(indexes, times_openmp_array, label='Array, OpenMP', s=60, c='g')
    plt.scatter(indexes, times_mpi_array, label='Array, MPI', s=60, c='#9467bd')
    plt.ylabel("t, sec")
    plt.legend()
    plt.grid(True)
    plt.show()

    s_openmp_array = [times_openmp_array[0] / time for time in times_openmp_array]
    s_mpi_array = [times_mpi_array[0] / time for time in times_mpi_array]

    print(s_openmp_array)
    print(s_mpi_array)

    plt.scatter(indexes, s_openmp_array, label='Array, OpenMP', s=60, c='g')
    plt.scatter(indexes, s_mpi_array, label='Array, MPI', s=60, c='#9467bd')
    plt.ylabel("Sp", fontsize=12)
    plt.legend()
    plt.grid(True)
    plt.show()

    e_openmp_array = [s_openmp_array[i] / threads[i] for i in range(len(s_openmp_array))]
    e_mpi_array = [s_mpi_array[i] / threads[i] for i in range(len(s_mpi_array))]

    print(e_openmp_array)
    print(e_mpi_array)

    plt.scatter(indexes, e_openmp_array, label='Array, OpenMP', s=60, c='g')
    plt.scatter(indexes, e_mpi_array, label='Array, MPI', s=60, c='#9467bd')
    plt.ylabel("Ep", fontsize=12)
    plt.legend()
    plt.grid(True)
    plt.show()
