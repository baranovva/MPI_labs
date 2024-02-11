import matplotlib.pyplot as plt
import numpy as np

part = 2
indexes = ['N=1', 'N=2', 'N=3', 'N=4']
threads = np.arange(1, 5)

if part == 1:
    times_openmp_integral = np.array([0.0382028, 0.0272019, 0.017011, 0.0134854])

    plt.scatter(indexes, times_openmp_integral, s=60, c='b')
    plt.ylabel("t, s")
    plt.grid(True)
    plt.show()

    s_openmp_integral = [times_openmp_integral[0] / time for time in times_openmp_integral]

    print(s_openmp_integral)

    plt.scatter(indexes, s_openmp_integral, s=60, c='b')
    plt.ylabel("Sp", fontsize=12)
    plt.grid(True)
    plt.show()

    e = [s_openmp_integral[i] / threads[i] for i in range(len(s_openmp_integral))]

    print(e)

    plt.scatter(indexes, e,  s=60, c='b')
    plt.ylabel("Ep", fontsize=12)
    plt.grid(True)
    plt.show()

elif part == 2:
    times_openmp_array = np.array([0.55303098, 0.29928988, 0.17527864, 0.142195923])

    plt.scatter(indexes, times_openmp_array, s=60, c='g')
    plt.ylabel("t, s")
    plt.grid(True)
    plt.show()

    s_openmp_array = [times_openmp_array[0] / time for time in times_openmp_array]

    print(s_openmp_array)

    plt.scatter(indexes, s_openmp_array, s=60, c='g')
    plt.ylabel("Sp", fontsize=12)
    plt.grid(True)
    plt.show()

    e_openmp_array = [s_openmp_array[i] / threads[i] for i in range(len(s_openmp_array))]

    print(e_openmp_array)

    plt.scatter(indexes, e_openmp_array, s=60, c='g')
    plt.ylabel("Ep", fontsize=12)
    plt.grid(True)
    plt.show()
