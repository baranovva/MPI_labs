import matplotlib.pyplot as plt
import numpy as np


def s(times):
    temp = [times[0] / time for time in times]
    print(temp)
    return temp


def e(s_):
    threads = np.arange(1, 5)
    temp = [s_[i] / threads[i] for i in range(len(s_))]
    print(temp)
    return temp


indexes = ['N=1', 'N=2', 'N=3', 'N=4']

times_openmp_j_41 = np.array([0.3504, 0.3307, 0.3041, 0.2805])
times_openmp_j_71 = np.array([1.5075, 1.0316, 0.9943, 0.9564])
times_openmp_j_101 = np.array([5.3585, 4.0301, 3.0320, 2.7862])
times_openmp_gz_41 = np.array([0.3708, 0.3210, 0.3029, 0.2845])
times_openmp_gz_71 = np.array([0.8820, 0.7320, 0.5602, 0.5438])
times_openmp_gz_101 = np.array([2.6631, 1.6667, 1.4895, 1.3768])
times_mpi_j_41 = np.array([0.91596, 0.50972, 0.380319, 0.311405])
times_mpi_j_71 = np.array([2.9545, 1.6109, 1.1865, 0.9511])
times_mpi_j_101 = np.array([5.9876, 3.2926, 2.4049, 1.9452])

Times = [times_openmp_j_41, times_openmp_j_71, times_openmp_j_101,
         times_openmp_gz_41, times_openmp_gz_71, times_openmp_gz_101,
         times_mpi_j_41, times_mpi_j_71, times_mpi_j_101]

s_openmp_j_41 = s(Times[0])
s_openmp_j_71 = s(Times[1])
s_openmp_j_101 = s(Times[2])
s_openmp_gz_41 = s(Times[3])
s_openmp_gz_71 = s(Times[4])
s_openmp_gz_101 = s(Times[5])
s_mpi_j_41 = s(Times[6])
s_mpi_j_71 = s(Times[7])
s_mpi_j_101 = s(Times[8])

print('=' * 10)

e_openmp_j_41 = e(s_openmp_j_41)
e_openmp_j_71 = e(s_openmp_j_71)
e_openmp_j_101 = e(s_openmp_j_101)
e_openmp_gz_41 = e(s_openmp_gz_41)
e_openmp_gz_71 = e(s_openmp_gz_71)
e_openmp_gz_101 = e(s_openmp_gz_101)
e_mpi_j_41 = e(s_mpi_j_41)
e_mpi_j_71 = e(s_mpi_j_71)
e_mpi_j_101 = e(s_mpi_j_101)

# 41x41
plt.scatter(indexes, Times[0], label='Jacobi, OpenMP, 41x41', s=60, c='b')
plt.scatter(indexes, Times[3], label='Gauss–Seidel, OpenMP, 41x41', s=60, c='r')
plt.scatter(indexes, Times[6], label='Jacobi, MPI, 41x41', s=60, c='g')
plt.ylabel("t, s")
plt.legend()
plt.grid(True)
plt.show()

plt.scatter(indexes, s_openmp_j_41, label='Jacobi, OpenMP, 41x41', s=60, c='b')
plt.scatter(indexes, s_openmp_gz_41, label='Gauss–Seidel, OpenMP, 41x41', s=60, c='r')
plt.scatter(indexes, s_mpi_j_41, label='Jacobi, MPI, 41x41', s=60, c='g')
plt.ylabel("Sp", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

plt.scatter(indexes, e_openmp_j_41, label='Jacobi, OpenMP, 41x41', s=60, c='b')
plt.scatter(indexes, e_openmp_gz_41, label='Gauss–Seidel, OpenMP, 41x41', s=60, c='r')
plt.scatter(indexes, e_mpi_j_41, label='Jacobi, MPI, 41x41', s=60, c='g')
plt.ylabel("Ep", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

# 71x71
plt.scatter(indexes, Times[1], label='Jacobi, OpenMP, 71x71', s=60, c='b')
plt.scatter(indexes, Times[4], label='Gauss–Seidel, OpenMP, 71x71', s=60, c='r')
plt.scatter(indexes, Times[7], label='Jacobi, MPI, 71x71', s=60, c='g')
plt.ylabel("t, s")
plt.legend()
plt.grid(True)
plt.show()

plt.scatter(indexes, s_openmp_j_71, label='Jacobi, OpenMP, 71x71', s=60, c='b')
plt.scatter(indexes, s_openmp_gz_71, label='Gauss–Seidel, OpenMP, 71x71', s=60, c='r')
plt.scatter(indexes, s_mpi_j_71, label='Jacobi, MPI, 71x71', s=60, c='g')
plt.ylabel("Sp", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

plt.scatter(indexes, e_openmp_j_71, label='Jacobi, OpenMP, 71x71', s=60, c='b')
plt.scatter(indexes, e_openmp_gz_71, label='Gauss–Seidel, OpenMP, 71x71', s=60, c='r')
plt.scatter(indexes, e_mpi_j_71, label='Jacobi, MPI, 71x71', s=60, c='g')
plt.ylabel("Ep", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

# 101x101
plt.scatter(indexes, Times[2], label='Jacobi, OpenMP, 101x101', s=60, c='b')
plt.scatter(indexes, Times[5], label='Gauss–Seidel, OpenMP, 101x101', s=60, c='r')
plt.scatter(indexes, Times[8], label='Jacobi, MPI, 101x101', s=60, c='g')
plt.ylabel("t, s")
plt.legend()
plt.grid(True)
plt.show()

plt.scatter(indexes, s_openmp_j_101, label='Jacobi, OpenMP, 101x101', s=60, c='b')
plt.scatter(indexes, s_openmp_gz_101, label='Gauss–Seidel, OpenMP, 101x101', s=60, c='r')
plt.scatter(indexes, s_mpi_j_101, label='Jacobi, MPI, 101x101', s=60, c='g')
plt.ylabel("Sp", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

plt.scatter(indexes, e_openmp_j_101, label='Jacobi, OpenMP, 101x101', s=60, c='b')
plt.scatter(indexes, e_openmp_gz_101, label='Gauss–Seidel, OpenMP, 101x101', s=60, c='r')
plt.scatter(indexes, e_mpi_j_101, label='Jacobi, MPI, 101x101', s=60, c='g')
plt.ylabel("Ep", fontsize=12)
plt.legend()
plt.grid(True)
plt.show()

cnt = np.array([1411, 6440, 13114, 825, 3580, 7271, 170, 467, 408])
cnt_sec = [[cnt[i] / time for time in Times[i]] for i in range(len(cnt))]
print([cnt_sec[i] for i in range(len(cnt_sec))])

labels = ['41x41', '71x71', '101x101']
bar_width = 0.25

fig, ax = plt.subplots()
bar1 = ax.bar([i - bar_width for i in range(3)], [cnt[i] for i in range(3)], bar_width, label='Jacobi, OpenMP')
bar2 = ax.bar([i for i in range(3)], [cnt[i] for i in range(3, 6)], bar_width, label='Gauss–Seidel, OpenMP')
bar3 = ax.bar([i + bar_width for i in range(3)], [cnt[i] for i in range(6, 9)], bar_width, label='Jacobi, MPI')
ax.set_xticks(range(3))
ax.set_xticklabels(labels)
ax.set_ylabel('Сетка')
ax.legend()
ax.set_yscale('log')
ax.set_ylabel('Количество итераций')
plt.grid(True)
plt.show()

fig, ax = plt.subplots()
bar1 = ax.bar([i - bar_width for i in range(3)], [cnt_sec[i][3] for i in range(3)], bar_width, label='Jacobi, OpenMP')
bar2 = ax.bar([i for i in range(3)], [cnt_sec[i][3] for i in range(3, 6)], bar_width, label='Gauss–Seidel, OpenMP')
bar3 = ax.bar([i + bar_width for i in range(3)], [cnt_sec[i][3] for i in range(6, 9)], bar_width, label='Jacobi, MPI')
ax.set_xticks(range(3))
ax.set_xticklabels(labels)
ax.set_ylabel('Сетка')
ax.legend()
ax.set_yscale('log')
ax.set_ylabel('Количество итераций в секунду')
plt.grid(True)
plt.show()
