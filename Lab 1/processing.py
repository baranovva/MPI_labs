import matplotlib.pyplot as plt
import numpy as np

processes = np.arange(1, 5)

times = np.array([1.2946848, 0.68593885, 0.47992434, 0.42186248])

plt.scatter(processes, times, s=60, c='b')
plt.ylabel("t, sec", fontsize=12)
plt.grid(True)
plt.show()

s = [times[0] / time for time in times]

plt.scatter(processes, s, s=60, c='b')
plt.ylabel("Sp", fontsize=12)
plt.grid(True)
plt.show()

e = [s[i] / processes[i] for i in range(len(s))]

plt.scatter(processes, e, s=60, c='b')
plt.ylabel("Ep", fontsize=12)
plt.grid(True)
plt.show()
