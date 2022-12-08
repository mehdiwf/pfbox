import matplotlib.pyplot as plt
import numpy as np

rust_time_list = np.array([0.636702614, 3.105668876, 6.165691205, 30.774588656, 61.782680750])
c_time_list = np.array([1.605196, 8.001673, 15.999283, 80.354180, 159.692652])
nb_step_list = np.array([10000, 50000, 100000, 500000, 1000000])

fig, axes= plt.subplots(2,1)

axes[0].scatter(nb_step_list, rust_time_list, c="orange",
                marker="+",
                label="rust")
axes[0].scatter(nb_step_list, c_time_list, c="grey",
                marker="+",
                label="c")
axes[0].set_xlabel("total step number")
axes[0].set_ylabel("execution time in seconds")
axes[0].set_xscale("log")
axes[0].legend()

axes[1].scatter(nb_step_list, c_time_list/rust_time_list,
                marker="+",
                c="red", label="c execution time / rust execution time")
axes[1].set_xlabel("total step number")
axes[1].set_ylabel("speed multiplication")
axes[1].set_xscale("log")
axes[1].set_ylim(2.45, 2.7)
axes[1].legend()
fig.savefig("rust_c_comparison.png")
plt.show()
