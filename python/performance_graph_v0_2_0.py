import matplotlib.pyplot as plt
import numpy as np

rust_personalPC = np.array([0.636702614, 3.105668876, 6.165691205, 30.774588656, 61.782680750])
c_personalPC = np.array([1.605196, 8.001673, 15.999283, 80.354180, 159.692652])
nb_step_list = np.array([10000, 50000, 100000, 500000, 1000000])
rust_jemac = np.array([1.363861436, 6.623272585, 13.228002208, 65.768119034, 131.536727645])
c_jemac = np.array([3.829171, 19.006658, 38.029916, 190.004333, 380.368056])

fig, axes= plt.subplots(2,1)

plot_jemac_personalPC_compar = False

if plot_jemac_personalPC_compar:
    info_personalPC = " personal PC"
    info_jemac = " jemac"
    ratio_max = 3.0
else:
    info_personalPC = ""
    info_jemac = ""
    ratio_max = 2.8

axes[0].scatter(nb_step_list, rust_personalPC, color=(236/255,103/255,28/255),
                marker="+",
                label=f"Rust{info_personalPC}")
axes[0].scatter(nb_step_list, c_personalPC, c="grey",
                marker="+",
                label=f"C{info_personalPC}")
if plot_jemac_personalPC_compar:
    axes[0].scatter(nb_step_list, rust_jemac, color=(236/255,103/255,28/255),
                    marker="d",
                    label=f"Rust{info_jemac}")
    axes[0].scatter(nb_step_list, c_jemac, c="grey",
                    marker="d",
                    label=f"C{info_jemac}")
axes[0].set_xlabel("total step number")
axes[0].set_ylabel("execution time in seconds")
axes[0].set_xscale("log")
axes[0].legend()

axes[1].scatter(nb_step_list, c_personalPC/rust_personalPC,
                marker="+",
                c="blue", label=f"C/Rust{info_personalPC}")
if plot_jemac_personalPC_compar:
    axes[1].scatter(nb_step_list, c_jemac/rust_jemac,
                    marker="+",
                    c="black", label=f"C/Rust{info_jemac}")
axes[1].set_xlabel("total step number")
axes[1].set_ylabel("execution time ratio")
axes[1].set_xscale("log")
axes[1].set_ylim(2.45, ratio_max)
axes[1].legend()
axes[0].set_title("Rust/C comparison on a Van Der Waals phase field code")
fig.savefig("rust_c_comparison.png")
plt.show()
