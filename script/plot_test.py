###################################################################################################
### Written by : Eric Tan
###
### Program: plot_test.py
### Purpose: Generate test plots outputed by test program. Used to verify the correctness of
###          of Ising, Clock, and XY Models.
###################################################################################################


import matplotlib.pyplot as plt
import numpy as np
import os
import shutil


src_dir = "../plot/data/"

for filename in os.listdir("../"):
    if filename.endswith(".txt"):
        shutil.move("../" + filename, src_dir)

# Load all 2D data. Note that T is the same for all plots so it needs to be loaded once
T, S_ising2_clean    = np.loadtxt(src_dir + "2D_ising_clean.txt", unpack=True)
S_ising2_disorder    = np.loadtxt(src_dir + "2D_ising_disorder.txt", usecols=(1,))
S_clock2_2_clean     = np.loadtxt(src_dir + "2D_clock_clean_q=2.txt", usecols=(1,))
S_clock2_2_disorder  = np.loadtxt(src_dir + "2D_clock_disorder_q=2.txt", usecols=(1,))
S_clock2_20_clean    = np.loadtxt(src_dir + "2D_clock_clean_q=20.txt", usecols=(1,))
S_clock2_20_disorder = np.loadtxt(src_dir + "2D_clock_disorder_q=20.txt", usecols=(1,))
S_xy2_clean          = np.loadtxt(src_dir + "2D_xy_clean.txt", usecols=(1,))
S_xy2_disorder       = np.loadtxt(src_dir + "2D_xy_disorder.txt", usecols=(1,))

fig1 = plt.figure(figsize=(20, 16), dpi=80, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size':12})

ax1 = fig1.add_subplot(221)
ax1.plot(T, S_ising2_clean)
ax1.plot(T, S_ising2_disorder)
ax1.set_xlim([0, 5])
ax1.set_ylim([0, 1])
ax1.set_xlabel(r"T [J / $k_B$]")
ax1.set_ylabel(r"S [$k_B$]")
ax1.set_title("Entropy of 2D Ising Model")
ax1.legend(["Clean system", r"Disorder system with $\Delta = 5$"])


ax2 = fig1.add_subplot(222)
ax2.plot(T, S_ising2_clean)
ax2.plot(T, S_ising2_disorder)
ax2.plot(T, S_clock2_2_clean)
ax2.plot(T, S_clock2_2_disorder)
ax2.set_ylim([0, 1])
ax2.set_xlim([0, 5])
ax2.set_xlabel(r"T [J / $k_B$]")
ax2.set_ylabel(r"S [$k_B$]")
ax2.set_title("Comparsion of Entropy of 2D Ising and Clock (2 spin)")
ax2.legend(["Ising Clean", r"Ising Disorder ($\Delta = 5$)",
         "Clock Clean (2 spins)", r"Clock Disorder ($\Delta = 5$)"])


ax3 = fig1.add_subplot(223)
ax3.plot(T, S_clock2_20_clean)
ax3.plot(T, S_clock2_20_disorder)
ax3.set_ylim([1.5, 3.0])
ax3.set_xlim([0, 5])
ax3.set_xlabel(r"T [J / $k_B$]")
ax3.set_ylabel(r"S [$k_B$]")
ax3.set_title("Entropy of 2D Clock Model With 20 Spins")
ax3.legend(["Clean ", r"Disorder With $\Delta = 5$"])


ax4 = fig1.add_subplot(224)
ax4.plot(T, S_xy2_clean)
ax4.plot(T, S_xy2_disorder)
ax4.set_ylim([2.5, 5.0])
ax4.set_xlim([0, 5])
ax4.set_xlabel(r"T [J / $k_B$]")
ax4.set_ylabel(r"S [$k_B$]")
ax4.set_title("Entropy of 2D XY Model")
ax4.legend(["Clean ", r"Disorder With $\Delta = 5$"])

fig1.savefig('../plot/test_plots.pdf', bbox_inches='tight')
