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

# Load all data. Note that T is the same for all plots so it needs to be loaded once
T, S_ising_clean   = np.loadtxt(src_dir + "2D_ising_clean.txt", unpack=True)
S_ising_disorder   = np.loadtxt(src_dir + "2D_ising_disorder.txt", usecols=1)
S_clock2_clean     = np.loadtxt(src_dir + "2D_clock_clean_q=2.txt", usecols=1)
S_clock2_disorder  = np.loadtxt(src_dir + "2D_clock_disorder_q=2.txt", usecols=1)
S_clock20_clean    = np.loadtxt(src_dir + "2D_clock_clean_q=20.txt", usecols=1)
S_clock20_disorder = np.loadtxt(src_dir + "2D_clock_disorder_q=20.txt", usecols=1)
S_xy_clean         = np.loadtxt(src_dir + "2D_xy_clean.txt", usecols=1)
S_xy_disorder      = np.loadtxt(src_dir + "2D_xy_disorder.txt", usecols=1)

plt.figure(figsize=(20, 16), dpi=80, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size':12})

plt.subplot(221)
plt.plot(T, S_ising_clean)
plt.plot(T, S_ising_disorder)
plt.ylim([0, 1])
plt.xlim([0, 5])
plt.xlabel(r"T [J / $k_B$]")
plt.ylabel(r"S [$k_B$]")
plt.title("Entropy of 2D Ising Model")
plt.legend(["Clean system", r"Disorder system with $\Delta = 5$"])


plt.subplot(222)
plt.plot(T, S_ising_clean)
plt.plot(T, S_ising_disorder)
plt.plot(T, S_clock2_clean)
plt.plot(T, S_clock2_disorder)
plt.ylim([0, 1])
plt.xlim([0, 5])
plt.xlabel(r"T [J / $k_B$]")
plt.ylabel(r"S [$k_B$]")
plt.title("Comparsion of Entropy of 2D Ising and Clock (2 spin)")
plt.legend(["Ising Clean", r"Ising Disorder ($\Delta = 5$)",
            "Clock Clean (2 spins)", r"Clock Disorder ($\Delta = 5$)"])

plt.subplot(223)
plt.plot(T, S_clock20_clean)
plt.plot(T, S_clock20_disorder)
plt.ylim([1.5, 3.0])
plt.xlim([0, 5])
plt.xlabel(r"T [J / $k_B$]")
plt.ylabel(r"S [$k_B$]")
plt.title("Entropy of 2D Clock Model With 20 Spins")
plt.legend(["Clean ", r"Disorder With $\Delta = 5$"])


plt.subplot(224)
plt.plot(T, S_xy_clean)
plt.plot(T, S_xy_disorder)
plt.ylim([2.5, 5.0])
plt.xlim([0, 5])
plt.xlabel(r"T [J / $k_B$]")
plt.ylabel(r"S [$k_B$]")
plt.title("Entropy of 2D XY Model")
plt.legend(["Clean ", r"Disorder With $\Delta = 5$"])

plt.savefig('../plot/test_plots.pdf', bbox_inches='tight')
