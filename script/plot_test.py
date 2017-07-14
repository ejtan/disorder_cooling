#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil


src_dir = "../plot/data/"

if not os.path.isdir("../plot"):
    os.mkdir("../plot")
    os.mkdir("../plot/data")

for filename in os.listdir("../"):
    if filename.endswith(".txt"):
        shutil.move("../" + filename, src_dir)

T, S_clean = np.loadtxt(src_dir + "2D_ising_clean.txt", unpack=True)
S_disorder = np.loadtxt(src_dir + "2D_ising_disorder.txt", usecols=1)

plt.plot(T, S_clean)
plt.plot(T, S_disorder)

plt.ylim([0, 1])
plt.xlim([0, 5])
plt.xlabel(r"T [J / $k_B$]")
plt.ylabel(r"S [$k_B$]")
plt.title("Entropy of 2D Ising system")
plt.legend(["Clean system", r"Disorder system with $\Delta = 5$"])
plt.show()
