import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

dof, L2error, H1error = np.loadtxt("../data/convergence.dat", unpack=True)

plt.figure(1)

plt.loglog(dof, L2error, 'b-', label="L2 error")
plt.loglog(dof, H1error, 'g-', label="H1 error")
plt.grid(True)
plt.xlabel("DoF")
plt.ylabel("error")
plt.title("Convergence Plot")
plt.legend()

plt.show()