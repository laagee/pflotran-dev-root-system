import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

path = []
path.append('.')

mc = np.genfromtxt('tracer_1D_MC-obs-sec-0.tec',skip_header=1)
frac = np.genfromtxt('tracer_1D_MC-obs-0.tec',skip_header=1)

f = plt.figure(figsize=(6,6))
ax=plt.subplot(1,1,1)
f.suptitle("Matrix Breakthrough Curves",fontsize=16)
plt.xlabel('Time [s]')
ax.set_xscale('log')
ax.set_yscale('log')
plt.ylabel('Concentration [mol/L]')

for i in range(1,6):
  lab = 'matrix ' + str(i)
  plt.plot(mc[:,0],mc[:,i],label=lab)

plt.plot(frac[:,0],frac[:,1],label='Fracture')
plt.xlim(1e4,1e6)
plt.legend(loc=3)
plt.show()


