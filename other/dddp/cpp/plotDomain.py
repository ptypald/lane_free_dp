import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
import sys

sys.path.append("outputs")
from domain_0 import *

font = {'family' : 'times',
        'weight' : 'bold',
        'size'   : 28}

plt.rc('font', **font)

x = domain['x']
xmin = np.transpose(domain['xmin'])
xmax = np.transpose(domain['xmax'])
x0 = np.transpose(domain['x0'])
y = domain['y']
ymin = np.transpose(domain['ymin'])
ymax = np.transpose(domain['ymax'])
y0 = np.transpose(domain['y0'])
vx = domain['vx']
vxmin = np.transpose(domain['vxmin'])
vxmax = np.transpose(domain['vxmax'])
vx0 = np.transpose(domain['vx0'])

idx = []
figure(figsize=(22, 12), dpi=80)
plt.subplot(2,2,1)
for k in range(domain['k']):
    idx.append(k*np.ones(len(x[k])))
    plt.scatter(idx[k], x[k], c = 'black')
plt.grid()
plt.plot(xmin, linestyle = '--', c = 'red')
plt.gca().set_prop_cycle(None)
plt.plot(xmax, linestyle = '--', c = 'red')
plt.plot(x0)
plt.xlabel("Time Steps")
plt.ylabel("Long. Position (m)")

idx = []
plt.subplot(2,2,2)
for k in range(domain['k']):
    idx.append(k*np.ones(len(y[k])))
    plt.scatter(idx[k], y[k], c = 'black')
plt.grid()
plt.plot(ymin, linestyle = '--', c = 'red')
plt.gca().set_prop_cycle(None)
plt.plot(ymax, linestyle = '--', c = 'red')
plt.plot(y0)
ax = plt.gca()
plt.xlabel("Time Steps")
plt.ylabel("Lateral Position (m)")

idx = []
plt.subplot(2,2,3)
for k in range(domain['k']):
    idx.append(k*np.ones(len(vx[k])))
    plt.scatter(idx[k], vx[k], c = 'black')
plt.grid()
plt.plot(vxmin, linestyle = '--', c = 'red')
plt.gca().set_prop_cycle(None)
plt.plot(vxmax, linestyle = '--', c = 'red')
plt.plot(vx0)
plt.xlabel("Time Steps")
plt.ylabel("Long. Speed (m/s)")

plt.tight_layout()
plt.suptitle('Feasible domain for iteration %i with Δux=Δvx=%.2f | Δuy=Δvy=%.2f | Δx=%.2f | Δy=%.2f'
    % (domain['iter'], domain['dux'], domain['duy'], domain['dx'], domain['dy']), y=1.0)

# plt.show()
# plt.savefig("test.png", bbox_inches = 'tight', pad_inches = 0.1)
plt.savefig("outputs/test1.png")

