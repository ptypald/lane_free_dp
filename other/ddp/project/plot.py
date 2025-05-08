import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
import sys
import glob
import re

font = {'family' : 'times',
        'weight' : 'bold',
        'size'   : 28}

plt.rc('font', **font)
# plt.rcParams['text.usetex'] = True
sys.path.append("outputs/domain")
sys.path.append("outputs/solution")

# scan files
domainPath = "outputs/domain/"
domainFiles = glob.glob(domainPath + "*.py")
for f in range(len(domainFiles)):
    domainFiles[f] = domainFiles[f].replace(domainPath,"")
    domainFiles[f] = domainFiles[f].replace(".py","")
# sort the list
domainFiles.sort()

solPath = "outputs/solution/"
solFiles = glob.glob(solPath + "*.py")
for f in range(len(solFiles)):
    solFiles[f] = solFiles[f].replace(solPath,"")
    solFiles[f] = solFiles[f].replace(".py","")
solFiles.sort()


for f in range(len(domainFiles)):
# for f in range(1):
    d = __import__(domainFiles[f])
    s = __import__(solFiles[f])

    # set data
    x = d.domain['x']
    opt_x = np.transpose(s.solution['x'])
    xmin = np.transpose(d.domain['xmin'])
    xmax = np.transpose(d.domain['xmax'])
    x0 = np.transpose(d.domain['x0'])
    y = d.domain['y']
    opt_y = np.transpose(s.solution['y'])
    ymin = np.transpose(d.domain['ymin'])
    ymax = np.transpose(d.domain['ymax'])
    y0 = np.transpose(d.domain['y0'])
    vx = d.domain['vx']
    opt_vx = np.transpose(s.solution['vx'])
    vxmin = np.transpose(d.domain['vxmin'])
    vxmax = np.transpose(d.domain['vxmax'])
    vx0 = np.transpose(d.domain['vx0'])

    opt_ux = np.transpose(s.solution['ux'])
    opt_uy = np.transpose(s.solution['uy'])

    # x data
    idx = []
    figure(figsize=(22, 12), dpi=80)
    plt.subplot(3,2,1)
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(x[k])))
        plt.scatter(idx[k], x[k], c = 'black')
    plt.plot(xmin, linestyle = '--', c = 'red')
    ax.set_prop_cycle(None)
    plt.plot(xmax, linestyle = '--', c = 'red')
    plt.plot(x0, label ='Initial')
    plt.plot(opt_x, label ='Optimal')
    plt.grid(False)
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Position (m)")
    plt.legend(fontsize="18", loc ="upper left") 

    # y data
    idx = []
    plt.subplot(3,2,2)
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(y[k])))
        plt.scatter(idx[k], y[k], c = 'black')
    plt.plot(ymin, linestyle = '--', c = 'red')
    ax.set_prop_cycle(None)
    plt.plot(ymax, linestyle = '--', c = 'red')
    plt.plot(y0)
    plt.plot(opt_y)
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([0, 10])
    plt.xlabel("Time Steps")
    plt.ylabel("Lat. Position (m)")

    # vx data
    idx = []
    plt.subplot(3,2,3)
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(vx[k])))
        plt.scatter(idx[k], vx[k], c = 'black')
    plt.plot(vxmin, linestyle = '--', c = 'red')
    ax.set_prop_cycle(None)
    plt.plot(vxmax, linestyle = '--', c = 'red')
    plt.plot(vx0)
    plt.plot(opt_vx)
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([10, 30])
    plt.xlabel("Time Steps")
    plt.ylabel("Lat. Speed")

    # ux data
    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(opt_ux, c = 'orange')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-2, 1])
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Accel. (m/s/s)")

    # uy data
    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(opt_uy, c = 'orange')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-1, 1])
    plt.xlabel("Time Steps")
    plt.ylabel("Lat. Speed (m/s)")

    plt.tight_layout()
    plt.suptitle('Feasible domain for iteration %i with Δux=Δvx=%.2f | Δuy=Δvy=%.2f | Δx=%.2f | Δy=%.2f'
        % (d.domain['iter'], d.domain['dux'], d.domain['duy'], d.domain['dx'], d.domain['dy']), y=1.0)

    plt.show()
    # plt.savefig("test.png", bbox_inches = 'tight', pad_inches = 0.1)
    # plt.savefig("../outputs/test1.png")


