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
sys.path.append("../outputs/dddp/domain")
sys.path.append("../outputs/dddp/solution")

# scan files
domainPath = "../outputs/dddp/domain/"
domainFiles = glob.glob(domainPath + "*.py")
for f in range(len(domainFiles)):
    domainFiles[f] = domainFiles[f].replace(domainPath,"")
    domainFiles[f] = domainFiles[f].replace(".py","")
# sort the list
domainFiles.sort()

solPath = "../outputs/dddp/solution/"
solFiles = glob.glob(solPath + "*.py")
for f in range(len(solFiles)):
    solFiles[f] = solFiles[f].replace(solPath,"")
    solFiles[f] = solFiles[f].replace(".py","")
solFiles.sort()


for f in range(len(domainFiles)):

    # domain data
    d = __import__(domainFiles[f])
    x = d.domain['x']
    xmin = np.transpose(d.domain['xmin'])
    xmax = np.transpose(d.domain['xmax'])
    x0 = np.transpose(d.domain['x0'])
    y = d.domain['y']
    ymin = np.transpose(d.domain['ymin'])
    ymax = np.transpose(d.domain['ymax'])
    y0 = np.transpose(d.domain['y0'])
    vx = d.domain['vx']
    vxmin = np.transpose(d.domain['vxmin'])
    vxmax = np.transpose(d.domain['vxmax'])
    vx0 = np.transpose(d.domain['vx0'])

    
    ux = d.domain['ux']
    uy = d.domain['uy']
    ux0 = np.transpose(d.domain['ux0'])
    uy0 = np.transpose(d.domain['uy0'])

    # solution data
    if (len(solFiles)):
        s = __import__(solFiles[f])
        opt_x = np.transpose(s.solution['x'])
        opt_y = np.transpose(s.solution['y'])
        opt_vx = np.transpose(s.solution['vx'])
        opt_ux = np.transpose(s.solution['ux'])
        opt_uy = np.transpose(s.solution['uy'])

    # x data
    idx = []
    figure(figsize=(22, 12), dpi=80)
    plt.plot()
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(x[k])))
        plt.scatter(idx[k], x[k], c = 'black')
    plt.plot(xmin, linestyle = '--', c = 'red')
    ax.set_prop_cycle(None)
    plt.plot(xmax, linestyle = '--', c = 'red')
    plt.plot(x0, label ='Initial')
    if (len(solFiles)):
        plt.plot(opt_x, label ='Optimal')
    plt.grid(False)
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Position (m)")
    plt.legend(fontsize="18", loc ="upper left") 
    plt.show()
    
    # y data
    idx = []
    figure(figsize=(22, 12), dpi=80)
    plt.subplot(2,2,1)
    # plt.plot()
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(y[k])))
        plt.scatter(idx[k], y[k], c = 'black')
    plt.plot(ymin, linestyle = '--', c = 'red')
    ax.set_prop_cycle(None)
    plt.plot(ymax, linestyle = '--', c = 'red')
    plt.plot(y0)
    if (len(solFiles)):
        plt.plot(opt_y)
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([0, 10])
    plt.xlabel("Time Steps")
    plt.ylabel("Lat. Position (m)")

    # vx data
    idx = []
    plt.subplot(2,2,2)
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(vx[k])))
        plt.scatter(idx[k], vx[k], c = 'black')
    plt.plot(vxmin, linestyle = '--', c = 'red')
    ax.set_prop_cycle(None)
    plt.plot(vxmax, linestyle = '--', c = 'red')
    plt.plot(vx0)
    if (len(solFiles)):
        plt.plot(opt_vx)
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([20, 35])
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Speed")

    # ux data
    idx = []
    # if (f > 0):
    plt.subplot(2,2,3)
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(ux[k])))
        plt.scatter(idx[k], ux[k], c = 'black')
    plt.plot(ux0)
    if (len(solFiles)):
        plt.plot(opt_ux, c = 'orange')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-2.5, 1.5])
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Accel. (m/s/s)")

    # uy data
    idx = []
    # if (f > 0):
    plt.subplot(2,2,4)
    ax = plt.gca()
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(uy[k])))
        plt.scatter(idx[k], uy[k], c = 'black')
    plt.plot(uy0)
    if (len(solFiles)):
        plt.plot(opt_uy, c = 'orange')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-1.5, 1.5])
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral. Speed. (m/s)")

    plt.tight_layout()
    plt.suptitle('Feasible domain for iteration %i \n Δux=Δvx=%.2f | Δuy=Δvy=%.2f | Δx=%.2f | Δy=%.2f'
        % (d.domain['iter'], d.domain['dux'], d.domain['duy'], d.domain['dx'], d.domain['dy']), y=1.0)

    plt.show()
    # plt.savefig("test.png", bbox_inches = 'tight', pad_inches = 0.1)
    # plt.savefig("../outputs/test1.png")


