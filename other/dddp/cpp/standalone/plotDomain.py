import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
import sys
import glob

font = {'family' : 'times',
        'weight' : 'bold',
        'size'   : 28}

plt.rc('font', **font)
sys.path.append("outputs/domain")
sys.path.append("outputs/solution")

# scan files
domainPath = "outputs/domain/"
domainFiles = glob.glob(domainPath + "*.py")
for f in range(len(domainFiles)):
    domainFiles[f] = domainFiles[f].replace(domainPath,"")
    domainFiles[f] = domainFiles[f].replace(".py","")
domainFiles.sort()

solPath = "outputs/solution/"
solFiles = glob.glob(solPath + "*.py")
for f in range(len(solFiles)):
    solFiles[f] = solFiles[f].replace(solPath,"")
    solFiles[f] = solFiles[f].replace(".py","")
solFiles.sort()


for f in range(len(domainFiles)):
    d = __import__(domainFiles[f])
    s = __import__(solFiles[f])

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

    idx = []
    figure(figsize=(22, 12), dpi=80)
    plt.subplot(3,2,1)
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(x[k])))
        plt.scatter(idx[k], x[k], c = 'black')
    plt.grid()
    plt.plot(xmin, linestyle = '--', c = 'red')
    plt.gca().set_prop_cycle(None)
    plt.plot(xmax, linestyle = '--', c = 'red')
    plt.plot(x0)
    plt.plot(opt_x)
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Position (m)")

    idx = []
    plt.subplot(3,2,2)
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(y[k])))
        plt.scatter(idx[k], y[k], c = 'black')
    plt.grid()
    plt.plot(ymin, linestyle = '--', c = 'red')
    plt.gca().set_prop_cycle(None)
    plt.plot(ymax, linestyle = '--', c = 'red')
    plt.plot(y0)
    plt.plot(opt_y)
    ax = plt.gca()
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Position (m)")

    idx = []
    plt.subplot(3,2,3)
    for k in range(d.domain['k']):
        idx.append(k*np.ones(len(vx[k])))
        plt.scatter(idx[k], vx[k], c = 'black')
    plt.grid()
    plt.plot(vxmin, linestyle = '--', c = 'red')
    plt.gca().set_prop_cycle(None)
    plt.plot(vxmax, linestyle = '--', c = 'red')
    plt.plot(vx0)
    plt.plot(opt_vx)
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Speed (m/s)")

    plt.subplot(3,2,5)
    plt.grid()
    plt.plot(opt_ux)
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Acceleration (m/s/s)")

    plt.subplot(3,2,6)
    plt.grid()
    plt.plot(opt_uy)
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Speed (m/s)")

    plt.tight_layout()
    plt.suptitle('Feasible domain for iteration %i with Δux=Δvx=%.2f | Δuy=Δvy=%.2f | Δx=%.2f | Δy=%.2f'
        % (d.domain['iter'], d.domain['dux'], d.domain['duy'], d.domain['dx'], d.domain['dy']), y=1.0)

    plt.show()
    # plt.savefig("test.png", bbox_inches = 'tight', pad_inches = 0.1)
    # plt.savefig("../outputs/test1.png")


