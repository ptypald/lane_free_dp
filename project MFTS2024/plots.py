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

if (len(sys.argv) == 2):
    method = sys.argv[1]
    print(method)
elif (len(sys.argv) == 3):
    method = sys.argv[1]
    init_method = sys.argv[2]
    print(method, init_method)
else:
    method = input("Please give the method: ")
    init_method = input("Please give the initial trajectory method: ")
    print(method)

path = "../outputs/" + method + "/solution/"
sys.path.append(path)
if method == "dddp":
    pathDomain = "../outputs/" + method + "/domain/"
    sys.path.append(pathDomain)

# scan files
solFiles = glob.glob(path + "*.py")
for f in range(len(solFiles)):
    solFiles[f] = solFiles[f].replace(path,"")
    solFiles[f] = solFiles[f].replace(".py","")
solFiles.sort()

if (method == "dddp"):
    domFiles = glob.glob(pathDomain + "*.py")
    for f in range(len(domFiles)):
        domFiles[f] = domFiles[f].replace(pathDomain,"")
        domFiles[f] = domFiles[f].replace(".py","")
    domFiles.sort()

xFirst = yFirst = vxFirst = vyFirst = uxFirst = uyFirst = None
for f in range(len(solFiles)):

    # solution data
    s = __import__(solFiles[f])
    opt_x = np.transpose(s.solution['x'])
    opt_y = np.transpose(s.solution['y'])
    opt_vx = np.transpose(s.solution['vx'])
    if method == "ddp":
        opt_vy = np.transpose(s.solution['vy'])
    opt_ux = np.transpose(s.solution['ux'])
    opt_uy = np.transpose(s.solution['uy'])

    x0 = np.transpose(s.solution['x0'])
    y0 = np.transpose(s.solution['y0'])
    vx0 = np.transpose(s.solution['vx0'])
    if method == "ddp":
        vy0 = np.transpose(s.solution['vy0'])
    ux0 = np.transpose(s.solution['ux0'])
    uy0 = np.transpose(s.solution['uy0'])

    if (method == "dddp"):
        # domain data
        d = __import__(domFiles[f])
        x = d.domain['x']
        xmin = np.transpose(d.domain['xmin'])
        xmax = np.transpose(d.domain['xmax'])
        y = d.domain['y']
        ymin = np.transpose(d.domain['ymin'])
        ymax = np.transpose(d.domain['ymax'])
        vx = d.domain['vx']
        vxmin = np.transpose(d.domain['vxmin'])
        vxmax = np.transpose(d.domain['vxmax'])       
        ux = d.domain['ux']
        uy = d.domain['uy']
       
    if xFirst is None:
        xFirst = x0
        yFirst = y0
        vxFirst = vx0
        vyFirst = vy0
        uxFirst = ux0
        uyFirst = uy0
        
    # x data
    figure(figsize=(22, 12), dpi=80)
    plt.subplot(3,2,1)
    ax = plt.gca()
    plt.plot(opt_x, label ='Optimal', c = 'orange')
    plt.plot(x0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Position \n(in m)")
    plt.legend(fontsize="18", loc ="upper left")

    # y data
    # figure(figsize=(22, 12), dpi=80)
    plt.subplot(3,2,2)
    # plt.plot()
    ax = plt.gca()
    plt.plot(opt_y, c = 'orange')
    plt.plot(y0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([0, 11])
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Position \n(in m)")

    # vx data
    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(opt_vx, c = 'orange')
    plt.plot(vx0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([20, 35])
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Speed \n(in m/s)")

    # vy data
    plt.subplot(3,2,4)
    ax = plt.gca()
    if method == "ddp":
        plt.plot(opt_vy, c = 'orange')
        plt.plot(vy0, label ='Initial', c = 'blue', linestyle = '--')
    else:
        plt.plot(opt_uy, c = 'orange')
        plt.plot(uy0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-2, 2])
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Speed \n(in m/s)")

    # ux data
    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(opt_ux, c = 'orange')
    plt.plot(ux0, label ='Initial', c = 'blue', linestyle = '--')
    plt.plot(np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
    plt.plot(-2.0*np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-2.5, 2.0])
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Acceleration \n(in m/s/s)")

    # uy data
    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(opt_uy, c = 'orange')
    plt.plot(uy0, label ='Initial', c = 'blue', linestyle = '--')
    plt.plot(np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
    plt.plot(-np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
    plt.grid(False)
    # ax.set_xlim()
    ax.set_ylim([-1.5, 1.5])
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Acceleration \n(in m/s/s)")

    plt.tight_layout()
    # plt.show()
    
    outFile = "../outputs/" + method + "/" + method + "_" + str(f+1) + "_" + init_method + ".png"
    plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)



# # plt.show()
# outFile = "../outputs/" + method + "/" + method + "_" + "last_first" + "_" + init_method + ".png"
# plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1)
    


