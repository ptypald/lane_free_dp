import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
import sys
import glob
import re

font = {'family' : 'times',
        'weight' : 'bold',
        'size'   : 82}

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

path = "../outputs/" + method + "/journal/solutions/ddp/temp/"
sys.path.append(path)
print(path)
# scan files
solFiles = glob.glob(path + "*.py")
for f in range(len(solFiles)):
    solFiles[f] = solFiles[f].replace(path,"")
    solFiles[f] = solFiles[f].replace(".py","")
solFiles.sort()

for f in range(len(solFiles)):

    # solution data
    s = __import__(solFiles[f])
    opt_x = np.transpose(s.solution['x'])
    opt_y = np.transpose(s.solution['y'])
    opt_vx = np.transpose(s.solution['vx'])
    opt_vy = np.transpose(s.solution['vy'])
    opt_ux = np.transpose(s.solution['ux'])
    opt_uy = np.transpose(s.solution['uy'])

    x0 = np.transpose(s.solution['x0'])
    y0 = np.transpose(s.solution['y0'])
    vx0 = np.transpose(s.solution['vx0'])
    vy0 = np.transpose(s.solution['vy0'])
    ux0 = np.transpose(s.solution['ux0'])
    uy0 = np.transpose(s.solution['uy0'])
           
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
    plt.subplot(3,2,2)
    ax = plt.gca()
    plt.plot(opt_y, c = 'orange')
    plt.plot(y0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
    ax.set_ylim([0, 11])
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Position \n(in m)")

    # vx data
    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(opt_vx, c = 'orange')
    plt.plot(vx0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
    ax.set_ylim([20, 35])
    plt.xlabel("Time Steps")
    plt.ylabel("Long. Speed \n(in m/s)")

    # vy data
    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(opt_vy, c = 'orange')
    plt.plot(vy0, label ='Initial', c = 'blue', linestyle = '--')
    plt.grid(False)
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
    ax.set_ylim([-5.5, 5.5])
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
    ax.set_ylim([-5.5, 5.5])
    plt.xlabel("Time Steps")
    plt.ylabel("Lateral Acceleration \n(in m/s/s)")

    # plt.tight_layout()
    plt.show()
    
    # outFile = "../outputs/" + method + "/" + method + "_" + str(f+1) + "_" + init_method + ".png"
    # plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)


## seperate plots
# # y data
# for f in range(len(solFiles)):
#     s = __import__(solFiles[f])
#     opt_y = np.transpose(s.solution['y'])
#     y0 = np.transpose(s.solution['y0'])

#     figure(figsize=(22, 12), dpi=80)    
#     # plt.subplot(3,2,f+1)
#     ax = plt.gca()
#     plt.plot(opt_y, label ='Optimal', c = 'orange', linewidth=6.0)
#     plt.plot(y0, label ='Initial', c = 'blue', linestyle = '--', linewidth=6.0)
#     plt.grid(False)
#     ax.set_ylim([0, 11])
#     plt.xlabel("Time Steps")
#     plt.ylabel("Lateral Position (m)")
#     plt.legend(fontsize="52", loc ="upper left")

#     # plt.tight_layout()
#     # plt.show()
#     outFile = "../outputs/" + method + "/" + "(" + str(f+1) + ") " + "y_" + method + "_" + init_method + ".png"
#     plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)

# # vy data
# for f in range(len(solFiles)):
#     s = __import__(solFiles[f])
#     opt_vy = np.transpose(s.solution['vy'])
#     vy0 = np.transpose(s.solution['vy0'])

#     figure(figsize=(22, 12), dpi=80)    
#     ax = plt.gca()
#     plt.plot(opt_vy, label ='Optimal', c = 'orange', linewidth=6.0)
#     plt.plot(vy0, label ='Initial', c = 'blue', linestyle = '--', linewidth=6.0)
#     plt.grid(False)
#     ax.set_ylim([-2, 2])
#     plt.xlabel("Time Steps")
#     plt.ylabel("Lateral Speed (m/s)")
#     plt.legend(fontsize="52", loc ="upper right")
    
#     # plt.tight_layout()
#     # plt.show()
#     outFile = "../outputs/" + method + "/" + "(" + str(f+1) + ") " + "vy_" + method + "_" + init_method + ".png"
#     plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)

# # uy data
# for f in range(len(solFiles)):
#     s = __import__(solFiles[f])
#     opt_uy = np.transpose(s.solution['uy'])
#     uy0 = np.transpose(s.solution['uy0'])

#     figure(figsize=(22, 12), dpi=80)    
#     ax = plt.gca()
#     plt.plot(opt_uy, label ='Optimal', c = 'orange', linewidth=6.0)
#     plt.plot(uy0, label ='Initial', c = 'blue', linestyle = '--', linewidth=6.0)
#     # plt.plot(np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
#     # plt.plot(-np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
#     plt.grid(False)
#     ax.set_ylim([-1.5, 1.5])
#     plt.xlabel("Time Steps")
#     plt.ylabel("Lateral Accel. (m/s/s)")
#     plt.legend(fontsize="52", loc ="upper right")

    
#     # plt.tight_layout()
#     # plt.show()
#     outFile = "../outputs/" + method + "/" + "(" + str(f+1) + ") " + "uy_" + method + "_" + init_method + ".png"
#     plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)





# # x data
# for f in range(len(solFiles)):
#     s = __import__(solFiles[f])
#     opt_x = np.transpose(s.solution['x'])
#     x0 = np.transpose(s.solution['x0'])

#     figure(figsize=(22, 12), dpi=80)    
#     # plt.subplot(3,2,f+1)
#     ax = plt.gca()
#     plt.plot(opt_x, label ='Optimal', c = 'orange', linewidth=6.0)
#     plt.plot(x0, label ='Initial', c = 'blue', linestyle = '--', linewidth=6.0)
#     plt.grid(False)
#     # ax.set_ylim([0, 11])
#     plt.xlabel("Time Steps")
#     plt.ylabel("Long. Position (m)")
#     plt.legend(fontsize="52", loc ="upper left")

#     # plt.tight_layout()
#     # plt.show()
#     outFile = "../outputs/" + method + "/" + "(" + str(f+1) + ") " + "x_" + method + "_" + init_method + ".png"
#     plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)

# vx data
# for f in range(len(solFiles)):
#     s = __import__(solFiles[f])
#     opt_vx = np.transpose(s.solution['vx'])
#     vx0 = np.transpose(s.solution['vx0'])

#     figure(figsize=(22, 12), dpi=80)    
#     ax = plt.gca()
#     plt.plot(opt_vx, label ='Optimal', c = 'orange', linewidth=6.0)
#     plt.plot(vx0, label ='Initial', c = 'blue', linestyle = '--', linewidth=6.0)
#     plt.grid(False)
#     ax.set_ylim([15, 35])
#     plt.xlabel("Time Steps")
#     plt.ylabel("Long. Speed (m/s)")
#     plt.legend(fontsize="52", loc ="lower right")
    
#     # plt.tight_layout()
#     # plt.show()
#     outFile = "../outputs/" + method + "/" + "(" + str(f+1) + ") " + "vx_" + method + "_" + init_method + ".png"
#     plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)

# # uy data
# for f in range(len(solFiles)):
#     s = __import__(solFiles[f])
#     opt_ux = np.transpose(s.solution['ux'])
#     ux0 = np.transpose(s.solution['ux0'])

#     figure(figsize=(22, 12), dpi=80)    
#     ax = plt.gca()
#     plt.plot(opt_ux, label ='Optimal', c = 'orange', linewidth=6.0)
#     plt.plot(ux0, label ='Initial', c = 'blue', linestyle = '--', linewidth=6.0)
#     # plt.plot(np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
#     # plt.plot(-np.ones(s.solution['k']), label ='Bounds', c = 'red', linestyle = '--')
#     plt.grid(False)
#     ax.set_ylim([-2, 1.5])
#     plt.xlabel("Time Steps")
#     plt.ylabel("Long. Accel. (m/s/s)")
#     plt.legend(fontsize="52", loc ="upper right")

    
#     # plt.tight_layout()
#     # plt.show()
#     outFile = "../outputs/" + method + "/" + "(" + str(f+1) + ") " + "ux_" + method + "_" + init_method + ".png"
#     plt.savefig(outFile, bbox_inches = 'tight', pad_inches = 0.1, dpi = 300)

