import numpy as np

dux = 1.0
duy = dux
uxDomain = [-dux, 0.0, dux]

numsteps = 10
T = 0.55

x0 = 5.0
y0 = 5.0
vx0 = 10.0
vy0 = 0.0
print(T*0.5)
x = [x0]
vx = [vx0]

for ux in uxDomain:
    for k in range(numsteps + 1):
        vx.append(vx[k] + ux*T)
        x.append(x[k] + vx[k]*T + 0.5*ux*pow(T, 2))
        print(x)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~")
    x = [x0]
    vx = [vx0]

