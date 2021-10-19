import matplotlib.pyplot as plt
from math import *
import csv

t = []
x = []
y = []
z = []
vx = []
vy = []
vz = []

with open('vis.dat', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        t.append(float(row[0]))
        x.append(float(row[1]))
        y.append(float(row[2]))
        z.append(float(row[3]))
        vx.append(float(row[4]))
        vy.append(float(row[5]))
        vz.append(float(row[6]))

r = [sqrt(x*x+y*y) for x,y in zip(x,y)]
vr = [(x*vx + y*vy)/sqrt(x*x+y*y) for x,y,vx,vy in zip(x,y,vx,vy)]


fig = plt.figure(figsize=(10,8)) 
ax = fig.add_subplot(2, 2, 1) 
ax2 = fig.add_subplot(2, 2, 2) 
ax3 = fig.add_subplot(2, 2, 3) 
ax4 = fig.add_subplot(2, 2, 4) 

# fs = [f - state.vfs*t for f,t in zip(path[:,1],path[:,4])]

# ax.plot(fs,path[:,0])
ax.set_xlabel('Surface Azumith $\phi$ (rad)');
ax.set_ylabel('Radius $r$ (m)');
# ax.axhline(state.rs, color='r')
ax2.plot(r,vr)
ax2.set_xlabel('Radius $r$ (m)');
ax2.set_ylabel('Radial Velocity $v_r$ (m/s)');
# ax2.axvline(state.rs, color='r')
ax2.axhline(0, color='r')
ax2.invert_xaxis()
ax2.invert_yaxis()
ax3.plot(t,r)
ax3.set_xlabel('Time $t$ (s)');
ax3.set_ylabel('Radius $r$ (m)');
# ax3.axhline(state.rs, color='r')

# vsrfs = [sqrt(vr*vr + (r*vf - state.rs* state.vfs)*(r*vf - state.rs* state.vfs)) for r,vr,vf in zip(path[:,0],path[:,2],path[:,3])]

# ax4.plot(path[:,0],vsrfs)
ax4.set_xlabel('Radius $r$ (m)');
ax4.set_ylabel('Surface Velocity $v_s$ (m/s)');
# ax4.axvline(state.rs, color='r')
ax4.axhline(0, color='r')
ax4.invert_xaxis()
plt.show() 
