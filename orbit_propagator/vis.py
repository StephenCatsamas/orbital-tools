import matplotlib.pyplot as plt
import numpy as np
from math import *
import csv

t = []
x = []
y = []
z = []
vx = []
vy = []
vz = []



with open('tmp/met.dat', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for i,row in enumerate(reader):
        if i == 1:
           body_name = row[0]
           mass = float(row[1])
           radius = float(row[2])
           landing_altitude = float(row[3])
           rotational_period = float(row[4])
           wx = float(row[5])
           wy = float(row[6])
           wz = float(row[7])


with open('tmp/vis.dat', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        t.append(float(row[0]))
        x.append(float(row[1]))
        y.append(float(row[2]))
        z.append(float(row[3]))
        vx.append(float(row[4]))
        vy.append(float(row[5]))
        vz.append(float(row[6]))

q = list(zip(x,y,z))
dq = list(zip(vx,vy,vz))

r = [np.linalg.norm(p) for p in q]
vr = [np.dot(p,dp)/np.linalg.norm(p) for p,dp in zip(q,dq)]

vvr = np.subtract( dq[0] , np.multiply(q[0],np.dot(q[0],dq[0])/np.dot(q[0],q[0])))


normal_vec = np.cross(q[0],vvr)

e1 = q[0]/np.linalg.norm(q[0])
e2 = vvr/np.linalg.norm(vvr)
e3 = normal_vec/np.linalg.norm(normal_vec)

norm_basis = np.array([e1,e2,e3])

q_nb = [np.matmul(norm_basis, p) for p in q]

MA = [atan2(y,x) for x,y,z in q_nb]


fig = plt.figure(figsize=(10,8)) 
ax = fig.add_subplot(2, 2, 1, projection='polar') 
ax2 = fig.add_subplot(2, 2, 2) 
ax3 = fig.add_subplot(2, 2, 3) 
ax4 = fig.add_subplot(2, 2, 4) 

# fs = [f - state.vfs*t for f,t in zip(path[:,1],path[:,4])]

ax.plot(MA,r)
ax.set_xlabel('Downrange Angle $\phi$ (rad)');
ax.set_ylabel('Radius $r$ (m)');
pi2 = [2*pi*x/1000 for x in range(1000)]
rh = [radius+landing_altitude for x in pi2]
ax.plot(pi2, rh, color='r')
ax2.plot(r,vr)
ax2.set_xlabel('Radius $r$ (m)');
ax2.set_ylabel('Radial Velocity $v_r$ (m/s)');
ax2.axvline(radius+landing_altitude, color='r')
ax2.axhline(0, color='r')
ax2.invert_xaxis()
ax2.invert_yaxis()
ax3.plot(t,r)
ax3.set_xlabel('Time $t$ (s)');
ax3.set_ylabel('Radius $r$ (m)');
ax3.axhline(radius+landing_altitude, color='r')



surf_vel= [(-wz*p[1],wz*p[0],0) for p in q]

rel_surfvel = [np.linalg.norm(np.subtract(dp, vs)) for dp,vs in zip(dq,surf_vel)]

ax4.plot(r,rel_surfvel)
ax4.set_xlabel('Radius $r$ (m)');
ax4.set_ylabel('Surface Velocity $v_s$ (m/s)');
ax4.axvline(radius+landing_altitude, color='r')
ax4.axhline(0, color='r')
ax4.invert_xaxis()
plt.show() 
