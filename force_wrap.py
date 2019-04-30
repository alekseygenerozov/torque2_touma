from __future__ import print_function
import ctypes 
import numpy as np
import rebound

import sys
import os


pre=os.path.join(os.path.dirname(__file__))
_force=ctypes.CDLL(pre+'/torque_integ.so')
_force.force.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
_force.force.restype = ctypes.POINTER(ctypes.c_double)

_force.force_tot_point.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,\
 ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
_force.force_tot_point.restype = ctypes.POINTER(ctypes.c_double)


def force(a, e, r1x, r1y, r1z, b):
	f1=_force.force(ctypes.c_double(a), ctypes.c_double(e), ctypes.c_double(r1x), ctypes.c_double(r1y), ctypes.c_double(r1z), ctypes.c_double(b))

	return np.array([f1[0], f1[1], f1[2]])

# def force_tot_point(a1, a2, p, e0, q, r1x, r1y, r1z, b):
def force_tot_point(params, pos):
	f1=_force.force_tot_point(ctypes.c_double(params['a1']), ctypes.c_double(params['a2']), ctypes.c_double(params['p']), ctypes.c_double(params['e0']), ctypes.c_double(params['q']),\
		ctypes.c_double(pos[0]), ctypes.c_double(pos[1]), ctypes.c_double(pos[2]), ctypes.c_double(params['b']))

	return np.array([f1[0], f1[1], f1[2]])

def leapfrog(params, pos, vel, delta_t, m_disk):
	a=m_disk*force_tot_point(params, pos)-pos/(pos[0]**2.+pos[1]**2.+pos[2]**2.)**1.5
	# a=-pos/(pos[0]**2.+pos[1]**2.+pos[2]**2.)**1.5
	vel=vel+0.5*delta_t*a
	pos=pos+delta_t*vel
	a=m_disk*force_tot_point(params, pos)-pos/(pos[0]**2.+pos[1]**2.+pos[2]**2.)**1.5
	# a=-pos/(pos[0]**2.+pos[1]**2.+pos[2]**2.)**1.5
	vel=vel+0.5*delta_t*a

	return pos,vel

# def yo8(params, pos, vel, delta_t, m_disk):
# 	d = [0.104242620869991e1, 0.182020630970714e1, 0.157739928123617e0, 0.244002732616735e1, -0.716989419708120e-2, -0.244699182370524e1, -0.161582374150097e1, -0.17808286265894516e1]
# 	for ii in range(6):
# 		pos, vel=leapfrog(params, pos, vel, d[ii]*delta_t, m_disk)
# 	pos,vel=leapfrog(params, pos, vel, d[7]*delta_t, m_disk)
# 	for ii in range(6):
# 		pos,vel=leapfrog(params, pos, vel, d[6-ii]*delta_t, m_disk)
# 	return pos,vel

def energy(pos, vel):
	return 0.5*(vel[0]**2.0+vel[1]**2.0+vel[2]**2.0)-1.0/(pos[0]**2.+pos[1]**2.+pos[2]**2.)**0.5

N=int(float(sys.argv[1]))
m_disk=float(sys.argv[2])
atest=float(sys.argv[3])
itest=float(sys.argv[4])
ptest=float(sys.argv[5])
qtest=-1.6
etest=0.9*(atest/1.0)**qtest

params0={'a1':1, 'a2':2, 'p':1.01, 'e0': etest, 'q':qtest, 'b':1e-4}
# print force_tot_point({'a1':1, 'a2':2, 'p':1.01, 'e0': 0.9, 'q':-1.6, 'b':1e-4}, np.array([0.1, 0.1, 0.1]))
sim=rebound.Simulation()
sim.G=1
sim.add(m=1, x=0, y=0, z=0, vx=0, vy=0, vz=0)
sim.add(m=1e-15, a=atest, e=etest, M=0.001, inc=itest, pomega=ptest, primary=sim.particles[0])
##Period of test orbit
per=2.0*np.pi*(atest**1.5)
sim.move_to_com()
pos,vel=np.array(sim.particles[1].xyz),np.array(sim.particles[1].vxyz)
pos0=np.copy(pos)
vel0=np.copy(vel)
en0=energy(pos0, vel0)


f=open('log_N{0}_m{1}_b{2}_a{3}_e{4:.2g}_i{5}_p{6}'.format(N, m_disk, params0['b'], atest, etest, itest, ptest), 'w')
f2=open('pos_N{0}_m{1}_b{2}_a{3}_e{4:.2g}_i{5}_p{6}'.format(N, m_disk, params0['b'], atest, etest, itest, ptest), 'w')
f3=open('tde_N{0}_m{1}_b{2}_a{3}_e{4:.2g}_i{5}_p{6}'.format(N, m_disk, params0['b'], atest, etest, itest, ptest), 'w')
norb=100
for i in range(norb*N):
	pos,vel=leapfrog(params0, pos, vel, per/N, m_disk)
	if ((i%1000)==0):
		f.write('{0:.10e} {1:.10e} \n'.format((energy(pos, vel)-en0)/en0, -1.0/2.0/energy(pos,vel)))
		f2.write('{0:.10e} {1:.10e} {2:.10e} {3:.10e} {4:.10e} {5:.10e}\n'.format(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))
	if np.linalg.norm([pos[0], pos[1], pos[2]])<1.0e-4:
		f3.write('{0:.10e} {1:.10e} {2:.10e} {3:.10e} {4:.10e} {5:.10e}\n'.format(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))

f.close()
f2.close()
f3.close()

sim=rebound.Simulation()
sim.G=1
sim.add(m=1, x=0, y=0, z=0, vx=0, vy=0, vz=0)
i=-1
sim.add(x=pos[0], y=pos[1], z=pos[2], vx=vel[0], vy=vel[1], vz=vel[2])
odat=sim.calculate_orbits(primary=sim.particles[0])
print('{0:.10e} {1:.10e} {2:.10e}\n'.format(odat[0].pomega, (odat[0].e-0.9)/0.9, (energy(pos, vel)-en0)/en0))