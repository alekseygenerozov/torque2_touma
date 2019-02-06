import numpy as np
from scipy.interpolate import griddata, interpn
import sys

from scipy.integrate import ode

loc='//home/aleksey/Dropbox/projects/disk_torque/torque_data_2/grid_disk_3/'

Ntrials=1
ang_test=np.arange(0., 91, 5)
ang_test_rad=ang_test*np.pi/180.0
a_test=np.arange(0.1, 1.01,0.1)
# a_test=[1.0]
e_test=[0.5, 0.6, 0.7, 0.8, 0.9]
# e_test[-1]=0.99
# e_test[0]=0.01
m=2.5e-7
mdisk=1000.0*m
t_sec=2.0*np.pi*(1.0/mdisk)
t_norm=t_sec/(2.0*np.pi)


def j(e,a):
	'''Specific angular momentum'''
	return (a*(1-e**2.))**0.5

def eccentricity(j, a):
	'''Eccentricity'''
	return (1-j**2.0/a)**0.5

jdot=np.zeros([len(e_test), len(a_test), len(ang_test), Ntrials])
idot=np.zeros([len(e_test), len(a_test), len(ang_test), Ntrials])
idot_avg=np.zeros([len(e_test), len(a_test), len(ang_test)])
jdot_avg=np.zeros([len(e_test), len(a_test), len(ang_test)])

for ii,e1 in enumerate(e_test):
	for jj,a1 in enumerate(a_test):
		for kk,ang in enumerate(ang_test[1:]):
				jdot[ii, jj, kk+1]=np.genfromtxt(loc+'tau_N1000_a_{0:g}_{1}_{2}_{3}'.format(e1, a1, ang))/m
				idot[ii, jj, kk+1]=np.genfromtxt(loc+'i_N1000_a_{0:g}_{1}_{2}_{3}'.format(e1, a1, ang))/m

idot_avg=idot
jdot_avg=jdot

def jdot_interp(e, a, omega):
	# ee, aa, oo=np.meshgrid(e_test, a_test, ang_test_rad)
	# ee_i, aa_i, oo_i=np.meshgrid(e, a, np.abs(omega))
	sign=1.0
	if omega<0:
		sign=-1.0
	return sign*interpn((e_test, a_test, ang_test_rad), jdot_avg, [e, a, abs(omega)], method='linear')


def idot_interp(e, a, omega):
	return interpn((e_test, a_test, ang_test_rad), idot_avg, [e, a, abs(omega)])


def rhs(t,y):
	j=y[0]
	a=y[1]
	omega=y[2]

	e=eccentricity(j, a)
	return [jdot_interp(e, a, omega), 0, idot_interp(e, a, omega)]

##Initial conditions for test particle 
e_part=0.7
a_part=1.0
j_part=j(e_part, a_part)
omega_part=0.1
##Time step
t_tot=100.0*t_norm
delta_t=2.0*np.pi
# delta_t=t_norm/100.0
t=0.0

r=ode(rhs)
y0=[j_part, a_part, omega_part]
r.set_initial_value(y0, t)

while r.successful() and r.t<t_tot:
	print r.t/(2.0*np.pi), eccentricity(r.y[0], r.y[1]), r.y[2]*180.0/np.pi
	r.integrate(r.t+delta_t)
