import numpy as np
from scipy.interpolate import griddata, interpn
import sys

from scipy.integrate import ode
import json

loc='//home/aleksey/Dropbox/projects/disk_torque/torque_data_2/dat_q-1.6/touma_compact/'



# Ntrials=1
ang_test=np.arange(0., 181, 5)
ang_test_rad=ang_test*np.pi/180.0
# a_test=np.arange(0.1, 1.01,0.1)
# a_test=[0.5]
a0=float(sys.argv[4])
e_test=np.arange(0., 0.91, 0.1)
e_test[0]=0.01
ecrit=0.999999
e_test_2=[0.99, ecrit, ecrit]
e_test=np.concatenate([e_test, e_test_2])
# e_test[-1]=0.99
# e_test[0]=0.01
m=2.5e-7
mdisk=1000.0*m
t_sec=2.0*np.pi*(1.0/mdisk)
t_norm=t_sec/(2.0*np.pi)
idot0=float(sys.argv[3])

print "#{0}".format((idot0)/m*t_norm)
f=open(loc+"prec_rate.json", "rb")
prec_rate=json.loads(json.load(f))

def renorm(omega):
	if omega>(np.pi):
		return -(2.0*np.pi-omega)
	elif omega<(-np.pi):
		return 2.0*np.pi+omega
	else:
		return omega

def j(e,a):
	'''Specific angular momentum'''
	return (a*(1-e**2.))**0.5

def eccentricity(j, a):
	'''Eccentricity'''
	return (1-j**2.0/a)**0.5

jdot=np.zeros([len(e_test), len(ang_test)])
idot=np.zeros([len(e_test), len(ang_test)])
idot_avg=np.zeros([len(e_test), len(ang_test)])
jdot_avg=np.zeros([len(e_test), len(ang_test)])

for ii,e1 in enumerate(e_test):
	for kk,ang in enumerate(ang_test):
			jdot[ii, kk]=prec_rate['tau_a_e{0}_a{1:.1g}_o{2:.1f}_q-1.6_ein0.9'.format(e1, a0, ang)]/m
			if np.isnan(jdot[ii, kk]):
				jdot[ii, kk]=prec_rate['tau_b_e{0}_a{1:.1g}_o{2:.1f}_q-1.6_ein0.9'.format(e1, a0, ang)]/m

			idot[ii, kk]=prec_rate['i_a_e{0}_a{1:.1g}_o{2:.1f}_q-1.6_ein0.9'.format(e1, a0, ang)]/m
			if np.isnan(idot[ii, kk]):
				idot[ii, kk]=prec_rate['i_b_e{0}_a{1:.1g}_o{2:.1f}_q-1.6_ein0.9'.format(e1, a0, ang)]/m
			
			if ii==len(e_test)-1:
				idot[ii,kk]=0.0


idot_avg=idot
jdot_avg=jdot
e_test[-1]=1.0

def jdot_interp(e, a, omega, j1):
	sign=1.0
	if omega<0:
		sign=-1.0
	omega=renorm(omega)	
	return sign*interpn(( (1.0-e_test)[::-1], ang_test_rad), jdot_avg[::-1,:], [1-e, abs(omega)], method='linear')


def idot_interp(e, a, omega, j1):
	omega=renorm(omega)

	return np.sign(j1)*interpn(( (1.0-e_test)[::-1], ang_test_rad), idot_avg[::-1,:], [1-e, abs(omega)], method='linear')-idot0/m


def rhs(t,y):
	j1=y[0]
	a=y[1]
	omega=y[2]
	omega=renorm(omega)

	e=eccentricity(j1, a)
	return [jdot_interp(e, a, omega, j1), 0, idot_interp(e, a, omega, j1)]

##Initial conditions for test particle 
e_part=float(sys.argv[1])
a_part=a0
j_part=j(e_part, a_part)
omega_part=float(sys.argv[2])
##Time step
t_tot=1000.0*t_norm
delta_t=2.0*np.pi
# delta_t=t_norm/100.0
t=0.0

r=ode(rhs)
y0=[j_part, a_part, omega_part]
r.set_initial_value(y0, t)

f=open('sol_pert_e{0}_om{1}_idot{2}_a{3}'.format(e_part, omega_part, idot0, a0), 'w')
while r.successful() and r.t<t_tot:
	last_e=eccentricity(r.y[0], a0)
	f.write('{0} {1} {2} {3}\n'.format(r.t/(2.0*np.pi), eccentricity(r.y[0], a0), r.y[2]*180.0/np.pi, r.y[0]))
	r.integrate(r.t+delta_t)
	r.y[2]=renorm(r.y[2])
	if ((last_e>ecrit) and (eccentricity(r.y[0], a0)>ecrit)):
		print 'lost @ time:{0}'.format(r.t)
f.close()