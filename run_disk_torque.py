from bash_command import bash_command as bc
import numpy as np
import sys

e1=sys.argv[1]
a1=sys.argv[2]
ang=sys.argv[3]
ang_rad=np.pi/180.0*float(ang)
dd=1.0
i=0

m=2.5e-7
No=10
while (dd>0.05) and (i<4):
	bc.bash_command('/projects/alge9397/code/c/torque2_touma/torque_integ {0} {1} {2} {3} {4}'.format(e1, a1, ang, No, 0))
	sys.stdout.flush()
	bc.bash_command('/projects/alge9397/code/c/torque2_touma/torque_integ {0} {1} {2} {3} {4}'.format(e1, a1, ang, No, 1))
	sys.stdout.flush()


	tdot1=np.genfromtxt('tau_a_{0:.1g}_{1:.1g}_{2:.1f}'.format(e1, a1, ang))
	tdot2=np.genfromtxt('tau_b_{0:.1g}_{1:.1g}_{2:.1f}'.format(e1, a1, ang))
	dd=abs((tdot1-tdot2)/tdot1)
	idot1=np.genfromtxt('i_a_{0:.1g}_{1:.1g}_{2:.1f}'.format(e1, a1, ang, idx))
	idot2=np.genfromtxt('i_b_{0:.1g}_{1:.1g}_{2:.1f}'.format(e1, a1, ang, idx))
	dd=max(dd, abs((idot1-idot2)/idot1))

	No=No*2
	i+=1
	print e1, a1, ang, i-1,No/2, abs((tdot1-tdot2)/tdot1), abs((idot1-idot2)/idot1)
	sys.stdout.flush()
