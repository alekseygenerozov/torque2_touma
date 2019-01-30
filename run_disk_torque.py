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
	bc.bash_command('~/code/c/torque2_touma/torque1 {0} {1} {2} {3} {4}'.format(e1, a1, ang, No, 0))
	sys.stdout.flush()
	bc.bash_command('~/code/c/torque2_touma/torque1 {0} {1} {2} {3} {4}'.format(e1, a1, ang, No, 1))
	sys.stdout.flush()


	tdot1=np.genfromtxt('tau_a_{0}_{1}_{2}'.format(e1, a1, ang))
	tdot2=np.genfromtxt('tau_b_{0}_{1}_{2}'.format(e1, a1, ang))
	dd=abs((tdot1-tdot2)/tdot1)

	No=No*2
	i+=1
	print e1, a1, ang, i-1,No/2, dd
	sys.stdout.flush()
