from bash_command import bash_command as bc
import numpy as np
import sys
import argparse
# import argparse

parser=argparse.ArgumentParser(description='Compute torque and precession rate on test orbit from END')
parser.add_argument('--etest',  default=0.7,
	help='eccentricity of test orbit', type=float)
parser.add_argument('--atest', default=0.99,
	help='sma of test orbit', type=float)
parser.add_argument('--ein', default=0.7,
	help='eccentricity at inner edge of END', type=float)
parser.add_argument('-q', default=0.0,
	help='specifies eccentricity gradient in disk', type=float)
parser.add_argument('-o','--pomega', default=0,
	help='specifies orientation ', type=float)
parser.add_argument('-b', default=1e-4,
	help='softening length', type=float)

args=parser.parse_args()

# e1=sys.argv[1]
# a1=sys.argv[2]
# ang=sys.argv[3]

dd=1.0
i=0
m=2.5e-7
No=1001
while (dd>0.05) and (i<11):
	# bc.bash_command('/projects/alge9397/code/c/torque2_touma/torque_integ --etest {0} --atest {1} -o {2} -n {3} -f {4} -b {5} -q {6} --ein {7}'\
	# 	.format(args.etest, args.atest, args.pomega, No, 0, args.b, args.q, args.ein))
	# sys.stdout.flush()
	# bc.bash_command('/projects/alge9397/code/c/torque2_touma/torque_integ --etest {0} --atest {1} -o {2} -n {3} -f {4} -b {5} -q {6} --ein {7}'\
	# 	.format(args.etest, args.atest, args.pomega, No, 1, args.b, args.q, args.ein))
	bc.bash_command('/projects/alge9397/code/c/torque2_touma/torque_integ --etest {0} --atest {1} -o {2} -n {3} -f {4} -b {5} -q {6} --ein {7}'\
		.format(args.etest, args.atest, args.pomega, No, 0, args.b, args.q, args.ein))
	sys.stdout.flush()
	bc.bash_command('/projects/alge9397/code/c/torque2_touma/torque_integ --etest {0} --atest {1} -o {2} -n {3} -f {4} -b {5} -q {6} --ein {7}'\
		.format(args.etest, args.atest, args.pomega, No, 1, args.b, args.q, args.ein))
	sys.stdout.flush()


	tdot1=np.genfromtxt('tau_a_e{0}_a{1}_o{2}_q{3}_ein{4}'.format(args.etest, args.atest, args.pomega, args.q, args.ein))
	tdot2=np.genfromtxt('tau_a_e{0}_a{1}_o{2}_q{3}_ein{4}'.format(args.etest, args.atest, args.pomega, args.q, args.ein))
	dd=abs((tdot1-tdot2)/tdot1)
	idot1=np.genfromtxt('i_a_e{0}_a{1}_o{2}_q{3}_ein{4}'.format(args.etest, args.atest, args.pomega, args.q, args.ein))
	idot2=np.genfromtxt('i_b_e{0}_a{1}_o{2}_q{3}_ein{4}'.format(args.etest, args.atest, args.pomega, args.q, args.ein))
	dd=abs((idot1-idot2)/idot1)

	No=No*2
	i+=1
	print args.etest, args.atest, args.pomega, i-1,No/2, abs((tdot1-tdot2)/tdot1), abs((idot1-idot2)/idot1)
	sys.stdout.flush()
