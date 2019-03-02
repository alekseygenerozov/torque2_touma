import numpy as np

a1s=np.arange(0., 1.01, 0.1)
a1s[0]=0.01
# a1s=[1.0]
angs=np.arange(0.0, 180.1, 5.)
e1s=np.arange(0, 0.91, 0.1)
e1s[0]=0.01
e2s=[0.99, 1.0-1.0e-6]
e1s=np.concatenate([e1s, e2s])
# e1s[-1]=0.99
# e1s[0]=0.01
# offsets=[0,1]

for e1 in e1s:
	for a1 in a1s:
		for ang in angs:
			temp=open('template.sh').read()
			temp=temp.replace('xx1', '{0}'.format(e1))
			temp=temp.replace('xx', '--etest {0}'.format(e1))
			temp=temp.replace('yy1', '{0}'.format(a1))
			temp=temp.replace('yy', '--atest {0:.1g}'.format(a1))
			temp=temp.replace('zz1', '{0}'.format(ang))
			temp=temp.replace('zz', '--omega {0:.1f} -q -1.6 --ein 0.9'.format(ang))
			f=open('e{0}_a{1:.1g}_ang{2:.1f}.sh'.format(e1, a1, ang), 'w')
			f.write(temp)
			f.close()