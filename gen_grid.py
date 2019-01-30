import numpy as np

a1s=np.arange(0., 1.01, 0.1)
a1s[0]=0.01
angs=np.arange(0.0, 91., 5.)
e1s=[0.7]
# e1s=np.arange(0.0, 1.01, 0.02)
# e1s[-1]=0.99
# e1s[0]=0.01
# offsets=[0,1]

for e1 in e1s:
	for a1 in a1s:
		for ang in angs:
			temp=open('template.sh').read()
			temp=temp.replace('xx', '{0:.2f}'.format(e1))
			temp=temp.replace('yy', '{0:.2f}'.format(a1))
			temp=temp.replace('zz', '{0:.0f}'.format(ang))
			f=open('e{0:.2f}_a{1:.2f}_ang{2:.0f}.sh'.format(e1, a1, ang), 'w')
			f.write(temp)
			f.close()