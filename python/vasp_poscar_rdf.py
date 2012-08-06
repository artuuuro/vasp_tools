#!/usr/bin/python2

# Author: Artur Tamm arturt@ut.ee

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with 
# this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import math

counter = 0
mass = []

r0 = 0.0
r1 = 10.0
bins = 100

mirrors = 1

if (len(sys.argv) > 2):
	try:
		r1 = float(sys.argv[1])
		bins = int(sys.argv[2])
	except ValueError:
		print "Some problem with input numbers"
		sys.exit(1)

dr = (r1-r0)/bins

for line in sys.stdin:
	mass.append(line)
	counter = counter + 1

scale = float(mass[1])

tmp=mass[2].split()

a1=float(tmp[0])
a2=float(tmp[1])
a3=float(tmp[2])

tmp=mass[3].split()

b1=float(tmp[0])
b2=float(tmp[1])
b3=float(tmp[2])

tmp=mass[4].split()

c1=float(tmp[0])
c2=float(tmp[1])
c3=float(tmp[2])

atomnr=[int(x) for x in mass[6].split()]

atoms1=[]
atoms2=[]
atoms3=[]
atoms4=[]

x=[]
y=[]
z=[]
q=[]

j = 0
k = 0

for i in xrange(0,sum(atomnr),1):
	tmp = mass[8+i].split()
	atoms1.append(float(tmp[0]))
	atoms2.append(float(tmp[1]))
	atoms3.append(float(tmp[2]))
	atoms4.append(j)
	k+=1
	if(k == atomnr[j]):
		j+=1
		k=0

for n in xrange(0,sum(atomnr),1):
	for i in xrange(-mirrors,mirrors+1,1):
		for j in xrange(-mirrors,mirrors+1,1):
			for k in xrange(-mirrors,mirrors+1,1):
				x.append((atoms1[n]*a1+i*a1+atoms2[n]*b1+j*b1+atoms3[n]*c1+c1*k)*scale)
				y.append((atoms1[n]*a2+i*a2+atoms2[n]*b2+j*b2+atoms3[n]*c2+c2*k)*scale)
				z.append((atoms1[n]*a3+i*a3+atoms2[n]*b3+j*b3+atoms3[n]*c3+c3*k)*scale)
				q.append(atoms4[n])

#now do the radial distribution thing
hist=[[0]*(bins+1)]*(len(atomnr)+1)

for n in xrange(0, sum(atomnr),1):
	x0 = (atoms1[n]*a1 + atoms2[n]*b1 + atoms3[n]*c1)*scale
	y0 = (atoms1[n]*a2 + atoms2[n]*b2 + atoms3[n]*c2)*scale
	z0 = (atoms1[n]*a3 + atoms2[n]*b3 + atoms3[n]*c3)*scale
	i0 = atoms4[n]
	
	for i in xrange(0,len(x),1):
		j = int(math.sqrt((x0-x[i])**2+(y0-y[i])**2+(z0-z[i])**2)/dr)
#		print n,i,j,q[i]
		
		if(j > bins):
			continue
		elif(j <= 0):
			continue

		hist[0][j] += 1
		if(q[i] == i0):
			hist[i0+1][j] += 1

integrate=0.0

print "#distance\tnr\tintegrate\tpairwise"
for i in xrange(0,bins+1,1):
	print str(dr*i) + "\t" + str(float(hist[0][i])/float(sum(atomnr))),
	integrate += float(hist[0][i])/float(sum(atomnr))
	print "\t" + str(integrate),

	for j in xrange(1,len(atomnr)+1,1):
		print "\t" + str(float(hist[j][i])/float(sum(atomnr))),
	
	print ""

