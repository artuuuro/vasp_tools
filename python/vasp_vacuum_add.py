#!/usr/bin/python

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

a = 0.0
b = 0.0
c = 0.0

if (len(sys.argv) < 4):
	print "Too few arguments"
	sys.exit(1)

try:
	a=float(sys.argv[1])
except ValueError:
	print "Vacuum size is not a number"
	sys.exit(1)

try:
	b=float(sys.argv[2])
except ValueError:
	print "Vacuum size is not a number"
	sys.exit(1)

try:
	c=float(sys.argv[3])
except ValueError:
	print "Vacuum size is not a number"
	sys.exit(1)

for line in sys.stdin:
	mass.append(line)
	counter = counter + 1

print mass[0],
print mass[1],

tmp=mass[2].split()

a1=float(tmp[0])
a2=float(tmp[1])
a3=float(tmp[2])

alength=math.sqrt(a1**2 + a2**2 + a3**2)
amul=(alength+a)/alength

tmp=mass[3].split()

b1=float(tmp[0])
b2=float(tmp[1])
b3=float(tmp[2])

blength=math.sqrt(b1**2 + b2**2 + b3**2)
bmul=(blength+b)/blength

tmp=mass[4].split()

c1=float(tmp[0])
c2=float(tmp[1])
c3=float(tmp[2])

clength=math.sqrt(c1**2 + c2**2 + c3**2)
cmul=(clength+c)/clength

#multiply vector a by factor a
print str(a1*amul) + "   " + str(a2*amul) + "   " + str(a3*amul) + "   !vector_a"
print str(b1*bmul) + "   " + str(b2*bmul) + "   " + str(b3*bmul) + "   !vector_a"
print str(c1*cmul) + "   " + str(c2*cmul) + "   " + str(c3*cmul) + "   !vector_a"

d=int(mass[5].split()[0])

print "   " + str(d) + "   !" + str(d) + " atoms in cell"
print mass[6],

atoms1=[]
atoms2=[]
atoms3=[]

for i in range(0,d,1):
	tmp = mass[7+i].split()
	atoms1.append(float(tmp[0]))
	atoms2.append(float(tmp[1]))
	atoms3.append(float(tmp[2]))

for i in range(0,d,1):
	print str(atoms1[i]/amul) + "   " + str(atoms2[i]/bmul) + "   " + str(atoms3[i]/cmul)

print 

