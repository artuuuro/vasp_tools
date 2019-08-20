#!/usr/bin/env python3

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

counter = 0
mass = []

if (len(sys.argv) < 2):
  print("Too few arguments")
  sys.exit(1)

try:
  a=float(sys.argv[1])
except ValueError:
  print("Multiplier is not a number")
  sys.exit(1)

if(a<0):
  print("Multiplier is negative")
  sys.exit(1)

# read POSCAR
for line in sys.stdin:
  mass.append(line)
  counter = counter + 1

print(mass[0], end='')
scale=float(mass[1])
print('{0:.9f}'.format(scale))

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
# end read POSCAR

## deformation matrix
d11 = 1.0 + a
d12 = 0.0
d13 = 0.0

d21 = 0.0
d22 = 1.0
d23 = 0.0

d31 = 0.0
d32 = 0.0
d33 = 1.0
## end deformation matrix

a1 = a1 * d11 + b1 * d12 + c1 * d13

#multiply vector a by factor a
print('{0:.9f}  {1:.9f}  {2:.9f}'.format(a1, a2, a3))
print('{0:.9f}  {1:.9f}  {2:.9f}'.format(b1, b2, b3))
print('{0:.9f}  {1:.9f}  {2:.9f}'.format(c1, c2, c3))

## this does not change
print(mass[5], end='')
print(mass[6], end='')
print(mass[7], end='')

## TODO this works only with Direct method probably

oldd=[int(x) for x in mass[6].split()]
for i in range(0,sum(oldd),1):
  [atoms1, atoms2, atoms3] = mass[8+i].split()
  print('{0:.9f}  {1:.9f}  {2:.9f}'.format(float(atoms1), float(atoms2), float(atoms3)))

print("") 

