#!/usr/bin/python2.7

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

# Randomisation is done by substituting atoms from original structure with new atoms
# Controllable parameters are:
# Number of substitutable atoms
# Radius where the next atom will be selected from
# Size of the displacement for the new atom
# Min distance between the atoms
# Also the name of atoms that will be substituted

from argparse import ArgumentParser
import math
import random

#some global things
scaling = 1.0
a = [0.0, 0.0, 0.0]
b = [0.0, 0.0, 0.0]
c = [0.0, 0.0, 0.0]
comment = ""
coord_type = ""

atom_types = []
atom_numbers = []

frame_atom_x = []
frame_atom_y = []
frame_atom_z = []

subst_atom_type = ""
subst_atom_number = 0
subst_atom_new_type = ""
subst_atom_delta = 0.0
subst_atom_radius = 1000.0

# command line option parser
parser = ArgumentParser()
parser.add_argument("-a", "--atom", dest="atom", help="Which atom to be substituted")
parser.add_argument("-e", "--element", dest="element", help="Name of the new atom that will be added", required = True)
parser.add_argument("-n", "--number", dest="number", help="How many atoms to be substituted", required=True)
parser.add_argument("-r", "--radius", dest="radius", help="Radius where the next atom will be put")
parser.add_argument("-d", "--displacement", dest="delta", help="Maximum change in new coordinates")

args = parser.parse_args()

#READ POSCAR FILE
input = open('POSCAR', 'r')

comment = input.readline()

scaling = float(input.readline().split()[0])

line = input.readline().split()
a[0] = float(line[0])
a[1] = float(line[1])
a[2] = float(line[2])

line = input.readline().split()
b[0] = float(line[0])
b[1] = float(line[1])
b[2] = float(line[2])

line = input.readline().split()
c[0] = float(line[0])
c[1] = float(line[1])
c[2] = float(line[2])

line = input.readline().split()

for atom in line:
	atom_types.append(atom)

line = input.readline().split()

for i in range(0,len(atom_types)):
	atom_numbers.append(int(line[i]))

coord_type = input.readline().split()[0]

scale_x = [scaling*a[0], scaling*a[1], scaling*a[2]]
scale_y = [scaling*b[0], scaling*b[1], scaling*b[2]]
scale_z = [scaling*c[0], scaling*c[1], scaling*c[2]]

for i in range(0,len(atom_types)):
	frame_atom_x.append([])
	frame_atom_y.append([])
	frame_atom_z.append([])

	for j in range(0,atom_numbers[i]):
		line = input.readline().split()
		ax = float(line[0])
		bx = float(line[1])
		cx = float(line[2])
		frame_atom_x[i].append( ax*scale_x[0] + bx*scale_y[0] + cx*scale_z[0])
		frame_atom_y[i].append( ax*scale_x[1] + bx*scale_y[1] + cx*scale_z[1])
		frame_atom_z[i].append( ax*scale_x[2] + bx*scale_y[2] + cx*scale_z[2])


#READ ARGUMENTS
subst_atom_number = int(args.number)
subst_atom_new_type = args.element

if(args.atom == None):
	subst_atom_type = atom_types[0]
else:
	subst_atom_type = args.atom

if(args.radius == None):
	subst_atom_radius = math.sqrt(math.pow(scale_x[0]+scale_y[0]+scale_z[0],2.0)+math.pow(scale_x[1]+scale_y[1]+scale_z[1],2.0)+math.pow(scale_x[2]+scale_y[2]+scale_z[2],2.0))*2.0
else:
	subst_atom_radius = float(args.radius)

if(args.delta == None):
	subst_atom_delta = 0.0
else:
	subst_atom_delta = float(args.delta)

input.close()

substitutable_atoms = 0
substitutable_index = 0

for i in range(0,len(atom_types)):
	if(atom_types[i] == subst_atom_type):
		substitutable_atoms = atom_numbers[i]
		substitutable_index = i
		break

if(subst_atom_number > substitutable_atoms):
	print "Number of substitute atoms larger than atoms available"

subst_atom_index = -1

for i in range(0,len(atom_types)):
	if(atom_types[i] == subst_atom_new_type):
		subst_atom_index = i
		break

if(subst_atom_index == -1):
	atom_types.append(subst_atom_new_type)
	atom_numbers.append(0)
	subst_atom_index = len(atom_numbers)-1
	frame_atom_x.append([])
	frame_atom_y.append([])
	frame_atom_z.append([])
	

center_x = 0.0
center_y = 0.0
center_z = 0.0

#here we select a random center around which new atoms will be introduced
i = random.randint(0, substitutable_atoms - 1)

center_x = frame_atom_x[substitutable_index][i]
center_y = frame_atom_y[substitutable_index][i]
center_z = frame_atom_z[substitutable_index][i]

tries = 10000

for i in range(0,subst_atom_number):
	for j in range(0, atom_numbers[substitutable_index]):
		l = random.randint(0, atom_numbers[substitutable_index] -1)

		#l = j
		distance = 0.0
		
		for o in range(-1,2):
			for p in range(-1,2):
				for r in range(-1,2):
					distance = math.pow(frame_atom_x[substitutable_index][l]-center_x+o*a[0]+p*b[0]+r*c[0],2)
					distance += math.pow(frame_atom_y[substitutable_index][l]-center_y+o*a[1]+p*b[1]+r*c[1],2)
					distance += math.pow(frame_atom_z[substitutable_index][l]-center_z+o*a[2]+p*b[2]+r*c[2],2)
					distance = math.sqrt(distance)

					if(distance <= subst_atom_radius):
						break

				if(distance <= subst_atom_radius):
					break

			if(distance <= subst_atom_radius):
				break

		if(distance <= subst_atom_radius):
			break

	frame_atom_x[substitutable_index][l] += random.uniform(-1.0,1.0)*subst_atom_delta
	frame_atom_y[substitutable_index][l] += random.uniform(-1.0,1.0)*subst_atom_delta
	frame_atom_z[substitutable_index][l] += random.uniform(-1.0,1.0)*subst_atom_delta

	frame_atom_x[subst_atom_index].append(frame_atom_x[substitutable_index][l])
	frame_atom_y[subst_atom_index].append(frame_atom_y[substitutable_index][l])
	frame_atom_z[subst_atom_index].append(frame_atom_z[substitutable_index][l])

	atom_numbers[subst_atom_index] += 1
	atom_numbers[substitutable_index] -= 1

	frame_atom_x[substitutable_index].pop(l)
	frame_atom_y[substitutable_index].pop(l)
	frame_atom_z[substitutable_index].pop(l)

	#if(j == atom_numbers[substitutable_index]):
	if(j == tries - 1):
		i = random.randint(0, atom_numbers[substitutable_index] - 1)

		center_x = frame_atom_x[substitutable_index][i]
		center_y = frame_atom_y[substitutable_index][i]
		center_z = frame_atom_z[substitutable_index][i]

#this is now the output
output = open('POSCAR_NEW', 'w')
output.write(comment)
output.write('   {0:10.6f}\n'.format(scaling))
output.write('   {0:10.6f}   {1:10.6f}   {2:10.6f}\n'.format(a[0], a[1], a[2]))
output.write('   {0:10.6f}   {1:10.6f}   {2:10.6f}\n'.format(b[0], b[1], b[2]))
output.write('   {0:10.6f}   {1:10.6f}   {2:10.6f}\n'.format(c[0], c[1], c[2]))

for i in atom_types:
	output.write('   {}'.format(i))
output.write('\n')

for i in atom_numbers:
	output.write('   {}'.format(i))
output.write('\n')

output.write('{}\n'.format(coord_type))

det = scale_x[0]*scale_y[1]*scale_z[2] + scale_x[1]*scale_y[2]*scale_z[0] + scale_x[2]*scale_y[0]*scale_z[1] - scale_x[2]*scale_y[1]*scale_z[0] - scale_x[1]*scale_y[0]*scale_z[2] - scale_x[0]*scale_y[2]*scale_z[1]

for i in range(0,len(atom_types)):
	for j in range(0, atom_numbers[i]):
		ax =  frame_atom_x[i][j]*(scale_y[1]*scale_z[2]-scale_y[2]*scale_z[1])
		ax += frame_atom_y[i][j]*(scale_y[2]*scale_z[0]-scale_y[0]*scale_z[2])
		ax += frame_atom_z[i][j]*(scale_y[0]*scale_z[1]-scale_y[1]*scale_z[0])
		ax =  ax/det
		
		bx =  frame_atom_x[i][j]*(scale_x[2]*scale_z[1]-scale_x[1]*scale_z[2])
		bx += frame_atom_y[i][j]*(scale_x[0]*scale_z[2]-scale_x[2]*scale_z[0])
		bx += frame_atom_z[i][j]*(scale_x[1]*scale_z[0]-scale_x[0]*scale_z[1])
		bx =  bx/det
		
		cx =  frame_atom_x[i][j]*(scale_x[1]*scale_y[2]-scale_x[2]*scale_y[1])
		cx += frame_atom_y[i][j]*(scale_x[2]*scale_y[0]-scale_x[0]*scale_y[2])
		cx += frame_atom_z[i][j]*(scale_x[0]*scale_y[1]-scale_x[1]*scale_y[0])
		cx =  cx/det
		
		output.write('{0:12.10f}  {1:12.10f}  {2:12.10f}\n'.format(ax, bx, cx))

output.close()

