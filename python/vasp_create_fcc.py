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

# first argument is lattice length
# second specifies name of the system
# third argument specifies element type (optional)
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-a", "--lattice", dest="lattice", help="Lattice constant")
parser.add_argument("-e", "--element", dest="element", help="Atom type")
parser.add_argument("-s", "--system", dest="system", help="System name")

args = parser.parse_args()

if args.lattice == None:
	lattice = "1.000"
else:
	lattice=args.lattice

if args.element == None:
	element = "X"
else:
	element=args.element

if args.system == None:
	system = "System"
else:
	system=args.system

print system
print "   " + lattice
#bcc lattice vectors
print "%.10f   %.10f   %.10f"%(1.0,0.0,0.0)
print "%.10f   %.10f   %.10f"%(0.0,1.0,0.0)
print "%.10f   %.10f   %.10f"%(0.0,0.0,1.0)
print "   " + element
print "   " + str(4)
print "Direct"
print "%.10f   %.10f   %.10f"%(0.0,0.0,0.0)
print "%.10f   %.10f   %.10f"%(0.5,0.5,0.0)
print "%.10f   %.10f   %.10f"%(0.5,0.0,0.5)
print "%.10f   %.10f   %.10f"%(0.0,0.5,0.5)
print

