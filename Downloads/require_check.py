#!/usr/bin/python

#    require_check.py
#    Copyright (C) 2014  LeukemiaNetworks

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License Version 2 
#    as published by the Free Software Foundation.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


# This program is intended to aid users in determining 
# which python packages need to be installed to run multidatanetanalysis.py

installed_list = []
need_list= []

# can it be installed
def mod (name):
	# if it can be imported, it is already installed, else output that it needs to be downloaded separately 
	try:
		__import__(name)
		installed_list.append(name)
	except ImportError:
		need_list.append(name)

# run mod method to determine what is installed
mod('os')
mod('shutil')
mod('timeit')
mod('gzip')
mod('zipfile')
mod('numpy')
mod('scipy')
mod('matplotlib')
mod('Tkinter')
mod('tkFileDialog')
mod('zipfile')

print "Already Installed\n-----------------"
for i in installed_list:
	print '\t',i
print"\nNeed to Install\n-----------------"
for j in need_list:
	print '\t',j
print "\n"
if len(need_list) ==0:
	print "No files need to be installed"
