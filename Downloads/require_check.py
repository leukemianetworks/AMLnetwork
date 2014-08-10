#!/usr/bin/python

installed_list = []
need_list= []
def mod (name):
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
