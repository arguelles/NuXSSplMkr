import matplotlib.pyplot as plt
import sys,os
import numpy as np
import operator
try:
	# Fast sparse-matrix implementation
	from photospline import spglam as glam
except ImportError:
	# Slow (yet readable) pure-python implementation
	from photospline.glam import glam
from photospline.glam.glam import grideval
from photospline.glam.bspline import bspline
from photospline import splinefitstable

if len(sys.argv) < 3:
    print sys.argv
    print "Usage: cross_section_checker.py PDFNAME neutype"

pdfname = sys.argv[1]
neutype = sys.argv[2]

# reading spline and table
variation = "central"
data_central = np.genfromtxt("../data/dsdxdy-"+neutype+"-N-cc-"+pdfname+"_"+variation+".dat")
variation = "minus"
data_minus = np.genfromtxt("../data/dsdxdy-"+neutype+"-N-cc-"+pdfname+"_"+variation+".dat")

shape = tuple(np.unique(data_central[:,i]).size for i in range(3))
data_central = data_central.reshape(shape + (data_central.size/reduce(operator.mul,shape),))

shape = tuple(np.unique(data_minus[:,i]).size for i in range(3))
data_minus = data_minus.reshape(shape + (data_minus.size/reduce(operator.mul,shape),))

q2_array = data_central[:,0,0,0]
x_array = data_central[0,:,0,1]
y_array = data_central[0,0,:,2]
dsigdxdy_array = data_central[:,:,:,3]

data_new = np.empty([len(q2_array)*len(x_array)*len(y_array),4])
#data_new = np.empty([[len(q2_array)*len(x_array)*len(y_array),4], dtype = float)

l = 0
for i,q2 in enumerate(q2_array):
    for j,x in enumerate(x_array):
        for k,y in enumerate(y_array):
            dsdxdy = data_minus[i,j,k][-1]
            data_new[l] = [q2,x,y,dsdxdy]
            if dsdxdy < 0:
                data_new[l][-1] = 1.0e-3*data_central[i,j,k][-1]
            l = l+1

    print i

np.savetxt("../data/dsdxdy-"+neutype+"-N-cc-"+pdfname+"_minus_fix.dat",data_new,delimiter = ' ')

