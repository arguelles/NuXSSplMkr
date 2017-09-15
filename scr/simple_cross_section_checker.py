import matplotlib.pyplot as plt
import sys,os
import numpy as np
import operator

from photospline import spglam as glam
from photospline.glam.glam import grideval
from photospline.glam.bspline import bspline
from photospline import splinefitstable

if len(sys.argv) < 4:
    print sys.argv
    print "Usage: cross_section_checker.py PDFNAME"

pdfname = sys.argv[1]
variation = sys.argv[2]
neutype = sys.argv[3]

# reading spline and table
print "reading spline"
result = splinefitstable.read('../fits/dsdxdy-'+neutype+'-N-cc-' + pdfname + '_' + variation + '.fits')
print "reading data"
datas = np.genfromtxt("../data/dsdxdy-"+neutype+"-N-cc-"+pdfname+"_"+variation+".dat")

shape = tuple(np.unique(datas[:,i]).size for i in range(3))
datas = datas.reshape(shape + (datas.size/reduce(operator.mul,shape),))

print "looping around"
q2_array = datas[:,0,0,0]
x_array = datas[0,:,0,1]
y_array = datas[0,0,:,2]
dsigdxdy_array = datas[:,:,:,3]

q2 = (1.0e3)**2
x = 1.0e-3
y = 1.0e-3

for i,q2 in enumerate(q2_array):
    for j,x in enumerate(x_array):
        for k,y in enumerate(y_array):
            #print i,j,k
            dsdxdy = datas[i,j,k,3]
            spl = np.power(10.,glam.grideval(result,[np.array([np.log10(q2)]),np.array([np.log10(x)]),np.array([np.log10(y)])])[0])[0][0]
            if np.abs(spl-dsdxdy)/dsdxdy > 1.0e-2:
                print dsdxdy,spl,spl/dsdxdy

