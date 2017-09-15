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

if len(sys.argv) < 4:
    print sys.argv
    print "Usage: cross_section_checker.py PDFNAME variation neutype"

pdfname = sys.argv[1]
variation = sys.argv[2]
neutype = sys.argv[3]

# reading spline and table
result = splinefitstable.read('../fits/dsdxdy-'+neutype+'-N-cc-' + pdfname + '_' + variation + '.fits')
datas = np.genfromtxt("../data/dsdxdy-"+neutype+"-N-cc-"+pdfname+"_"+variation+".dat")

shape = tuple(np.unique(datas[:,i]).size for i in range(3))
datas = datas.reshape(shape + (datas.size/reduce(operator.mul,shape),))

q2_array = datas[:,0,0,0]
x_array = datas[0,:,0,1]
y_array = datas[0,0,:,2]
dsigdxdy_array = datas[:,:,:,3]

for i,q2 in enumerate(q2_array):
    if not (i%20 == 0):
        continue
    for j,x in enumerate(x_array):
        if not (j%20 == 0):
            continue
        fig = plt.figure(figsize = (8,6))
        dsdxdy = datas[i,j,:,3]

        plt.plot(y_array,dsdxdy, label = "CalculationValues", color = "blue", lw = 2)
        spllist = np.power(10.,glam.grideval(result,[np.array([np.log10(q2)]),np.array([np.log10(x)]),np.log10(y_array)])[0])
        #print spllist[0]
        plt.plot(y_array,spllist[0], label = "InterValues", color = "red", lw =2, ls = "dashed")
        plt.semilogx()
        plt.title(r"$Q2 = "+ str(q2) + "," + " x = " + str(x) + "$")
        plt.legend(loc = "lower left")
        plt.xlabel("y")
        plt.ylabel("dsigma/dxdy")

        plt.savefig("./plots/dsdxdy_"+pdfname+"_"+neutype+"_"+str(q2)+"_"+str(x)+"_"+variation+".png")
        print "Saved: ./plots/dsdxdy_"+pdfname+"_"+neutype+"_"+str(q2)+"_"+str(x)+"_"+variation+".png"
        plt.close(fig)

