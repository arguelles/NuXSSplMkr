import numpy
import operator
import re
import warnings

# Fast sparse-matrix implementation
#from icecube.photospline import spglam as glam 
#from icecube.photospline.glam.glam import grideval
#from icecube.photospline.glam.bspline import bspline
#from icecube.photospline import splinefitstable

# Fast sparse-matrix implementation
from photospline import spglam as glam
from photospline.glam.glam import grideval
from photospline.glam.bspline import bspline
from photospline import splinefitstable
import os

"""
Script to make structure functions spline tables.
C.A. Arguelles Delgado - aug.03.14
"""

#warnings.filterwarnings('error')

def ModLog10(x):
    if x <= 0. :
        print x
        return -50
    else:
        return numpy.log10(x)

ModLog10 = numpy.vectorize(ModLog10)

def SplineFitMaker1D(filename, scale = 'lin', prefix = '', skip_header = 0, column = 1, N = 50, outname = "", oscale = 'lin'):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x), asummes x to be the first column.
    """
    if(column < 1):
        print "Error: column < 1."
        exit()

    datas = numpy.loadtxt(filename, skiprows = skip_header)

    f = lambda x : x;
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print "Error: unknown scale."
        exit()

    of = lambda x : x;
    if oscale == "log":
        of = lambda x : ModLog10(x)
    elif oscale == "lin":
        pass
    else:
        print "Error: unknown scale."
        exit()

    #shape = tuple(numpy.unique(datas[:,i]).size for i in range(2))
    #datas = datas.reshape(shape + (datas.size/reduce(operator.mul,shape),))

    x = f(datas[:,0])
    z = of(datas[:,column])

    knots = [numpy.linspace(x.min()-1,x.max()+1,N,endpoint = True)]
    order = 2
    smooth = 1.0e-15

    weight = numpy.ones(z.shape)
    result = glam.fit(z,weight,[x],knots,order,smooth)

    # creatiing new filename and saving
    if outname == "":
        nfilename = filename.split('.').pop() + ".fits"
    else :
        nfilename = outname
    # save the result to a FITS table
    if os.path.exists(prefix+nfilename):
        os.unlink(prefix+nfilename)
    # this file can later be loaded and evaluated with the photospline C library
    splinefitstable.write(result, prefix+nfilename)
    print "Done. Generated :"  + prefix+nfilename

def SplineFitMaker2D(filename, scale = 'lin', prefix = '', skip_header = 0, column = 2, N = 50, outname = ""):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y), asummes x y to be the first two columns.
    """
    if(column < 2):
        print "Error: column < 2."
        exit()

    datas = numpy.loadtxt(filename, skiprows = skip_header)

    f = lambda x : x;
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print "Error: unknown scale."
        exit()

    shape = tuple(numpy.unique(datas[:,i]).size for i in range(2))
    datas = datas.reshape(shape + (datas.size/reduce(operator.mul,shape),))

    x = f(datas[:,0,0])
    y = f(datas[0,:,1])
    z = datas[:,:,column]

    knots = [numpy.linspace(x.min()-1,x.max()+1,N,endpoint = True),numpy.linspace(y.min()-1,y.max()+1,N,endpoint = True)]
    #knots = [x,y]
    order = 2
    #smooth = 1.0e-5
    smooth = 1.0e-15

    weight = numpy.ones(z.shape)
    #weight = 1+zz
    result = glam.fit(z,weight,[x,y],knots,order,smooth)

    # creatiing new filename and saving
    if outname == "":
        nfilename = filename.split('.').pop() + ".fits"
    else :
        nfilename = outname
    # save the result to a FITS table
    if os.path.exists(prefix+nfilename):
        os.unlink(prefix+nfilename)
    # this file can later be loaded and evaluated with the photospline C library
    splinefitstable.write(result, prefix+nfilename)
    print "Done. Generated :"  + prefix+nfilename

def SplineFitMaker3D(filename, scale = 'lin', prefix = '', skip_header = 0, column = 2, N = 50, outname = "", oscale = 'lin'):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y,w), asummes x/y/w to be the first/second/third column.
    """
    if(column < 3):
        print "Error: column < 3."
        exit()

    datas = numpy.loadtxt(filename, skiprows = skip_header)

    f = lambda x : x;
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print "Error: unknown scale."
        exit()

    of = lambda x : x;
    if oscale == "log":
        of = lambda x : ModLog10(x)
    elif oscale == "lin":
        pass
    else:
        print "Error: unknown scale."
        exit()

    shape = tuple(numpy.unique(datas[:,i]).size for i in range(3))
    datas = datas.reshape(shape + (datas.size/reduce(operator.mul,shape),))

    x = f(datas[:,0,0,0])
    y = f(datas[0,:,0,1])
    w = f(datas[0,0,:,2])
    z = of(datas[:,:,:,column])

    knots = [numpy.linspace(x.min()-1,x.max()+1,N,endpoint = True),
             numpy.linspace(y.min()-1,y.max()+1,N,endpoint = True),
             numpy.linspace(w.min()-1,w.max()+1,N,endpoint = True)]
    #knots = [x,y]
    order = 2
    #smooth = 1.0e-5
    smooth = 1.0e-15

    weight = numpy.ones(z.shape)
    #weight = 1+zz
    result = glam.fit(z,weight,[x,y,w],knots,order,smooth)

    # creatiing new filename and saving
    if outname == "":
        nfilename = filename.split('.').pop() + ".fits"
    else :
        nfilename = outname
    # save the result to a FITS table
    if os.path.exists(prefix+nfilename):
        os.unlink(prefix+nfilename)
    # this file can later be loaded and evaluated with the photospline C library
    splinefitstable.write(result, prefix+nfilename)
    print "Done. Generated :"  + prefix+nfilename

if __name__ == "__main__":
    inpath = "/home/carguelles/NuXSSplMkr/data/newxs_2017/"
    outpath = "/home/carguelles/NuXSSplMkr/fits/newxs_2017/"

    neutrino_type = ['numu','numubar']

    pdf_list = ['CT10nlo_central','CT10nlo_minus','CT10nlo_plus',
                'HERAPDF15NLO_EIG_central','HERAPDF15NLO_EIG_minus','HERAPDF15NLO_EIG_plus',
                'NNPDF23_nlo_as_0118_central','NNPDF23_nlo_as_0118_minus','NNPDF23_nlo_as_0118_plus']

    #pdf_list = ['NNPDF23_nlo_as_0118_central','NNPDF23_nlo_as_0118_minus','NNPDF23_nlo_as_0118_plus']
    pdf_list = ['HERAPDF15NLO_EIG_central']

    #pdf = "HERAPDF15NLO_EIG_central"
    #filename = "dsdxdy-numu-N-cc-"+pdf
    #SplineFitMaker3D(inpath + filename + ".dat", outname = filename + ".fits",
    #        scale = 'log',prefix = outpath, N = 50, column = 3, oscale = 'log')
#
#    filename = "dsdxdy-numubar-N-cc-"+pdf
#    SplineFitMaker3D(inpath + filename + ".dat", outname = filename + ".fits",
#            scale = 'log',prefix = outpath, N = 50, column = 3, oscale = 'log' )

#    quit()
    for int_type in ["cc","nc"]:
        for pdf in pdf_list:
            for neutype in neutrino_type:
                filename = "sigma-"+neutype+"-N-"+int_type+"-"+pdf
                print "processing: "+filename
                SplineFitMaker1D(inpath + filename + ".dat", outname = filename + ".fits",
                        scale = 'log',prefix = outpath, N = 65, column = 1, oscale = 'log')

    exit()

    for int_type in ["cc","nc"]:
        for pdf in pdf_list:
            for neutype in neutrino_type:
                filename = "dsdxdy-"+neutype+"-N-"+int_type+"-"+pdf
                print "processing: "+filename
                SplineFitMaker3D(inpath + filename + ".dat", outname = filename + ".fits",
                        scale = 'log',prefix = outpath, N = 65, column = 3, oscale = 'log')


