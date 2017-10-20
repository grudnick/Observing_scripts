import numpy as np
import math as mt
import matplotlib.pyplot as plt

def dither_batch():
    xsize = 60.0
    ysize = 60.0
    npts = 117
    mindist = 25.0
    ntol = 3
    outfile = "dither_Jband.txt"
    dither_make(xsize, ysize, npts, mindist, ntol, outfile)
   

def dither_make(xsize, ysize, npts, mindist, ntol, outfile):

    '''Written by Gregory Rudnick 16 October 2017

    PURPOSE:

    Make a set of x and y coordinates that are random but have a
    minimum distance to a set of preceeding coordinates.

    This isn't being done in an especially clever way.  Just start at
    the first point and keep generating random following points, only
    accepting them if they satisfy the minimum distance requirement

    INPUT PARAMETERS:

    xsize,ysize: size of dither box in arcsec

    npts: number of dither positions

    mindist: The minimum distance between successive dithers in arcsec

    ntol: the number of dithers in the past for which the minimum
    distance is enforced.  Should roughly be the half of the number of
    dithers you will use to estimate the sky.

    outfile: the file with the output positions

    '''

    #initialize arrays with initial point
    x = np.array([np.random.rand(1) * xsize - xsize/2.],dtype=float)   
    y = np.array([np.random.rand(1) * ysize - ysize/2.],dtype=float)   

    i = 1
    while i<npts:
        #make random position within the box
        xtest = np.random.rand(1) * xsize - xsize/2.
        ytest = np.random.rand(1) * ysize - ysize/2.
        itest = i - 1
        
        #this will stay true as long as the preceeding points aren't
        #closer than mindist
        distcheck = True          
        while ((itest >= 0) and (i-itest <= ntol) and distcheck):
            dtest = mt.sqrt( (xtest - x[itest])**2 + (ytest - y[itest])**2)
            
            #if it satisfies the condition then go one further back.
            #If it doesn't then mark as a failure`
            if dtest>=mindist:
                #print("true",i,itest,dtest)
                itest = itest -1
            else:
                distcheck = False
                #print("false",i,itest,dtest)
                
        #if this point satisfied the proximity condition then keep
        #this point.  Otherwise repeat again
        if distcheck:
            x = np.append(x,xtest)
            y = np.append(y,ytest)
            i += 1

    #now go through and calculate relative shifts.  
    xrel = np.array(x[0])
    yrel = np.array(y[0])
    i = 1
    while i < npts:
        xrel = np.append(xrel, x[i] - x[i - 1])
        yrel = np.append(yrel, y[i] - y[i - 1])
        i += 1

    #the last cycle should put the telescope back at the original pointing
    xrel = np.append(xrel, -x[npts - 1])
    yrel = np.append(yrel, -y[npts - 1])

    print(x,y)
    print("")
    print(xrel,yrel)
    sumx = np.sum(xrel)
    sumy = np.sum(yrel)
    print("")
    print("these should be zero if the shifts move you back to the origin at the end")
    print(sumx,sumy)
    
            
    plt.clf()
    fig, dithplot = plt.subplots()
    dithplot.plot(x,y)
    dithplot.plot(x,y,'gs')
    #dithplot.plot(xrel,yrel,'ro',markersize=5)
    dithplot.axis([-50., 50., -50., 50.])
    dithplot.set_xlabel('x')
    dithplot.set_ylabel('y')
    plotfile = outfile + '.pdf'
    plt.savefig(plotfile)
    #plt.show()
    
    #need to write to output file in x,y format

