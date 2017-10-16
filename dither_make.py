import numpy as np
import math as mt
import matplotlib.pyplot as plt

def dither_make(xsize, ysize, npts, mindist, ntol, outfile):

    '''Written by Gregory Rudnick 16 October 2017

    PURPOSE:

    Make a set of x and y coordinates that are random but have a
    minimum distance to a set of preceeding coordinates.

    This isn't being done in an especially clever way.  Just start at
    the first point and keep generating random following points, only
    accepting them if they satisfy the minimum distance requirement

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
        print(itest,i-itest,distcheck)
        while ((itest >= 0) and (i-itest <= ntol) and distcheck):
            dtest = mt.sqrt( (xtest - x[itest])**2 + (ytest - y[itest])**2)
            
            #if it satisfies the condition then go one further back.
            #If it doesn't then mark as a failure`
            if dtest>=mindist:
                print("true",i,itest,dtest)
                itest = itest -1
            else:
                distcheck = False
                print("false",i,itest,dtest)
                
        #if this point satisfied the proximity condition then keep
        #this point.  Otherwise repeat again
        if distcheck:
            x = np.append(x,xtest)
            y = np.append(y,ytest)
            i = i + 1

    #print(x,y)
    print("hello")
    plt.clf()
    fig, dithplot = plt.subplots()
    dithplot.plot(x,y)
    dithplot.plot(x,y,'gs')
    dithplot.axis([-50., 50., -50., 50.])
    dithplot.set_xlabel('x')
    dithplot.set_ylabel('y')
    plt.savefig('test.pdf')
    #plt.show()
    

