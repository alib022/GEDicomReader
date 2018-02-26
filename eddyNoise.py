import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse
from clint.textui import colored
from vtk.util import numpy_support
from scipy.ndimage.filters import uniform_filter
from rolling_window import rolling_window


import matplotlib.pyplot as plt

def eddyCurrentCorrection(UOrg, VOrg, WOrg, randNoiseThreshold, eddyCurrentThreshold, eddyOrder):


    USTD = numpy.zeros(UOrg.shape)
    VSTD = numpy.zeros(VOrg.shape)
    WSTD = numpy.zeros(WOrg.shape)
    

    for tIter in range(UOrg.shape[3]):
        for kIter in range(UOrg.shape[2]):
#            URoll = pandas.Series(UOrg[:,:,kIter,tIter])
#            VRoll = pandas.Series(VOrg[:,:,kIter,tIter])
#            WRoll = pandas.Series(WOrg[:,:,kIter,tIter])
            
#            USTD[:,:,kIter,tIter]= URoll.rolling(3, center=True, win_type='gaussian', min_periods=(3//2)).std().values
#            VSTD[:,:,kIter,tIter]= VRoll.rolling(3, center=True, win_type='gaussian', min_periods=(3//2)).std().values
#            WSTD[:,:,kIter,tIter]= WRoll.rolling(3, center=True, win_type='gaussian', min_periods=(3//2)).std().values

            USTD[:,:,kIter,tIter] = numpy.std(rolling_window(UOrg[:,:,kIter,tIter], 3), -1)
            VSTD[:,:,kIter,tIter] = numpy.std(rolling_window(VOrg[:,:,kIter,tIter], 3), -1)
            WSTD[:,:,kIter,tIter] = numpy.std(rolling_window(WOrg[:,:,kIter,tIter], 3), -1)
    

    #USTD = numpy.std(rolling_window(UOrg, 3), -1)
    #VSTD = numpy.std(rolling_window(VOrg, 3), -1)
    #WSTD = numpy.std(rolling_window(WOrg, 3), -1)
    
    print("Ustd Max: ")
    print(USTD.max())
    print("Vstd Max: ")
    print(VSTD.max())
    print("Wstd Max: ")
    print(WSTD.max())

    print("Ustd Min: ")
    print(USTD.min())
    print("Vstd Min: ")
    print(VSTD.min())
    print("Wstd Min: ")
    print(WSTD.min())

    print("USTD size: ")
    print(USTD.shape)

    print(USTD[:,:,1,1])

    
    

    eddyMaskIndicesU = numpy.where(USTD < (eddyCurrentThreshold*USTD.max()/100))
    eddyMaskIndicesV = numpy.where(VSTD < (eddyCurrentThreshold*VSTD.max()/100))
    eddyMaskIndicesW = numpy.where(WSTD < (eddyCurrentThreshold*WSTD.max()/100))

    print(eddyMaskIndicesU)
 
    

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.imshow(UOrg[:,:,20,1])
    ax2.imshow(eddyMaskIndicesU[:,:,20,1])
    plt.show()

    
    sys.exit()
   # print(eddyMaskIndicesU.shape)

    del Ustd, Vstd, Wstd  
    
    flowCorrected = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2], 3, UOrg.shape[3]])
   
       
    xInit = numpy.linspace(0, UOrg.shape[0], UOrg.shape[0])
    yInit = numpy.linspace(0, UOrg.shape[1], UOrg.shape[1])
    
    X, Y = numpy.meshgrid(xInit, yInit, sparse=False, indexing='ij')
    
    X2=X*X
    Y2=Y*Y
    XY=X*Y
    
    
    plainU = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2], UOrg.shape[3]])
    plainV = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2], UOrg.shape[3]])
    plainW = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2], UOrg.shape[3]])
    
    if eddyOrder == 1:
        # best-fit linear plane
        for kIter in range(UOrg.shape[3]):
            for iIter in range(UOrg.shape[2]):
            
                BU = UOrg[:, :, iIter, kIter].ravel()[eddyMaskIndicesU[:,:,iIter,kIter]]
                BV = VOrg[:, :, iIter, kIter].ravel()[eddyMaskIndicesV[:,:,iIter,kIter]]
                BW = WOrg[:, :, iIter, kIter].ravel()[eddyMaskIndicesW[:,:,iIter,kIter]]
                    
                D = numpy.ones((len(eddyMaskIndicesU[:,:,iIter,kIter]), 3))
                D[:, 0] = X.ravel()[eddyMaskIndicesU[:,:,iIter,kIter]]
                D[:, 1] = Y.ravel()[eddyMaskIndicesU[:,:,iIter,kIter]]
            
       
                CU,_,_,_ = scipy.linalg.lstsq(D, BU)    # coefficients
                CV,_,_,_ = scipy.linalg.lstsq(D, BV)    # coefficients
                CW,_,_,_ = scipy.linalg.lstsq(D, BW)
        
                # evaluate it on grid
                plainU[:,:,iIter, kIter] = CU[0]*X + CU[1]*Y + CU[2]
                plainV[:,:,iIter, kIter] = CV[0]*X + CV[1]*Y + CV[2]
                plainW[:,:,iIter, kIter] = CW[0]*X + CW[1]*Y + CW[2]
            
            
        
    
    elif eddyOrder == 2:
        # best-fit quadratic curve
        
        for iIter in range(plainU.shape[2]):
            
            eddyMaskTemp = eddyMask[:,:,iIter]
            indEddyMask = numpy.argwhere(eddyMaskTemp.ravel()).ravel()
            BU = UOrgMean[:, :, iIter].ravel()[indEddyMask]
            BV = VOrgMean[:, :, iIter].ravel()[indEddyMask]
            BW = WOrgMean[:, :, iIter].ravel()[indEddyMask]
                    
            D = numpy.ones((len(indEddyMask), 6))
            D[:, 0] = X.ravel()[indEddyMask]
            D[:, 1] = Y.ravel()[indEddyMask]
            D[:, 2] = XY.ravel()[indEddyMask]
            D[:, 3] = X2.ravel()[indEddyMask]
            D[:, 4] = Y2.ravel()[indEddyMask]
            
       
            CU,_,_,_ = scipy.linalg.lstsq(D, BU)    # coefficients
            CV,_,_,_ = scipy.linalg.lstsq(D, BV)    # coefficients
            CW,_,_,_ = scipy.linalg.lstsq(D, BW)
        
            # evaluate it on grid
            plainU[:,:,iIter] = CU[0]*X + CU[1]*Y + CU[2]*XY + CU[3]*X2 + CU[4]*Y2 + CU[5]
            plainV[:,:,iIter] = CV[0]*X + CV[1]*Y + CV[2]*XY + CV[3]*X2 + CV[4]*Y2 + CV[5]
            plainW[:,:,iIter] = CW[0]*X + CW[1]*Y + CW[2]*XY + CW[3]*X2 + CW[4]*Y2 + CW[5]
            
    
    #for k in range(UOrg.shape[3]):
        
    flowCorrected[:, :, :,0] = UOrg - plainU
    flowCorrected[:, :, :,1] = VOrg - plainV
    flowCorrected[:, :, :,2] = WOrg - plainW
    
    
    return flowCorrected

def randNoise(args, UOrg, VOrg, WOrg, randThre, load=1, save=0):

    SDU = stdWindow(UOrg, 5)
    SDV = stdWindow(VOrg, 5)
    SDW = stdWindow(WOrg, 5)
    if save:
       numpy.save(args.output + '/SDU', SDU)
       numpy.save(args.output + '/SDV', SDV)        
       numpy.save(args.output + '/SDW', SDW)
            

    
    flowCorrected = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2],3,UOrg.shape[3]])


    UOrg[numpy.where(SDU < randThre*numpy.max(SDU))] = 0
    VOrg[numpy.where(SDV < randThre*numpy.max(SDV))] = 0
    WOrg[numpy.where(SDW < randThre*numpy.max(SDW))] = 0

    flowCorrected[:,:,:,0] = UOrg
    flowCorrected[:,:,:,1] = VOrg
    flowCorrected[:,:,:,2] = WOrg

    return flowCorrected





