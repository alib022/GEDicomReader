import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse, math
from clint.textui import colored
from vtk.util import numpy_support

def eddyCurrentCorrection(UOrg, VOrg, WOrg, randNoiseThreshold, eddyCurrentThreshold, eddyOrder):


    Ustd = numpy.std(UOrg, axis=3)
    Vstd = numpy.std(VOrg, axis=3)
    Wstd = numpy.std(WOrg, axis=3)



    #### Random noise correction

    #randNoiseThreshold = 60.0
    noiseMask = numpy.logical_and(numpy.logical_and(Ustd > randNoiseThreshold, Vstd > randNoiseThreshold), Wstd > randNoiseThreshold)
    #noiseMask[:, :cut.start] = False
    #noiseMask[:, cut.stop:] = False

    #### Eddy currents correction

    #eddyMask = Ustd < 20  
    eddyMask = Ustd <   eddyCurrentThreshold

    
    del Ustd, Vstd, Wstd  
    
    flowCorrected = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2],3,UOrg.shape[3]])
    shape = (UOrg.shape[0], UOrg.shape[1])
       
    xInit = numpy.linspace(0, shape[0], shape[0])
    yInit = numpy.linspace(0, shape[1], shape[1])
    
    X, Y = numpy.meshgrid(xInit, yInit, sparse=False, indexing='ij')
    
    X2=X*X
    Y2=Y*Y
    XY=X*Y
    
    
    plainU = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2]])
    plainV = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2]])
    plainW = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2]])
    
    #flowMean = numpy.mean(flowData, axis=4)
    UOrgMean = numpy.mean(UOrg, axis=3)
    VOrgMean = numpy.mean(VOrg, axis=3)
    WOrgMean = numpy.mean(WOrg, axis=3)
    
    if eddyOrder == 1:
        # best-fit linear plane
        for iIter in range(plainU.shape[2]):
            
            #maskDataTemp = maskData[:,:,iIter ]
            eddyMaskTemp = eddyMask[:,:,iIter]
            indEddyMask = numpy.argwhere(eddyMaskTemp.ravel()).ravel()
            # noiseMaskInd = np.argwhere(noiseMask.ravel()).ravel()
            
            
            BU = UOrgMean[:, :, iIter].ravel()[indEddyMask]
            BV = VOrgMean[:, :, iIter].ravel()[indEddyMask]
            BW = WOrgMean[:, :, iIter].ravel()[indEddyMask]
                    
            D = numpy.ones((len(indEddyMask), 3))
            D[:, 0] = X.ravel()[indEddyMask]
            D[:, 1] = Y.ravel()[indEddyMask]
            
       
            CU,_,_,_ = scipy.linalg.lstsq(D, BU)    # coefficients
            CV,_,_,_ = scipy.linalg.lstsq(D, BV)    # coefficients
            CW,_,_,_ = scipy.linalg.lstsq(D, BW)
        
            # evaluate it on grid
            plainU[:,:,iIter] = CU[0]*X + CU[1]*Y + CU[2]
            plainV[:,:,iIter] = CV[0]*X + CV[1]*Y + CV[2]
            plainW[:,:,iIter] = CW[0]*X + CW[1]*Y + CW[2]
            
            
        
    
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
            
    
    for k in range(UOrg.shape[3]):
        
        flowCorrected[:, :, :,0, k] = UOrg[:, :, :, k] - plainU
        flowCorrected[:, :, :,1, k] = VOrg[:, :, :, k] - plainV
        flowCorrected[:, :, :,2, k] = WOrg[:, :, :, k] - plainW
    
    
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


    UOrg[numpy.where(SDU > randThre*numpy.max(SDU))] = 0
    VOrg[numpy.where(SDV > randThre*numpy.max(SDV))] = 0
    WOrg[numpy.where(SDW > randThre*numpy.max(SDW))] = 0

    flowCorrected[:,:,:,0] = UOrg
    flowCorrected[:,:,:,1] = VOrg
    flowCorrected[:,:,:,2] = WOrg

    return flowCorrected

def stdWindow(UInput, WW=5):
    UShape = UInput.shape
    print(UShape)
    maskNoise = numpy.zeros([UShape[0], UShape[1], UShape[2]], dtype=bool)
    SD = numpy.zeros((UShape[0]-WW+1, UShape[1]-WW+1,UShape[2], UShape[3]))
    for kIter in range(UShape[2]):
        print(kIter)
        Mean = numpy.zeros((UShape[0]-WW+1, UShape[1]-WW+1, UShape[3]))
        
        for iIter in range(WW):
            for jIter in range(WW):
                Mean += UInput[iIter:UShape[0]-WW+1+iIter, jIter:UShape[1]-WW+1+jIter, kIter]
                SD[:,:,kIter] += UInput[iIter:UShape[0]-WW+1+iIter, jIter:UShape[1]-WW+1+jIter, kIter]**2 
        
        Mean /= WW**2
        SD[:,:,kIter] = numpy.sqrt(SD[:,:,kIter]/WW**2 - Mean**2)
    
    
    return SD


 
