import scipy.io, numpy, saveVTK, warnings
#from clint.textui import colored
#from vtk.util import numpy_support
#from scipy.ndimage.filters import uniform_filter
from rolling_window import rolling_window



import matplotlib.pyplot as plt
#from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


warnings.simplefilter(action='ignore', category=FutureWarning)

def eddyCurrentCorrection(UOrg, VOrg, WOrg, magData, eddyCurrentThreshold=15, eddyOrder=2, STDPower=2,   plotEddyPlane=1, plotPlain=20):


    USTD = numpy.zeros((UOrg.shape[0],UOrg.shape[1]))
    VSTD = numpy.zeros(USTD.shape)
    WSTD = numpy.zeros(USTD.shape)
    
    flowCorrected = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2], 3, UOrg.shape[3]])
          
    xInit = numpy.linspace(0, UOrg.shape[0], UOrg.shape[0])
    yInit = numpy.linspace(0, UOrg.shape[1], UOrg.shape[1])
    
    X, Y = numpy.meshgrid(xInit, yInit, sparse=False, indexing='ij')
    
    
    X2=X*X
    Y2=Y*Y
    XY=X*Y
    
    
    plainU = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2]])
    plainV = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2]])
    plainW = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2]])
    
    if eddyOrder == 1:
        
        for iIter in range(UOrg.shape[2]):
            
                # best-fit linear plane
                UFit = UOrg[:, :, iIter, -1].copy()
                VFit = VOrg[:, :, iIter, -1].copy()
                WFit = WOrg[:, :, iIter, -1].copy()
                
                
                magDataSelected = magData[:,:,iIter].copy()
                
                USTD[1:(UOrg.shape[0]-1),1:(UOrg.shape[1]-1)] = numpy.std(rolling_window(UOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
                VSTD[1:(VOrg.shape[0]-1),1:(VOrg.shape[1]-1)] = numpy.std(rolling_window(VOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
                WSTD[1:(WOrg.shape[0]-1),1:(WOrg.shape[1]-1)] = numpy.std(rolling_window(WOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))

                USTDSelectInd = numpy.where(USTD > (eddyCurrentThreshold/100) * USTD.max(), USTD, 0)
                VSTDSelectInd = numpy.where(VSTD > (eddyCurrentThreshold/100) * VSTD.max(), VSTD, 0)
                WSTDSelectInd = numpy.where(WSTD > (eddyCurrentThreshold/100) * WSTD.max(), WSTD, 0)
    
                USTDSelected = USTD[USTDSelectInd].copy()
                VSTDSelected = VSTD[VSTDSelectInd].copy()
                WSTDSelected = WSTD[WSTDSelectInd].copy()
            
                with numpy.errstate(invalid='ignore'):
                    weightU = magDataSelected / (USTD) ** STDPower
                    weightV = magDataSelected / (VSTD) ** STDPower
                    weightW = magDataSelected / (WSTD) ** STDPower
           
                weightU[numpy.isnan(weightU)] = 0
                weightV[numpy.isnan(weightV)] = 0
                weightW[numpy.isnan(weightW)] = 0
                weightU[numpy.isinf(weightU)] = 0
                weightV[numpy.isinf(weightV)] = 0
                weightW[numpy.isinf(weightW)] = 0
    
                     
                notZeroIndU = numpy.where(numpy.logical_and(weightU > 0, USTDSelectInd))
                notZeroIndV = numpy.where(numpy.logical_and(weightV > 0, VSTDSelectInd))
                notZeroIndW = numpy.where(numpy.logical_and(weightW > 0, WSTDSelectInd))
                
                BU = UFit[notZeroIndU]
                BV = VFit[notZeroIndV]
                BW = WFit[notZeroIndW]              
                
                DU = numpy.c_[X[notZeroIndU], Y[notZeroIndU], numpy.ones(len(BU))]
                DV = numpy.c_[X[notZeroIndV], Y[notZeroIndV], numpy.ones(len(BV))]
                DW = numpy.c_[X[notZeroIndW], Y[notZeroIndW], numpy.ones(len(BW))]
                
                WU = numpy.sqrt(weightU[notZeroIndU])
                WV = numpy.sqrt(weightV[notZeroIndV])
                WW = numpy.sqrt(weightW[notZeroIndW])
                
                
                DUW = (WU * DU.T).T
                DVW = (WV * DV.T).T
                DWW = (WW * DW.T).T
                
                BUW = BU * WU
                BVW = BV * WV
                BWW = BW * WW

                CU,_,_,_ = scipy.linalg.lstsq(DUW, BUW)    # coefficients
                CV,_,_,_ = scipy.linalg.lstsq(DVW, BVW)    # coefficients
                CW,_,_,_ = scipy.linalg.lstsq(DWW, BWW)
        
                # evaluate it on grid
                plainU[:,:,iIter] = CU[0]*X + CU[1]*Y + CU[2]
                plainV[:,:,iIter] = CV[0]*X + CV[1]*Y + CV[2]
                plainW[:,:,iIter] = CW[0]*X + CW[1]*Y + CW[2]


    elif eddyOrder == 2:
        # best-fit quadratic curve
        
        for iIter in range(plainU.shape[2]):

            # best-fit linear plane
            UFit = UOrg[:, :, iIter, -1].copy()
            VFit = VOrg[:, :, iIter, -1].copy()
            WFit = WOrg[:, :, iIter, -1].copy()
                
            magDataSelected = magData[:,:,iIter].copy()
                
            USTD[1:(UOrg.shape[0]-1),1:(UOrg.shape[1]-1)] = numpy.std(rolling_window(UOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
            VSTD[1:(VOrg.shape[0]-1),1:(VOrg.shape[1]-1)] = numpy.std(rolling_window(VOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
            WSTD[1:(WOrg.shape[0]-1),1:(WOrg.shape[1]-1)] = numpy.std(rolling_window(WOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))


            USTDSelectInd = numpy.where(USTD > (eddyCurrentThreshold/100) * USTD.max(), USTD, 0)
            VSTDSelectInd = numpy.where(VSTD > (eddyCurrentThreshold/100) * VSTD.max(), VSTD, 0)
            WSTDSelectInd = numpy.where(WSTD > (eddyCurrentThreshold/100) * WSTD.max(), WSTD, 0)
    
            
            with numpy.errstate(invalid='ignore'):
                 weightU = magDataSelected / (USTD) ** STDPower
                 weightV = magDataSelected / (VSTD) ** STDPower
                 weightW = magDataSelected / (WSTD) ** STDPower
           
            weightU[numpy.isnan(weightU)] = 0
            weightV[numpy.isnan(weightV)] = 0
            weightW[numpy.isnan(weightW)] = 0
            weightU[numpy.isinf(weightU)] = 0
            weightV[numpy.isinf(weightV)] = 0
            weightW[numpy.isinf(weightW)] = 0
    
                     
            notZeroIndU = numpy.where(numpy.logical_and(weightU > 0, USTDSelectInd))
            notZeroIndV = numpy.where(numpy.logical_and(weightV > 0, VSTDSelectInd))
            notZeroIndW = numpy.where(numpy.logical_and(weightW > 0, WSTDSelectInd))
                
            BU = UFit[notZeroIndU]
            BV = VFit[notZeroIndV]
            BW = WFit[notZeroIndW]
                    
                
            DU = numpy.c_[X[notZeroIndU], Y[notZeroIndU], XY[notZeroIndU], X2[notZeroIndU] , Y2[notZeroIndU] , numpy.ones(len(BU))]
            DV = numpy.c_[X[notZeroIndV], Y[notZeroIndV], XY[notZeroIndV], X2[notZeroIndV] , Y2[notZeroIndV] , numpy.ones(len(BV))]
            DW = numpy.c_[X[notZeroIndW], Y[notZeroIndW], XY[notZeroIndW], X2[notZeroIndW] , Y2[notZeroIndW] , numpy.ones(len(BW))]

            WU = numpy.sqrt(weightU[notZeroIndU])
            WV = numpy.sqrt(weightV[notZeroIndV])
            WW = numpy.sqrt(weightW[notZeroIndW])
                            
            DUW = (WU * DU.T).T
            DVW = (WV * DV.T).T
            DWW = (WW * DW.T).T
                
            BUW = BU * WU
            BVW = BV * WV
            BWW = BW * WW
                
            CU,_,_,_ = scipy.linalg.lstsq(DUW, BUW)    # coefficients
            CV,_,_,_ = scipy.linalg.lstsq(DVW, BVW)    # coefficients
            CW,_,_,_ = scipy.linalg.lstsq(DWW, BWW)
        
            # evaluate it on grid
            plainU[:,:,iIter] = CU[0]*X + CU[1]*Y + CU[2]*XY + CU[3]*X2 + CU[4]*Y2 + CU[5]
            plainV[:,:,iIter] = CV[0]*X + CV[1]*Y + CV[2]*XY + CV[3]*X2 + CV[4]*Y2 + CV[5]
            plainW[:,:,iIter] = CW[0]*X + CW[1]*Y + CW[2]*XY + CW[3]*X2 + CW[4]*Y2 + CW[5]
            
    
    
    if plotEddyPlane:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d') 
        # Plot the surface.
        ax.plot_surface(X, Y, plainU[:,:,plotPlain], linewidth=0, antialiased=False)
        #fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()
        numpy.save("PlaneU", plainU)
        numpy.save("PlaneV", plainV)
        numpy.save("PlaneW", plainW)

   


    for k in range(UOrg.shape[3]):
        
        flowCorrected[:, :, :,0, k] = UOrg[:,:,:,k] - plainU
        flowCorrected[:, :, :,1, k] = VOrg[:,:,:,k] - plainV
        flowCorrected[:, :, :,2, k] = WOrg[:,:,:,k] - plainW

    maskT = weightU * weightV * weightW    
    mask = numpy.where(maskT > 0)
    flowCorrected[mask,:] = 0
    print("end of correction")
    
    return flowCorrected

def randNoise(UOrg, VOrg, WOrg, randThre=25):

    USTD = numpy.zeros((UOrg.shape[0],UOrg.shape[1],UOrg.shape[2]))
    VSTD = numpy.zeros(USTD.shape)
    WSTD = numpy.zeros(USTD.shape)

    

    for kIter in range(UOrg.shape[2]):

        USTD[1:(UOrg.shape[0]-1),1:(UOrg.shape[1]-1), kIter] = numpy.std(rolling_window(UOrg[:,:,kIter,:], (3,3,0), toend=False), axis=(4,3,1))
        VSTD[1:(VOrg.shape[0]-1),1:(VOrg.shape[1]-1), kIter] = numpy.std(rolling_window(VOrg[:,:,kIter,:], (3,3,0), toend=False), axis=(4,3,1))
        WSTD[1:(WOrg.shape[0]-1),1:(WOrg.shape[1]-1), kIter] = numpy.std(rolling_window(WOrg[:,:,kIter,:], (3,3,0), toend=False), axis=(4,3,1))
            
    
    flowCorrected = numpy.zeros([UOrg.shape[0], UOrg.shape[1], UOrg.shape[2],3,UOrg.shape[3]])

    noiseMask = numpy.zeros(USTD.shape)
    noiseMask[(USTD > (randThre*USTD.max()/100)) & (VSTD > (randThre*VSTD.max()/100) ) & (WSTD > (randThre*WSTD.max()/100))] = 1

    UOrg[(USTD > (randThre*USTD.max()/100)) & (VSTD > (randThre*VSTD.max()/100) ) & (WSTD > (randThre*WSTD.max()/100))] = 0
    VOrg[(USTD > (randThre*USTD.max()/100)) & (VSTD > (randThre*VSTD.max()/100) ) & (WSTD > (randThre*WSTD.max()/100))] = 0
    WOrg[(USTD > (randThre*USTD.max()/100)) & (VSTD > (randThre*VSTD.max()/100) ) & (WSTD > (randThre*WSTD.max()/100))] = 0
    
    flowCorrected[:,:,:,0] = UOrg
    flowCorrected[:,:,:,1] = VOrg
    flowCorrected[:,:,:,2] = WOrg


    return flowCorrected



def randNoiseV2(magData, UOrg, VOrg, WOrg, randThre=0.2, plotBool=1, plotPlain=20):


    noiseMask = numpy.zeros(magData.shape)
    noiseThre = (randThre) * (magData.max()- magData.min())
    noiseMask[magData > noiseThre] = 1
    print("mag data max")
    print(magData.max())
    print("mag data min")
    print(magData.min())    
    print("noise threshold")
    print(noiseThre)
    print("noise mask shape")
    print(noiseMask.shape)
    
    print(noiseMask.max())
    PixelSize = numpy.array([0.70, 0.70, 0.4])
    
    saveVTK.saveVTKSeg(noiseMask,False,False, PixelSize, 0, "../")
    
            
    
    return noiseThre

