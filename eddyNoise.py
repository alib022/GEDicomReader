import scipy.io, numpy, saveVTK, warnings
#from clint.textui import colored
#from vtk.util import numpy_support
#from scipy.ndimage.filters import uniform_filter
from rolling_window import rolling_window



import matplotlib.pyplot as plt
#from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


warnings.simplefilter(action='ignore', category=FutureWarning)

def eddyCurrentCorrection(UOrg, VOrg, WOrg, magData, STDPower=2,  eddyCurrentThreshold=15, eddyOrder=2, plotBool=0, verbous=1, plotEddyPlane=1, plotPlain=20):


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
                
                UFit = UFit.ravel()
                VFit = VFit.ravel()
                WFit = WFit.ravel()
                
                magDataSelected = magData[:,:,iIter].copy()
                
                USTD[1:(UOrg.shape[0]-1),1:(UOrg.shape[1]-1)] = numpy.std(rolling_window(UOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
                VSTD[1:(VOrg.shape[0]-1),1:(VOrg.shape[1]-1)] = numpy.std(rolling_window(VOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
                WSTD[1:(WOrg.shape[0]-1),1:(WOrg.shape[1]-1)] = numpy.std(rolling_window(WOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))

                USTDSelectInd = numpy.where(USTD > (eddyCurrentThreshold/100) * USTD.max())
                VSTDSelectInd = numpy.where(VSTD > (eddyCurrentThreshold/100) * VSTD.max())
                WSTDSelectInd = numpy.where(WSTD > (eddyCurrentThreshold/100) * WSTD.max())
    
                USTDSelected = USTD[USTDSelectInd].copy()
                VSTDSelected = VSTD[VSTDSelectInd].copy()
                WSTDSelected = WSTD[WSTDSelectInd].copy()
            
                with numpy.errstate(invalid='ignore'):
                      weightU = magDataSelected[USTDSelectInd] / (USTDSelected) ** STDPower
                      weightV = magDataSelected[VSTDSelectInd] / (VSTDSelected) ** STDPower
                      weightW = magDataSelected[WSTDSelectInd] / (WSTDSelected) ** STDPower
            
                weightU[numpy.isnan(weightU)] = 0
                weightV[numpy.isnan(weightV)] = 0
                weightW[numpy.isnan(weightW)] = 0
                weightU[numpy.isinf(weightU)] = 0
                weightV[numpy.isinf(weightV)] = 0
                weightW[numpy.isinf(weightW)] = 0
    
                     
                notZeroIndU = numpy.where(weightU > 0)
                notZeroIndV = numpy.where(weightV > 0)
                notZeroIndW = numpy.where(weightW > 0)
                
                BU = UFit[notZeroIndU]
                BV = VFit[notZeroIndV]
                BW = WFit[notZeroIndW]
                    
                Xravel = X.ravel()
                Yravel = Y.ravel()
                
                
                DU = numpy.c_[Xravel[notZeroIndU], Yravel[notZeroIndU], numpy.ones(len(BU))]
                DV = numpy.c_[Xravel[notZeroIndV], Yravel[notZeroIndV], numpy.ones(len(BV))]
                DW = numpy.c_[Xravel[notZeroIndW], Yravel[notZeroIndW], numpy.ones(len(BW))]
                
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
                
            UFit = UFit.ravel()
            VFit = VFit.ravel()
            WFit = WFit.ravel()
                
            magDataSelected = magData[:,:,iIter].copy()
                
            USTD[1:(UOrg.shape[0]-1),1:(UOrg.shape[1]-1)] = numpy.std(rolling_window(UOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
            VSTD[1:(VOrg.shape[0]-1),1:(VOrg.shape[1]-1)] = numpy.std(rolling_window(VOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))
            WSTD[1:(WOrg.shape[0]-1),1:(WOrg.shape[1]-1)] = numpy.std(rolling_window(WOrg[:,:,iIter,:], (3,3,0), toend=False), axis=(4,3,1))

            USTDSelectInd = numpy.where(USTD > (eddyCurrentThreshold/100) * USTD.max())
            VSTDSelectInd = numpy.where(VSTD > (eddyCurrentThreshold/100) * VSTD.max())
            WSTDSelectInd = numpy.where(WSTD > (eddyCurrentThreshold/100) * WSTD.max())
    
            USTDSelected = USTD[USTDSelectInd].copy()
            VSTDSelected = VSTD[VSTDSelectInd].copy()
            WSTDSelected = WSTD[WSTDSelectInd].copy()
            
            with numpy.errstate(invalid='ignore'):
                 weightU = magDataSelected[USTDSelectInd] / (USTDSelected) ** STDPower
                 weightV = magDataSelected[VSTDSelectInd] / (VSTDSelected) ** STDPower
                 weightW = magDataSelected[WSTDSelectInd] / (WSTDSelected) ** STDPower
           
            weightU[numpy.isnan(weightU)] = 0
            weightV[numpy.isnan(weightV)] = 0
            weightW[numpy.isnan(weightW)] = 0
            weightU[numpy.isinf(weightU)] = 0
            weightV[numpy.isinf(weightV)] = 0
            weightW[numpy.isinf(weightW)] = 0
    
                     
            notZeroIndU = numpy.where(weightU > 0)
            notZeroIndV = numpy.where(weightV > 0)
            notZeroIndW = numpy.where(weightW > 0)
                
            BU = UFit[notZeroIndU]
            BV = VFit[notZeroIndV]
            BW = WFit[notZeroIndW]
                    
            Xravel = X.ravel()
            Yravel = Y.ravel()
            XYravel = XY.ravel()
            X2ravel = X2.ravel()
            Y2ravel = Y2.ravel()
                
            DU = numpy.c_[Xravel[notZeroIndU], Yravel[notZeroIndU], XYravel[notZeroIndU], X2ravel[notZeroIndU] , Y2ravel[notZeroIndU] , numpy.ones(len(BU))]
            DV = numpy.c_[Xravel[notZeroIndV], Yravel[notZeroIndV], XYravel[notZeroIndV], X2ravel[notZeroIndV] , Y2ravel[notZeroIndV] , numpy.ones(len(BV))]
            DW = numpy.c_[Xravel[notZeroIndW], Yravel[notZeroIndW], XYravel[notZeroIndW], X2ravel[notZeroIndW] , Y2ravel[notZeroIndW] , numpy.ones(len(BW))]
                
            CU,_,_,_ = scipy.linalg.lstsq(DU, BU)    # coefficients
            CV,_,_,_ = scipy.linalg.lstsq(DV, BV)    # coefficients
            CW,_,_,_ = scipy.linalg.lstsq(DW, BW)
        
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

   


    for k in range(UOrg.shape[3]):
        
        flowCorrected[:, :, :,0, k] = UOrg[:,:,:,k] - plainU
        flowCorrected[:, :, :,1, k] = VOrg[:,:,:,k] - plainV
        flowCorrected[:, :, :,2, k] = WOrg[:,:,:,k] - plainW
    
    
    return flowCorrected

def randNoise(UOrg, VOrg, WOrg, randThre=25, plotBool=1, plotPlain=20):

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
    
    PixelSize = numpy.array([0.70, 0.70, 0.4])
    saveVTK.saveVTKSeg(noiseMask,False,False, PixelSize, 0, "../")
    
    with open("randNoisestatictissue.mat", 'wb') as matlabFile:
        scipy.io.savemat(matlabFile, mdict={'UOrg': UOrg})
        scipy.io.savemat(matlabFile, mdict={'VOrg': VOrg})
        scipy.io.savemat(matlabFile, mdict={'WOrg': WOrg})



    flowCorrected[:,:,:,0] = UOrg
    flowCorrected[:,:,:,1] = VOrg
    flowCorrected[:,:,:,2] = WOrg

    
    if plotBool:
        vmax = numpy.max([UOrg[:,:,plotPlain,1].max(), VOrg[:,:,plotPlain,1].max(), WOrg[:,:,plotPlain,1].max()])
        vmin = numpy.min([UOrg[:,:,plotPlain,1].min(), VOrg[:,:,plotPlain,1].min(), WOrg[:,:,plotPlain,1].min()])
        vmax = numpy.max([vmax, -vmin])
        vmin = -vmax

        # plot with various axes scales
        plt.figure(1)

       
        plt.subplot(131)
        plt.imshow(UOrg[:,:,plotPlain,1], vmin=vmin, vmax=vmax, cmap='seismic')
        plt.title('U Org')


       
        plt.subplot(132)
        plt.imshow(VOrg[:,:,plotPlain, 1], vmin=vmin, vmax=vmax, cmap='seismic')
        plt.title('V Org')

        plt.subplot(133)
        plt.imshow(WOrg[:,:,plotPlain, 1], vmin=vmin, vmax=vmax, cmap='seismic')
        plt.title('W Org')

        plt.show()
    

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

