import scipy.io, numpy, saveVTK
from clint.textui import colored
#from vtk.util import numpy_support
#from scipy.ndimage.filters import uniform_filter
from rolling_window import rolling_window



import matplotlib.pyplot as plt
#from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def eddyCurrentCorrection(UOrg, VOrg, WOrg, eddyCurrentThreshold=8, eddyOrder=2, plotBool=0, verbous=0, plotEddyPlane=1, plotPlain=20):


    USTD = numpy.zeros((UOrg.shape[0],UOrg.shape[1],UOrg.shape[2]))
    VSTD = numpy.zeros(USTD.shape)
    WSTD = numpy.zeros(USTD.shape)
    

    for kIter in range(UOrg.shape[2]):

        USTD[1:(UOrg.shape[0]-1),1:(UOrg.shape[1]-1), kIter] = numpy.std(rolling_window(UOrg[:,:,kIter,:], (3,3,0), toend=False), axis=(4,3,1))
        VSTD[1:(VOrg.shape[0]-1),1:(VOrg.shape[1]-1), kIter] = numpy.std(rolling_window(VOrg[:,:,kIter,:], (3,3,0), toend=False), axis=(4,3,1))
        WSTD[1:(WOrg.shape[0]-1),1:(WOrg.shape[1]-1), kIter] = numpy.std(rolling_window(WOrg[:,:,kIter,:], (3,3,0), toend=False), axis=(4,3,1))

    
    if verbous:    
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


    

    staticTissueU = UOrg[:,:,:,-1].copy()
    staticTissueV = VOrg[:,:,:,-1].copy()
    staticTissueW = WOrg[:,:,:,-1].copy()

    staticTissueU[(USTD > (eddyCurrentThreshold*USTD.max()/100)) & (VSTD > (eddyCurrentThreshold*VSTD.max()/100) ) & (WSTD > (eddyCurrentThreshold*WSTD.max()/100))] = 0
    staticTissueV[(USTD > (eddyCurrentThreshold*USTD.max()/100)) & (VSTD > (eddyCurrentThreshold*VSTD.max()/100) ) & (WSTD > (eddyCurrentThreshold*WSTD.max()/100))] = 0
    staticTissueW[(USTD > (eddyCurrentThreshold*USTD.max()/100)) & (VSTD > (eddyCurrentThreshold*VSTD.max()/100) ) & (WSTD > (eddyCurrentThreshold*WSTD.max()/100))] = 0
    
    with open("eddyNoisestatictissue.mat", 'wb') as matlabFile:
        scipy.io.savemat(matlabFile, mdict={'staticTissueU': staticTissueU})
        scipy.io.savemat(matlabFile, mdict={'staticTissueV': staticTissueV})
        scipy.io.savemat(matlabFile, mdict={'staticTissueW': staticTissueW})

    
    
    if plotBool:

        vmax = numpy.max([UOrg[:,:,plotPlain,1].max(), staticTissueU[:,:,plotPlain].max()])
        vmin = numpy.min([UOrg[:,:,plotPlain,1].min(), staticTissueU[:,:,plotPlain].min()])
        vmax = numpy.max([vmax, -vmin])
        vmin = -vmax
    
        #for i in range(1,3):
        # plot with various axes scales
        plt.figure(plotPlain)

       
        plt.subplot(121)
        plt.imshow(UOrg[:,:,plotPlain,1], vmin=vmin, vmax=vmax, cmap='seismic')
        plt.title('Org data')


           
        plt.subplot(122)
        plt.imshow(staticTissueU[:,:,plotPlain], vmin=vmin, vmax=vmax, cmap='seismic')
        plt.title('static tissue')

        plt.show()
            

    del USTD, VSTD, WSTD

    
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
        # best-fit linear plane
        for iIter in range(UOrg.shape[2]):
            
                BUInd = staticTissueU[:,:,iIter].copy()
                BVInd = staticTissueV[:,:,iIter].copy()
                BWInd = staticTissueW[:,:,iIter].copy()

                notZeroIndU = numpy.nonzero(BUInd)
                notZeroIndV = numpy.nonzero(BVInd)
                notZeroIndW = numpy.nonzero(BWInd)

                BU = BUInd[notZeroIndU].ravel()
                BV = BVInd[notZeroIndV].ravel()
                BW = BWInd[notZeroIndW].ravel()
                
                print("BU shape")
                print(BU.shape)
                print("BV shape")
                print(BV.shape)
                print("BW shape")
                print(BW.shape)
                
                
                DU = numpy.c_[X[notZeroIndU].ravel(), Y[notZeroIndU].ravel(), numpy.ones(len(BU))]
                DV = numpy.c_[X[notZeroIndV].ravel(), Y[notZeroIndV].ravel(), numpy.ones(len(BV))]
                DW = numpy.c_[X[notZeroIndW].ravel(), Y[notZeroIndW].ravel(), numpy.ones(len(BW))]
            
                print("DU shape")
                print(DU.shape)
                print("DV shape")
                print(DV.shape)
                print("DW shape")
                print(DW.shape)
                
                CU,_,_,_ = scipy.linalg.lstsq(DU, BU)    # coefficients
                CV,_,_,_ = scipy.linalg.lstsq(DV, BV)    # coefficients
                CW,_,_,_ = scipy.linalg.lstsq(DW, BW)
        
                # evaluate it on grid
                plainU[:,:,iIter] = CU[0]*X + CU[1]*Y + CU[2]
                plainV[:,:,iIter] = CV[0]*X + CV[1]*Y + CV[2]
                plainW[:,:,iIter] = CW[0]*X + CW[1]*Y + CW[2]


    elif eddyOrder == 2:
        # best-fit quadratic curve
        
        for iIter in range(plainU.shape[2]):

            BUInd = staticTissueU[:,:,iIter].copy()
            BVInd = staticTissueV[:,:,iIter].copy()
            BWInd = staticTissueW[:,:,iIter].copy()

            notZeroIndU = numpy.nonzero(BUInd)
            notZeroIndV = numpy.nonzero(BVInd)
            notZeroIndW = numpy.nonzero(BWInd)

            BU = BUInd[notZeroIndU].ravel()
            BV = BVInd[notZeroIndV].ravel()
            BW = BWInd[notZeroIndW].ravel()
            
            print("BU shape")
            print(BU.shape)
            print("BV shape")
            print(BV.shape)
            print("BW shape")
            print(BW.shape)
            
                
            DU = numpy.c_[X[notZeroIndU].ravel(), Y[notZeroIndU].ravel(), XY[notZeroIndU].ravel(), X2[notZeroIndU].ravel() , Y2[notZeroIndU].ravel() , numpy.ones(len(BU))]
            DV = numpy.c_[X[notZeroIndV].ravel(), Y[notZeroIndV].ravel(), XY[notZeroIndV].ravel(), X2[notZeroIndV].ravel() , Y2[notZeroIndV].ravel() , numpy.ones(len(BV))]
            DW = numpy.c_[X[notZeroIndW].ravel(), Y[notZeroIndW].ravel(), XY[notZeroIndW].ravel(), X2[notZeroIndW].ravel() , Y2[notZeroIndW].ravel() , numpy.ones(len(BW))]

            print("DU shape")
            print(DU.shape)
            print("DV shape")
            print(DV.shape)
            print("DW shape")
            print(DW.shape)
                
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

