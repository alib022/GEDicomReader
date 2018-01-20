import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse, math
from clint.textui import colored
from vtk.util import numpy_support

def readGETOF(args, PatientDataStruc):
          
    print( colored.green("\nLooking for TOF data... \n"))
    MagPathStr = args.input
    
        
    
    filesListTEMP = glob.glob(MagPathStr + "/*") 
    ds = dicom.read_file(filesListTEMP[0])
    ConstDimsTemp = (int(ds.Rows), int(ds.Columns), int(ds.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
 #   ConstDimsTemp = (int(ds.Rows), int(ds.Columns), int(math.ceil(len(filesListTEMP)/ int(ds.CardiacNumberOfImages))), int(ds.CardiacNumberOfImages))

    dXY = ds.PixelSpacing
    dZ = ds.SpacingBetweenSlices
    pixel_spc = (dXY[0],dXY[1],dZ)
    if "GE" in ds.Manufacturer:
        proceed = True
    else:
        proceed = False
        print(colored.red("FatalError: We currently can not load files from " + ds.Manufacturer + "."))
        sys.exit()
          
 #   MagPathStr = str(FolderPath)
    PathList=MagPathStr.split("/")
    basePath = MagPathStr.replace(PathList[-1],"")


    if proceed:
        sliceLocation = []
        # flow files list
        lstFilesDCM = []
        
         
            ################## Reading time of flight files
            # listing magnitude files
        for dirName, subdirList, fileList in os.walk( args.input + "/"):
            for filename in fileList:
            # if ".dcm" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))
                ds = dicom.read_file(lstFilesDCM[-1])
                sliceLocation.append(ds.SliceLocation)
                
                
            # Get ref file
        RefDs = dicom.read_file(lstFilesDCM[0])
            
            
        #triggerTimeTemp = sorted(set(triggerTime), key=float)
        #print(triggerTimeTemp)
        #sliceLocationTemp = set(sliceLocation)       
        sliceLocationTemp = sorted(set(sliceLocation), key=float)
        #print(len(sliceLocationTemp))


       

            
        ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns),math.ceil(len(filesListTEMP)/ int(ds.CardiacNumberOfImages)))
        ReadData = numpy.zeros(ConstPixelDims, dtype=numpy.double)

        for iFile in lstFilesDCM:
            dsTemp = dicom.read_file (iFile)
            #print(dsTemp.TriggerTime)
            ReadData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation)]= dsTemp.pixel_array.astype('float')
                    
    #print(ReadData.shape)
    magDataTemp = ReadData
    if args.mat:
       scipy.io.savemat(args.output + "/TOF.mat", mdict={'TOF': magDataTemp})
            
        
        
        
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
    magSize = magDataTemp.shape
    totalNodes = magSize[0] * magSize[1] * magSize[2]

    if (args.vtk == False and args.mat == False):
        print(colored.yellow("We will not save any file since you didnt select your preference! (VTK or MAT)"))
        
    
    if args.vtk:
        
        saveVTKSeg(magDataTemp, False,True, pixel_spc, totalNodes, args.output)
        
    if args.mat:
        with open(args.output + "/TOFData.mat", 'wb') as matlabFile:
            scipy.io.savemat(matlabFile, mdict={'TOF': magDataTemp})
            
        
    
    
 #   return RefDs
 #   printReport(outPath, RefDs)
def readGEcMRA(args, PatientDataStruc):
          
    print( colored.green("\nLooking for cMRA data... \n"))
    MagPathStr = args.input
    
        
    
    filesListTEMP = glob.glob(MagPathStr + "/*") 
    ds = dicom.read_file(filesListTEMP[0])
    ConstDimsTemp = (int(ds.Rows), int(ds.Columns), int(ds.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
    dXY = ds.PixelSpacing
    dZ = ds.SpacingBetweenSlices
    pixel_spc = (dXY[0],dXY[1],dZ)
    if "GE" in ds.Manufacturer:
        proceed = True
    else:
        proceed = False
        print(colored.red("FatalError: We currently can not load files from " + ds.Manufacturer + "."))
        sys.exit()
          
 #   MagPathStr = str(FolderPath)
    PathList=MagPathStr.split("/")
    basePath = MagPathStr.replace(PathList[-1],"")


    if proceed:
        sliceLocation = []
        # flow files list
        lstFilesDCM = []
        triggerTime = []
         
            ################## Reading time of flight files
            # listing magnitude files
        for dirName, subdirList, fileList in os.walk( args.input + "/"):
            for filename in fileList:
            # if ".dcm" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))
                ds = dicom.read_file(lstFilesDCM[-1])
                sliceLocation.append(ds.SliceLocation)
                triggerTime.append(ds.TriggerTime)
                
            # Get ref file
        RefDs = dicom.read_file(lstFilesDCM[0])
            
            
        triggerTimeTemp = sorted(set(triggerTime), key=float)
        #print(triggerTimeTemp)
        #sliceLocationTemp = set(sliceLocation)       
        sliceLocationTemp = sorted(set(sliceLocation), key=float)
        #print(len(sliceLocationTemp))


       

            
        ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), int(len(sliceLocationTemp)), int(len(triggerTimeTemp)))
        print(ConstPixelDims)
        ReadData = numpy.zeros(ConstPixelDims, dtype=numpy.double)

        for iFile in lstFilesDCM:
            dsTemp = dicom.read_file (iFile)
            #print(dsTemp.TriggerTime)
            ReadData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation), triggerTimeTemp.index(dsTemp.TriggerTime)]= dsTemp.pixel_array.astype('float')
                    
    #print(ReadData.shape)
    magDataTemp = ReadData.mean(3)
    if args.mat:
       scipy.io.savemat(args.output + "/cMRA.mat", mdict={'cMRA': magDataTemp})
            
        
        
        
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
    magSize = magDataTemp.shape
    totalNodes = magSize[0] * magSize[1] * magSize[2]

    if (args.vtk == False and args.mat == False):
        print(colored.yellow("We will not save any file since you didnt select your preference! (VTK or MAT)"))
        
    
    if args.vtk:
        
        saveVTKSeg(magDataTemp,True, False, pixel_spc, totalNodes, args.output)
        
    if args.mat:
        with open(args.output + "/cMRAData.mat", 'wb') as matlabFile:
            scipy.io.savemat(matlabFile, mdict={'cMRA': magDataTemp})
            
        
    
    
 #   return RefDs
 #   printReport(outPath, RefDs)

def readGEFlow(inputFlags, PatientDataStruc):
          
    print( colored.green("\nLooking for flow data... \n"))
    print(colored.cyan("Mag folder is: " + PatientDataStruc.MagPath))
    print(colored.cyan("R/L folder is: " + PatientDataStruc.FlowPathRL))
    print(colored.cyan("A/P folder is: " + PatientDataStruc.FlowPathAP))
    print(colored.cyan("S/I folder is: " + PatientDataStruc.FlowPathSI))
                            

    for folderNumber in range(0,4):
        sliceLocation = []
        # flow files list
        lstFilesDCM = []
        triggerTime = []
            
        if folderNumber == 0:
            folderPath = PatientDataStruc.MagPath
            print(colored.cyan("Reading the Magnitude files."))
           
        if folderNumber == 1:
            if inputFlags.segmentation:
                break
            folderPath = PatientDataStruc.FlowPathRL
            print(colored.cyan("Reading the flow files (R/L)."))
           
        if folderNumber == 2:
            if inputFlags.segmentation:
                break
            folderPath = PatientDataStruc.FlowPathAP
            print(colored.cyan("Reading the flow files (A/P)."))
            
        if folderNumber == 3:
            if inputFlags.segmentation:
                break
            folderPath = PatientDataStruc.FlowPathSI
            print(colored.cyan("Reading the flow files (S/I)."))
           
            
          
            ################## Reading time of flight files
            # listing magnitude files
        for dirName, subdirList, fileList in os.walk( folderPath + "/"):
            for filename in fileList:
            # if ".dcm" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))
                ds = dicom.read_file(lstFilesDCM[-1])
                sliceLocation.append(ds.SliceLocation)
                triggerTime.append(ds.TriggerTime)
                
            # Get ref file
        RefDs = dicom.read_file(lstFilesDCM[0])
            
            
        triggerTimeTemp = sorted(set(triggerTime), key=float)
        #print(triggerTimeTemp)
        #sliceLocationTemp = set(sliceLocation)       
        sliceLocationTemp = sorted(set(sliceLocation), key=float)
        #print(sliceLocationTemp)

       

            
        if folderNumber == 0:
                #ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns),66, int(ds.CardiacNumberOfImages))
                #print(ConstPixelDims)
                ReadData = numpy.zeros(PatientDataStruc.MagVecSize, dtype=numpy.double)

                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)
                    ReadData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation),triggerTimeTemp.index(dsTemp.TriggerTime)]= dsTemp.pixel_array.astype('float')
                    
                magDataTemp = ReadData.mean(3)
                if inputFlags.mat:
                    scipy.io.savemat(inputFlags.output + "/mag.mat", mdict={'magDataTemp': magDataTemp})
 #               numpy.save(args.output +"/mag", magDataTemp) 
                
        else:
                if folderNumber == 1:
 #                   ConstFlowPixelDims = (int(RefDs.Rows), int(RefDs.Columns), 66, 3, int(ds.CardiacNumberOfImages))
                    flowData = numpy.zeros(PatientDataStruc.FlowVecSize, dtype=numpy.double)
                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)        
                    flowData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation), folderNumber-1,triggerTimeTemp.index(dsTemp.TriggerTime)]= dsTemp.pixel_array.astype('float')
    
                if inputFlags.mat:
                    scipy.io.savemat(inputFlags.output + "/vel.mat", mdict={'flowData': flowData})
                #print(flowData.shape)


    if inputFlags.segmentation is False:
        ### The combination of +x +y and -z and permuted x and y is working. Ali Aug24 2017
        UOrg = inputFlags.velocitysign[0] * (flowData[:, :, :, args.velocityorder[0]].squeeze())
        VOrg = inputFlags.velocitysign[1] * (flowData[:, :, :, args.velocityorder[1]].squeeze())
        WOrg = inputFlags.velocitysign[2] * (flowData[:, :, :, args.velocityorder[2]].squeeze())
        
        flowCorrected = numpy.zeros([flowData.shape[0], flowData.shape[1], flowData.shape[2],3,flowData.shape[4]])
        
          
        if inputFlags.eddycurrent:
            flowCorrected = eddyCurrentCorrection(UOrg, VOrg, WOrg, inputFlags.randomnoise, inputFlags.eddythreshold, inputFlags.eddyplane)

        elif inputFlags.randomnoise is not None:
            flowCorrected = randNoise(inputFlags, UOrg, VOrg, WOrg, inputFlags.randomnoise/100, 0, 1)
            
            
        else:
            flowCorrected[:, :, :,0, :] = UOrg
            flowCorrected[:, :, :,1, :] = VOrg
            flowCorrected[:, :, :,2, :] = WOrg
            
        
        
        
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
    magSize = magDataTemp.shape
    totalNodes = magSize[0] * magSize[1] * magSize[2]

    if (inputFlags.vtk == False and inputFlags.mat == False):
        print(colored.yellow("We will ONLY save in npy format, since you didnt select your preference! (VTK or MAT)"))
        
    if not inputFlags.segmentation:
        numpy.save(inputFlags.output +"/FlowData", flowCorrected)    
    if inputFlags.vtk:
        if inputFlags.segmentation:
            
            saveVTKSeg(magDataTemp,False,False, PatientDataStruc.PixelSize, totalNodes, inputFlags.output)
        else:
            
            saveVTK(magDataTemp, flowCorrected,  PatientDataStruc.PixelSize, totalNodes, inputFlags.output)
    
    
    if inputFlags.mat:
        if inputFlags.segmentation:
            with open(inputFlags.output + "/FlowData.mat", 'wb') as matlabFile:
                scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
            
        else:
            with open(inputFlags.output + "/FlowData.mat", 'wb') as matlabFile:
                scipy.io.savemat(matlabFile, mdict={'velocity': flowData})
                scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
    
    
    
 #   return RefDs
 #   printReport(outPath, RefDs)

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
 
def saveVTK(magDataTemp, flowCorrected,pixel_spc, totalNodes, outPath):
        
        for timeIter in range(0,flowCorrected.shape[4]):
        
            #   magData = numpy_support.numpy_to_vtk(num_array=magDataTemp.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
            #   flowDataInVTK = numpy_support.numpy_to_vtk(num_array=flowData[:,:,:,:,timeIter].ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        
            magDataOut = numpy.reshape(magDataTemp, (totalNodes, 1), order='F')
        
            flowDataOut = numpy.reshape(flowCorrected[:,:,:,:,timeIter], (totalNodes, 3), order='F')
        
        #        #########################
        #        # Convert the VTK array to vtkImageData
        
            img_vtk = vtk.vtkImageData()
            img_vtk.SetDimensions(magDataTemp.shape)
            img_vtk.SetSpacing(pixel_spc)
            img_vtk.SetOrigin(0,0,0)
    
        
            magArray = vtk.vtkDoubleArray()
            magArray.SetNumberOfComponents(1)
            magArray.SetNumberOfTuples(img_vtk.GetNumberOfPoints())
            magArray.SetName("Magnitude")
    
            for i in range(0,totalNodes):
                magArray.SetValue(i,magDataOut[i])
        
            flowArray = vtk.vtkDoubleArray()
            flowArray.SetNumberOfComponents(3)
            flowArray.SetNumberOfTuples(img_vtk.GetNumberOfPoints())
            flowArray.SetName("Velocity")
            for i in range(0,totalNodes):
                flowArray.SetTuple3(i,flowDataOut[i,0], flowDataOut[i,1], flowDataOut[i,2])
        
            img_vtk.GetPointData().AddArray(magArray)
            img_vtk.GetPointData().AddArray(flowArray)
        
        
            writer = vtk.vtkXMLImageDataWriter()
            writer.SetFileName(outPath + '/FlowData_%d.vti' % timeIter)
            writer.SetInputData(img_vtk)
            ## This is set so we can see the data in a text editor.
            writer.SetDataModeToAscii()
            writer.Write()

        

def saveVTKSeg(magDataTemp, cMRA,TOF, pixel_spc, totalNodes, outPath):
        
        
            #   magData = numpy_support.numpy_to_vtk(num_array=magDataTemp.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
            #   flowDataInVTK = numpy_support.numpy_to_vtk(num_array=flowData[:,:,:,:,timeIter].ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        
 #           magDataOut = numpy.reshape(magDataTemp, (totalNodes, 1), order='F')
            MagVTK = numpy_support.numpy_to_vtk(num_array=magDataTemp.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
                    
            #        #########################
            #        # Convert the VTK array to vtkImageData
            img_vtk = vtk.vtkImageData()
            img_vtk.SetDimensions(magDataTemp.shape)
            img_vtk.AllocateScalars(vtk.VTK_FLOAT,1)
            img_vtk.SetSpacing(pixel_spc)
            img_vtk.SetOrigin(0,0,0)
            img_vtk.GetPointData().SetScalars(MagVTK)        
        
        #        #########################
        #        # Convert the VTK array to vtkImageData
        
 #           img_vtk = vtk.vtkImageData()
#            img_vtk.SetDimensions(magDataTemp.shape)
#            img_vtk.SetSpacing(pixel_spc)
#            img_vtk.SetOrigin(0,0,0)

#            dims = imageData.GetDimensions()
 
            # Fill every entry of the image data with "2.0"
#            for z in range(dims[2]):
#                for y in range(dims[1]):
#                    for x in range(dims[0]):
#                        img_vtk.SetScalarComponentFromDouble(x, y, z, 0, 2.0)
#            for i in range(0,totalNodes):
#                 img_vtk.SetValue(i,magDataOut[i])
        

        
 #           magArray = vtk.vtkDoubleArray()
#            magArray.SetNumberOfComponents(1)
#            magArray.SetNumberOfTuples(img_vtk.GetNumberOfPoints())
#            magArray.SetName("Magnitude")
    
                    
#            img_vtk.GetPointData().AddArray(magArray)
                
            writer = vtk.vtkXMLImageDataWriter()

            if cMRA:
                writer.SetFileName(outPath + '/cMRAData.vti')
            elif TOF:
                writer.SetFileName(outPath + '/TOFData.vti')                
            else:
                writer.SetFileName(outPath + '/MagData.vti')

            writer.SetInputData(img_vtk)
            ## This is set so we can see the data in a text editor.
            writer.SetDataModeToAscii()
            writer.Write()

            return 0

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

