import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse
from clint.textui import colored


def readGEFlow(args):
          
    print( colored.green("\nLooking for flow data... \n"))
    MagPathStr = args.input
    foldersList = [os.path.join(MagPathStr,o) for o in os.listdir(MagPathStr) if os.path.isdir(os.path.join(MagPathStr,o))]
        
    if not foldersList:
        filesListTEMP = glob.glob(MagPathStr + "/*") 
            
        ds = dicom.read_file(filesListTEMP[0])
        if "GE" in ds.Manufacturer: 
                print("It's GE sequence!")
        else:
                print("We currently can not load files from " + ds.Manufacturer + ".")
    else:
            
                for dirName in foldersList:
                    filesListTEMP = glob.glob(dirName + "/*") 
                    
                    ds = dicom.read_file(filesListTEMP[0])
                    if "GE" in ds.Manufacturer:
                        
                        proceed = True
                        if 100 <= int(ds.SeriesNumber) <= 199:
                            PathFlowDataMAG = dirName
                            print(colored.cyan("Mag folder is: " + PathFlowDataMAG))
                        if 200 <= int(ds.SeriesNumber) <= 299:
                            PathFlowDataRL = dirName
                            print(colored.cyan("R/L folder is: " + PathFlowDataRL))
                            ConstDimsTemp = (int(ds.Rows), int(ds.Columns), int(ds.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
                            dXY = ds.PixelSpacing
                            dZ = ds.SpacingBetweenSlices
                            pixel_spc = (dXY[0],dXY[1],dZ)
                            #print(pixel_spc)
                        if 300 <= int(ds.SeriesNumber) <= 399:
                            PathFlowDataAP = dirName
                            print(colored.cyan("A/P folder is: " + PathFlowDataAP))
                        if 400 <= int(ds.SeriesNumber) <= 499:
                            PathFlowDataSI = dirName
                            print(colored.cyan("S/I folder is: " + PathFlowDataSI))
                            
                    else:
                        proceed = False
                        print(colored.red("FatalError: We currently can not load files from " + ds.Manufacturer + "."))
                        sys.exit()
          
 #   MagPathStr = str(FolderPath)
    PathList=MagPathStr.split("/")
    basePath = MagPathStr.replace(PathList[-1],"")


    if proceed:
      flowData = None 
      for folderNumber in range(0,4):
        sliceLocation = []
        # flow files list
        lstFilesDCM = []
        triggerTime = []
            
        if folderNumber == 0:
            folderPath = PathFlowDataMAG
            print(colored.cyan("Reading the Magnitude files."))
           
        if folderNumber == 1:
            folderPath = PathFlowDataRL
            print(colored.cyan("Reading the flow files (R/L)."))
           
        if folderNumber == 2:
            folderPath = PathFlowDataAP
            print(colored.cyan("Reading the flow files (A/P)."))
            
        if folderNumber == 3:
            folderPath = PathFlowDataSI
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
                ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), int(RefDs.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
                ReadData = numpy.zeros(ConstPixelDims, dtype=numpy.double)

                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)
                    ReadData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation),triggerTimeTemp.index(dsTemp.TriggerTime)]= dsTemp.pixel_array.astype('float')
                    
                magDataTemp = ReadData.mean(3)
                if args.mat:
                    scipy.io.savemat(args.output + "/mag.mat", mdict={'magDataTemp': magDataTemp})

        else:
                if flowData is None:
                    ConstFlowPixelDims = (int(RefDs.Rows), int(RefDs.Columns), int(RefDs.ImagesInAcquisition), 3, int(ds.CardiacNumberOfImages))
                    flowData = numpy.zeros(ConstFlowPixelDims, dtype=numpy.double)
                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)        
                    flowData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation), folderNumber-1,triggerTimeTemp.index(dsTemp.TriggerTime)]= dsTemp.pixel_array.astype('float')
    
                if args.mat:
                    scipy.io.savemat(args.output + "/vel.mat", mdict={'flowData': flowData})
                #print(flowData.shape)


    
    ### The combination of +x +y and -z and permuted x and y is working. Ali Aug24 2017
    UOrg = args.velocitysign[0] * (flowData[:, :, :, args.velocityorder[0]].squeeze())
    VOrg = args.velocitysign[1] * (flowData[:, :, :, args.velocityorder[1]].squeeze())
    WOrg = args.velocitysign[2] * (flowData[:, :, :, args.velocityorder[2]].squeeze())
    
    flowCorrected = numpy.zeros([flowData.shape[0], flowData.shape[1], flowData.shape[2],3,flowData.shape[4]])
    
      
    if args.eddycurrent:
        flowCorrected = eddyCurrentCorrection(UOrg, VOrg, WOrg, args)
        
        
    else:
        flowCorrected[:, :, :,0, :] = UOrg
        flowCorrected[:, :, :,1, :] = VOrg
        flowCorrected[:, :, :,2, :] = WOrg
        
    
    
    
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
    magSize = magDataTemp.shape
    totalNodes = magSize[0] * magSize[1] * magSize[2]
    if args.vtk and args.mat is False:
        print(colore.yellow("We will not save any file since you didnt select your preference. (VTK or MAT)"))
    
    if args.vtk:
        saveVTK(magDataTemp, flowCorrected,pixel_spc, totalNodes, args.output)
    
    if args.mat:
        with open(outPath + "/FlowData.mat", 'wb') as matlabFile:
            scipy.io.savemat(matlabFile, mdict={'velocity': flowData})
            scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
    
    
    
    return RefDs
 #   printReport(outPath, RefDs)
    
