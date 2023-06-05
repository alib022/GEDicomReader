import pydicom,os, glob, hdf5storage, numpy, vtk, sys, datetime, argparse, math, saveVTK
from clint.textui import colored
from vtk.util import numpy_support

def readGETOF(args, PatientDataStruc):
          
    print( colored.green("\nLooking for TOF data... \n"))
    MagPathStr = args.input
    
        
    
    filesListTEMP = glob.glob(MagPathStr + "/*") 
    ds = pydicom.read_file(filesListTEMP[0])
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
                ds = pydicom.read_file(lstFilesDCM[-1])
                sliceLocation.append(ds.SliceLocation)
                
                
            # Get ref file
        RefDs = pydicom.read_file(lstFilesDCM[0])
            
            
        #triggerTimeTemp = sorted(set(triggerTime), key=float)
        #print(triggerTimeTemp)
        #sliceLocationTemp = set(sliceLocation)       
        sliceLocationTemp = sorted(set(sliceLocation), key=float)
        #print(len(sliceLocationTemp))


       

        
        if int(RefDs.CardiacNumberOfImages) == 0:
            ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(filesListTEMP))
        else:
            ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns),math.ceil(len(filesListTEMP)/ int(RefDs.CardiacNumberOfImages)))
        
        ReadData = numpy.zeros(ConstPixelDims, dtype=numpy.double)

        for iFile in lstFilesDCM:
            dsTemp = pydicom.read_file (iFile)
            #print(dsTemp.TriggerTime)
            ReadData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation)]= dsTemp.pixel_array.astype('float')
                    
    #print(ReadData.shape)
    magDataTemp = ReadData
    if args.mat:
       hdf5storage.savemat(args.output + "/TOF.mat", {'TOF': magDataTemp}, format='7.3', oned_as='column', store_python_metadata=True)     
        
        
        
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
    magSize = magDataTemp.shape
    totalNodes = magSize[0] * magSize[1] * magSize[2]

    if (args.vtk == False and args.mat == False):
        print(colored.yellow("We will not save any file since you didnt select your preference! (VTK or MAT)"))
        
    
    if args.vtk:
        
        saveVTK.saveVTKSeg(magDataTemp, False,True, pixel_spc, totalNodes, args.output)
        
    if args.mat:
        hdf5storage.savemat(args.output + "/TOFData.mat", {'TOF': magDataTemp}, format='7.3', oned_as='column', store_python_metadata=True)
        
    
    
 #   return RefDs
 #   printReport(outPath, RefDs)
def readGEcMRA(args, PatientDataStruc):
          
    print( colored.green("\nLooking for cMRA data... \n"))
    MagPathStr = args.input
    
        
    
    filesListTEMP = glob.glob(MagPathStr + "/*") 
    ds = pydicom.read_file(filesListTEMP[0])
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
                ds = pydicom.read_file(lstFilesDCM[-1])
                sliceLocation.append(ds.SliceLocation)
                triggerTime.append(ds.TriggerTime)
                
            # Get ref file
        RefDs = pydicom.read_file(lstFilesDCM[0])
            
            
        triggerTimeTemp = sorted(set(triggerTime), key=float)
        #print(triggerTimeTemp)
        #sliceLocationTemp = set(sliceLocation)       
        sliceLocationTemp = sorted(set(sliceLocation), key=float)
        #print(len(sliceLocationTemp))


       

            
        ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), int(len(sliceLocationTemp)), int(len(triggerTimeTemp)))
        print(ConstPixelDims)
        ReadData = numpy.zeros(ConstPixelDims, dtype=numpy.double)

        for iFile in lstFilesDCM:
            dsTemp = pydicom.read_file (iFile)
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
        hdf5storage.savemat(args.output + "/cMRAData.mat", {'cMRA': magDataTemp}, format='7.3', oned_as='column', store_python_metadata=True)   
        
    
    
 #   return RefDs
 #   printReport(outPath, RefDs)
