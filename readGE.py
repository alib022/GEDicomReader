import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse
from clint.textui import colored
from readGEFlow import readGEFlow
''' This function reads GE Flow data '''



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
    
    flowCorrected[noiseMask] = 0
    print("After Eddy correction")
    
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

def readPatientInfo(FolderPath):
    
    MagPathStr = str(FolderPath)
    foldersList = [os.path.join(MagPathStr,o) for o in os.listdir(MagPathStr) if os.path.isdir(os.path.join(MagPathStr,o))]
        
    if not foldersList:
        filesListTEMP = glob.glob(MagPathStr + "/*") 
            
        ds = dicom.read_file(filesListTEMP[0])
        if "GE" in ds.Manufacturer: 
                print("It's GE sequence!")
        else:
                print("We currently can not load files from " + ds.Manufacturer + ".")
                sys.exit()
    else:
            
                for dirName in foldersList:
                    filesListTEMP = glob.glob(dirName + "/*") 
                    
                    ds = dicom.read_file(filesListTEMP[0])
                    if "GE" in ds.Manufacturer:
                        proceed = True
                        if 100 <= int(ds.SeriesNumber) <= 199:
                            PathFlowDataMAG = dirName
                            
                        if 200 <= int(ds.SeriesNumber) <= 299:
                            PathFlowDataRL = dirName
                            
                            ConstDimsTemp = (int(ds.Rows), int(ds.Columns), int(ds.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
                            dXY = ds.PixelSpacing
                            dZ = ds.SpacingBetweenSlices
                            pixel_spc = (dXY[0],dXY[1],dZ)
                            #print(pixel_spc)
                        if 300 <= int(ds.SeriesNumber) <= 399:
                            PathFlowDataAP = dirName
                            
                        if 400 <= int(ds.SeriesNumber) <= 499:
                            PathFlowDataSI = dirName
                            
                            
                    else:
                        proceed = False
                        print("We currently can not load files from " + ds.Manufacturer + ".")
                        sys.exit()
          
    MagPathStr = str(FolderPath)
    PathList=MagPathStr.split("/")
    basePath = MagPathStr.replace(PathList[-1],"")


    flowData = None 
    folderPath = PathFlowDataMAG
           
    lstFilesDCM = []
          
            ################## Reading time of flight files
            # listing magnitude files
    for dirName, subdirList, fileList in os.walk( folderPath + "/"):
        for filename in fileList:
            # if ".dcm" in filename.lower():  # check whether the file's DICOM
            lstFilesDCM.append(os.path.join(dirName,filename))
                
                
            # Get ref file
    RefDs = dicom.read_file(lstFilesDCM[0])
    print(colored.blue("\t Patient ID: " + RefDs.PatientID ))
    print(colored.blue("\t Manufacturer Name: " + RefDs.Manufacturer ))
    #print("M: " + RefDs.SoftwareVersion )
            
        

def printReport(outPath, RefDs):
    # file-output.py
    today = datetime.date.today()
    dXY = RefDs.PixelSpacing
    dZ = RefDs.SpacingBetweenSlices
    pixel_spc = (dXY[0],dXY[1],dZ)
    f = open(outPath + "/readMe",'w')
    f.write('This is the report for reading GE produced DICOM files. \n In case any problems contact: ali.bakhshinejad@gmail.com \n Produced at' + str(today))
    f.write('--'*20)
    f.write('\n Patient information')
    f.write('\n Patient Name: ' + RefDs.PatientName)
    f.write('\n Patient ID: ' + RefDs.PatientID)
    f.write('\n Patient Position: ' + RefDs.PatientPosition)
    f.write('--'*5)
    f.write('\n Image information:')
  #  f.write('\n Image Orientation Position: ' + RefDs.ImageOrientationPosition)
    f.write('\n Resolution: ' + str(pixel_spc))
    f.close() 
    
def main():

    parser = argparse.ArgumentParser(description="GE 4D Flow DICOM reader developed by Ali Bakhshinejad. contact: ali.bakhshinejad@gmail.com")

    parser.add_argument("-i", "--input", help="Path to the main folder.")
    parser.add_argument("-v", "--velocityorder", help="The order of reading velocity files, default value is [1,0,2] which reresents [y,x,z]")
    parser.add_argument("-s", "--velocitysign", help="Sign for each velocity component, default value is [1,1,-1]")
    parser.add_argument("-e", "--eddycurrent", action="store_true", help="Activating Eddy current correction function")
    parser.add_argument("-p", "--eddyplane", type=int, help="The plane order to fit on the static tissue. Currently we support 1st and second order (value: 1 or 2)")
    parser.add_argument("-t", "--eddythreshold", type=int, help="The threshold value to generate static tissue mask (default value is standard deviation less than 20)")
    parser.add_argument("-n", "--randomnoise", type=int, help="Threshold for random noise correction.(Default is 60)")
    parser.add_argument("-ol", "--output", help="Output location")
    parser.add_argument("--vtk", action="store_true", help="save in VTK format")
    parser.add_argument("--mat", action="store_true", help="save in MAT format")

    args = parser.parse_args()


    if args.input is None:
        print(colored.red("FatalError: Input location is missing."))
        sys.exit()
    else:
        print(colored.green("We are looking to read data from: "))
        readPatientInfo(args.input)

    if args.velocityorder is None:
        args.velocityorder = numpy.array([1,0,2])
    else:
        args.velocityorder = numpy.array(args.velocityorder)

    if args.velocitysign is None:
        args.velocitysign = numpy.array([1,1,-1])
    else:
        args.velocitysign = numpy.array(args.velocitysign)

    if args.output is None:
        print(colored.red("FatalError: output location is missing."))
        sys.exit()
    else:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
   
    
    RefDs =  readGEFlow(args)
    printReport(args.output, RefDs)

    print(colored.magenta("Done!"))
    
    
main()
    

