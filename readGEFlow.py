import dicom,os, glob, scipy.io, numpy, vtk, sys, saveVTK, math, eddyNoise
from clint.textui import colored
from vtk.util import numpy_support



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
 #                   print(PatientDataStruc.FlowVecSize)
                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)        
                    flowData[:,:,sliceLocationTemp.index(dsTemp.SliceLocation), folderNumber-1,triggerTimeTemp.index(dsTemp.TriggerTime)]= dsTemp.pixel_array.astype('float')
    
                if inputFlags.mat:
                    scipy.io.savemat(inputFlags.output + "/vel.mat", mdict={'flowData': flowData})
                #print(flowData.shape)


    if inputFlags.segmentation is False:
        ### The combination of -x +y and -z and permuted x and y is working. Ali Aug24 2017
        UOrg = inputFlags.velocitysign[0] * (flowData[:, :, :, inputFlags.velocityorder[0]].squeeze())
        VOrg = inputFlags.velocitysign[1] * (flowData[:, :, :, inputFlags.velocityorder[1]].squeeze())
        WOrg = inputFlags.velocitysign[2] * (flowData[:, :, :, inputFlags.velocityorder[2]].squeeze())
        
        flowCorrected = numpy.zeros([flowData.shape[0], flowData.shape[1], flowData.shape[2],3,flowData.shape[4]])

        
          
        if inputFlags.eddycurrent:
            flowCorrected = eddyNoise.eddyCurrentCorrection(UOrg, VOrg, WOrg, inputFlags.randomnoise, inputFlags.eddythreshold, inputFlags.eddyplane)

        elif inputFlags.randomnoise is not None:
            flowCorrected = eddyNoise.randNoise(inputFlags, UOrg, VOrg, WOrg, int(inputFlags.randomnoise)/100, 0, 0)
            
            
        else:
            flowCorrected[:, :, :,0, :] = UOrg
            flowCorrected[:, :, :,1, :] = VOrg
            flowCorrected[:, :, :,2, :] = WOrg
            
        
        vTemp = numpy.amax(flowCorrected, axis=4)
        #flowCorrected.mean(4)
        vTemp2 = numpy.sqrt(vTemp[:,:,:,0]**2 + vTemp[:,:,:,1]**2+ vTemp[:,:,:,2]**2)
        vCMRA = numpy.multiply(vTemp2, magDataTemp)
        
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
    magSize = magDataTemp.shape
    totalNodes = magSize[0] * magSize[1] * magSize[2]

    if (inputFlags.vtk == False and inputFlags.mat == False):
        print(colored.yellow("We will ONLY save in npy format, since you didnt select your preference! (VTK or MAT)"))
        numpy.save(inputFlags.output +"/FlowData", flowCorrected)
        
   # if not inputFlags.segmentation:
   #     numpy.save(inputFlags.output +"/FlowData", flowCorrected)    
    if inputFlags.vtk:
        if inputFlags.segmentation:
            
            saveVTK.saveVTKSeg(magDataTemp,False,False, PatientDataStruc.PixelSize, totalNodes, inputFlags.output)
            
        else:
            
            saveVTK.saveVTK(vCMRA, flowCorrected,  PatientDataStruc.PixelSize, totalNodes, inputFlags.output)
            saveVTK.saveVTKSeg(vCMRA,True,False, PatientDataStruc.PixelSize, totalNodes, inputFlags.output)
            saveVTK.saveVTKSeg(magDataTemp,False,False, PatientDataStruc.PixelSize, totalNodes, inputFlags.output)
    
    
    if inputFlags.mat:
        if inputFlags.segmentation:
            with open(inputFlags.output + "/FlowData.mat", 'wb') as matlabFile:
                scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
            
        else:
            with open(inputFlags.output + "/FlowData.mat", 'wb') as matlabFile:
                scipy.io.savemat(matlabFile, mdict={'velocity': flowData})
                scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
    
    


        
