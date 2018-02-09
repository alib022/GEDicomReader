import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse
from clint.textui import colored
from vtk.util import numpy_support
from saveVTKSeg import saveVTKSeg
from saveVTK import saveVTK

def readSIEMENSFlow(args):
          
    print( colored.green("\nLooking for flow data... \n"))
    MagPathStr = args.input
    foldersList = [os.path.join(MagPathStr,o) for o in os.listdir(MagPathStr) if os.path.isdir(os.path.join(MagPathStr,o))]
        
    if not foldersList:
        filesListTEMP = glob.glob(MagPathStr + "/*") 
            
        ds = dicom.read_file(filesListTEMP[0])
        if "SIEMENS" in ds.Manufacturer: 
                print("It's SIEMENS sequence!")
        else:
                print("We currently can not load files from " + ds.Manufacturer + ".")
    else:
            
                for dirName in foldersList:
                    filesListTEMP = glob.glob(dirName + "/*") 
                    
                    ds = dicom.read_file(filesListTEMP[0])
                    if "SIEMENS" in ds.Manufacturer:
                        
                        proceed = True
                        if "magnitude"  in ds.ImageComments:
                            PathFlowDataMAG = dirName
                            print(colored.cyan("Mag folder is: " + PathFlowDataMAG))
                        if "phase" in ds.ImageComments:
                            PathFlowData = dirName
                            print(colored.cyan("Flow folder is: " + PathFlowData))
                            #ConstDimsTemp = (int(ds.Rows), int(ds.Columns), len(filesListTEMP))
                            #dXY = ds.PixelSpacing
                            #dZ = ds.SpacingBetweenSlices
                            #pixel_spc = (dXY[0],dXY[1],dZ)
                            #print(pixel_spc)
                        
                            
                    else:
                        proceed = False
                        print(colored.red("FatalError: We currently can not load files from " + ds.Manufacturer + "."))
                        sys.exit()
          
 #   MagPathStr = str(FolderPath)
   # PathList=MagPathStr.split("/")
   # basePath = MagPathStr.replace(PathList[-1],"")


    if proceed:
      flowData = None 
      for folderNumber in range(0,2):
    #    sliceLocation = []
        # flow files list
        lstFilesDCM = []
     #   triggerTime = []
            
        if folderNumber == 0:
            folderPath = PathFlowDataMAG
            print(colored.cyan("Reading the Magnitude files."))
           
        if folderNumber == 1:
            if args.segmentation:
                break
            folderPath = PathFlowData
            print(colored.cyan("Reading the flow files ."))
           
       
           
            
          
            ################## Reading time of flight files
            # listing magnitude files
        for dirName, subdirList, fileList in os.walk( folderPath + "/"):
            for filename in fileList:
            # if ".dcm" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))
      #          ds = dicom.read_file(lstFilesDCM[-1])
      #          sliceLocation.append(ds.SliceLocation)
      #          triggerTime.append(ds.TriggerTime)
                
            # Get ref file
        RefDs = dicom.read_file(lstFilesDCM[0])
            
            
      #  triggerTimeTemp = sorted(set(triggerTime), key=float)
        #print(triggerTimeTemp)
        #sliceLocationTemp = set(sliceLocation)       
       # sliceLocationTemp = sorted(set(sliceLocation), key=float)
        #print(sliceLocationTemp)

       

            
        if folderNumber == 0:
                ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
                dXY = ds.PixelSpacing
                dZ = ds.SliceThickness
                pixel_spc = (dXY[0],dXY[1],dZ)
                MagData = numpy.zeros(ConstPixelDims, dtype=numpy.double)
                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)
                    MagData[:,:,int(dsTemp.InstanceNumber)-1]= dsTemp.pixel_array.astype('float')
             
                #print(MagData.shape)
               # print(RefDs)
                MagDataAve = numpy.reshape(MagData, (int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM)/int(RefDs.CardiacNumberofImages), RefDs.CardiacNumberofImages), order='F')
                #print(MagDataAve.shape)
                magDataTemp = MagDataAve.mean(3)
                
               # magDataRevX = numpy.flip(magDataTemp, 0)
               # magDataRevZ = numpy.flip(magDataRevX, 2)
                
                magSize = magDataTemp.shape
                totalNodes = magSize[0] * magSize[1] * magSize[2]
                
                if args.vtk:
                    saveVTKSeg(magDataTemp,False,False, pixel_spc, totalNodes, args.output)
                if args.mat:
                    scipy.io.savemat(args.output + "/mag.mat", mdict={'magDataTemp': magDataTemp})
             #   numpy.save(args.output +"/mag", magDataTemp) 
                
        else :
                if flowData is None:
                    ConstFlowPixelDims = (int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM))
                    flowDataTemp = numpy.zeros(ConstFlowPixelDims, dtype=numpy.double)
                    sliceLocationTemp = numpy.zeros([len(lstFilesDCM),2], dtype=numpy.double)
                   # icounter = 0
                for iFile in lstFilesDCM:
                    dsTemp = dicom.read_file (iFile)        
                    flowDataTemp[:,:,int(dsTemp.InstanceNumber)-1]= dsTemp.pixel_array.astype('int')
                    sliceLocationTemp[int(dsTemp.InstanceNumber)-1] = numpy.array([int(dsTemp.InstanceNumber)-1, float(dsTemp.SliceLocation)])
                   # icounter += 1
                    
                    
                print(sliceLocationTemp.shape)
                print(sliceLocationTemp)
                sys.exit()
                #print(flowDataTemp.shape)
                flowDataTemp = numpy.reshape(flowDataTemp, (int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM)/3, 3), order='F')
                
                UAll = flowDataTemp[:,:,:,args.velocityorder[0]].squeeze()
                VAll = flowDataTemp[:,:,:,args.velocityorder[1]].squeeze()
                WAll = flowDataTemp[:,:,:,args.velocityorder[2]].squeeze()
                
                UOrg = numpy.zeros([int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM)/30,10], dtype=numpy.double)
                VOrg = numpy.zeros([int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM)/30,10], dtype=numpy.double)
                WOrg = numpy.zeros([int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM)/30,10], dtype=numpy.double)
                
                UOrg[:,:,:,0] = UAll[:,:,0:30]
                UOrg[:,:,:,1] = UAll[:,:,30:60]
                
                VOrg[:,:,:,0] = VAll[:,:,0:30]
                VOrg[:,:,:,1] = VAll[:,:,30:60]
                
                WOrg[:,:,:,0] = WAll[:,:,0:30]
                WOrg[:,:,:,1] = WAll[:,:,30:60]
                
               # flowDataTemp = numpy.reshape(flowDataTemp, (int(RefDs.Rows), int(RefDs.Columns),len(lstFilesDCM)/3/int(RefDs.CardiacNumberofImages), 3, int(RefDs.CardiacNumberofImages)), order='F')
                #sliceLocationTemp = numpy.reshape(sliceLocationTemp, (len(lstFilesDCM)/3,2, 3), order='F')
                
                #flowDataRevX = numpy.flip(flowDataTemp, 0)
                #flowDataRevZ = numpy.flip(flowDataRevX, 2)
#                sliceLocations = numpy.unique(sliceLocationTemp[:,1])
#                print(len(sliceLocations))
#                sliceLocationsIndeces = (numpy.array(numpy.where(sliceLocationTemp[:,1] == sliceLocations[0])))
#                print(numpy.array(numpy.where(sliceLocationTemp[:,1] == sliceLocations[0])).shape)
#                sys.exit()
                
                
#                for jj in range(0,len(sliceLocations)-1):
#                    #for ii in range(0,int(RefDs.CardiacNumberofImages)-1):
#                        if sliceLocationsIndeces[1,jj] < 300:
#                            UOrg[:,:,jj,ii] = flowDataTemp[:,:,jj,0].squeeze()
#                    
                        
                 #USlice = numpy.sort(sliceLocationTemp[0:300], axis=1)
                  #      print(USlice)
                      
                        #UOrg = flowDataTemp[:,:,ii,0].squeeze()
                        #USliceLocation = sliceLocationTemp[ii,:]
                        #USliceLocationSort = numpy.sort(USliceLocation, axis=1)
             #   sliceUnique = numpy.unique(sliceLocationTemp[0:300,:], axis = 0)
             #   print(sliceUnique)
             #           
                #    elif 299 < ii < 600:
                 #       VOrg = flowDataTemp[:,:,ii,1].squeeze()
                  #      VOrgSliceLocation = sliceLocationTemp[ii,:]
                        
                  #  elif 599 < ii < 900:
                  #      WOrg = flowDataTemp[:,:,ii,2].squeeze()
                  #      WOrgSliceLocation = sliceLocationTemp[ii,:]
                        
                        
#                selector = 0
#                if args.segmentation is False:
#                     ### The combination of x +y and -z and permuted x and y is working. Ali Aug24 2017
#                     UOrg = args.velocitysign[0] * (flowDataTemp[:, :, :, args.velocityorder[0]].squeeze())
#                     VOrg = args.velocitysign[1] * (flowDataTemp[:, :, :, args.velocityorder[1]].squeeze())
#                     WOrg = args.velocitysign[2] * (flowDataTemp[:, :, :, args.velocityorder[2]].squeeze())
#        
                flowCorrected = numpy.zeros([flowDataTemp.shape[0], flowDataTemp.shape[1], flowDataTemp.shape[2]/10,3,10])
#                #if args.eddycurrent:
#                #    flowCorrected = eddyCurrentCorrection(UOrg, VOrg, WOrg, args.randomnoise, args.eddythreshold, args.eddyplane)

                #elif args.randomnoise is not None:
                #    flowCorrected = randNoise(args, UOrg, VOrg, WOrg, args.randomnoise/100, 0, 1)
            
            
                #else:
                flowCorrected[:, :, :,0, :] = UOrg
                flowCorrected[:, :, :,1, :] = VOrg
                flowCorrected[:, :, :,2, :] = WOrg
        
                
                if args.mat:
                    scipy.io.savemat(args.output + "/vel.mat", mdict={'flowData': flowData})
                #print(flowData.shape)
                
                if args.vtk:
    #    if args.segmentation:
            
     #       saveVTKSeg(magDataTemp,False,False, pixel_spc, totalNodes, args.output)
     #   else:
            
                    saveVTK(magDataTemp, flowCorrected,pixel_spc, totalNodes, args.output)
    
                
                


   
          
    #    if args.eddycurrent:
    #        flowCorrected = eddyCurrentCorrection(UOrg, VOrg, WOrg, args.randomnoise, args.eddythreshold, args.eddyplane)

    #    elif args.randomnoise is not None:
    #        flowCorrected = randNoise(args, UOrg, VOrg, WOrg, args.randomnoise/100, 0, 1)
            
            
    #    else:
     #       flowCorrected[:, :, :,0, :] = UOrg
     #       flowCorrected[:, :, :,1, :] = VOrg
     #       flowCorrected[:, :, :,2, :] = WOrg
            
        
        
        
    print(colored.green("\nGetting ready to write files... This takes a little bit of time"))
    
   # magSize = magDataTemp.shape
   #totalNodes = magSize[0] * magSize[1] * magSize[2]

    #if (args.vtk == False and args.mat == False):
     #   print(colored.yellow("We will ONLY save in npy format, since you didnt select your preference! (VTK or MAT)"))
        
    #if not args.segmentation:
     #   numpy.save(args.output +"/FlowData", flowCorrected)    
    #if args.vtk:
    #    if args.segmentation:
            
     #       saveVTKSeg(magDataTemp,False,False, pixel_spc, totalNodes, args.output)
     #   else:
            
      #      saveVTK(magDataTemp, flowCorrected,pixel_spc, totalNodes, args.output)
    
    
    #if args.mat:
    #    if args.segmentation:
    #        with open(args.output + "/FlowData.mat", 'wb') as matlabFile:
    #            scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
            
    #    else:
     #       with open(args.output + "/FlowData.mat", 'wb') as matlabFile:
      #          scipy.io.savemat(matlabFile, mdict={'velocity': flowData})
       #         scipy.io.savemat(matlabFile, mdict={'magnitude': magDataTemp})
    
    
    
    return RefDs
 #   printReport(outPath, RefDs)
