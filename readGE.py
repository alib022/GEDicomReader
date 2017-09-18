import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse, timeit
from clint.textui import colored
from readGEFlow import readGEFlow

''' This function reads GE Flow data '''


        


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
    print(colored.yellow("\t"+"**"*20))
    print(colored.yellow("\t\tGE 4D Flow DICOM reader. \n\t\tDeveloped by: Ali Bakhshinejad \n\t\tali.bakhshinejad@gmail.com"))
    print(colored.yellow("\t"+"**"*20))
    
    
    
    print(colored.green("Reading data for case:"))
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

    start_time = timeit.default_timer()

    parser = argparse.ArgumentParser(description="GE 4D Flow DICOM reader developed by Ali Bakhshinejad. contact: ali.bakhshinejad@gmail.com")

    parser.add_argument("-i", "--input", help="Path to the main folder.")
    parser.add_argument("-v", "--velocityorder", help="The order of reading velocity files, default value is [1,0,2] which reresents [y,x,z]")
    parser.add_argument("-si", "--velocitysign", help="Sign for each velocity component, default value is [1,1,-1]")
    parser.add_argument("-e", "--eddycurrent", action="store_true", help="Activating Eddy current correction function")
    parser.add_argument("-p", "--eddyplane", type=int, help="The plane order to fit on the static tissue. Currently we support 1st and second order (value: 1 or 2)")
    parser.add_argument("-t", "--eddythreshold", type=int, help="The threshold value to generate static tissue mask (default value is standard deviation less than 20)")
    parser.add_argument("-n", "--randomnoise", type=int, help="Threshold for random noise correction.(Default is 60)")
    parser.add_argument("-ol", "--output", help="Output location")
    parser.add_argument("--vtk", action="store_true", help="save in VTK format")
    parser.add_argument("--mat", action="store_true", help="save in MAT format")
    parser.add_argument("-se", "--segmentation",  action="store_true", help="Only save magnitude file to be used for segmentation purposes.")

    args = parser.parse_args()
    

    if args.input is None:
        print(colored.red("FatalError: Input location is missing."))
        sys.exit()
    else:
        #print(colored.green("We are looking to read data from: "))
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

    if args.eddythreshold is None:
        args.eddythreshold = 20

    if args.randomnoise is None:
        args.randomnoise = 60

    
    
   
    #print(args)
    RefDs =  readGEFlow(args)
    printReport(args.output, RefDs)
    # code you want to evaluate
    elapsed = timeit.default_timer() - start_time
    print(colored.yellow("Execuation time: " + str(elapsed)+ " s \nDone!"))
    
    
main()
    

