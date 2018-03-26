import pydicom,os, glob, sys, GEReadInfo, DICOMClasses
from clint.textui import colored


def readPatientInfo(FolderPath, cmra, tof):
    Version = "0.1"
    MagPathStr = str(FolderPath)
    foldersList = [os.path.join(MagPathStr,o) for o in os.listdir(MagPathStr) if os.path.isdir(os.path.join(MagPathStr,o))]
 
    
    if not foldersList:
        
        PatientDataStruc = DICOMClasses.PatientData()        
        filesListTEMP = glob.glob(MagPathStr + "/*") 
            
        ds = pydicom.read_file(filesListTEMP[0])

        dXY = ds.PixelSpacing
        dZ = ds.SpacingBetweenSlices
        PatientDataStruc.PixelSize = (dXY[0],dXY[1],dZ)

        PatientDataStruc.PatientID = ds.PatientID
        PatientDataStruc.Manufacturer = ds.Manufacturer

        if not "GE" in ds.Manufacturer: 
                print("This is GE DICOM reader, it doesnt open " + ds.Manufacturer + ".")
                sys.exit()
    else:

        PatientDataStruc = GEReadInfo.main(foldersList)
                            
    print(colored.yellow("\t"+"**"*20))
    print(colored.yellow("\t\tGE 4D Flow DICOM reader. \n\t\tDeveloped by: Ali Bakhshinejad \n\t\tali.bakhshinejad@gmail.com \n\t\t Version="+ Version))
    print(colored.yellow("\t"+"**"*20))
    
    print(colored.green("Reading data for case:"))
    print(colored.blue("\t Patient ID: " + PatientDataStruc.PatientID ))
    print(colored.blue("\t Manufacturer Name: " + PatientDataStruc.Manufacturer ))
    print(colored.blue("\t Scan resolution: " + str(PatientDataStruc.PixelSize) ))

    return PatientDataStruc, Version
