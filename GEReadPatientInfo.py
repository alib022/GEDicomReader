import dicom,os, glob, sys, GEReadInfo
from clint.textui import colored


def readPatientInfo(FolderPath, cmra, tof):
    Version = "0.0.2"
    MagPathStr = str(FolderPath)
    foldersList = [os.path.join(MagPathStr,o) for o in os.listdir(MagPathStr) if os.path.isdir(os.path.join(MagPathStr,o))]
 
    
    if not foldersList:
        filesListTEMP = glob.glob(MagPathStr + "/*") 
            
        ds = dicom.read_file(filesListTEMP[0])
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
