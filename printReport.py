import sys, datetime


def printReport(outPath, PatientDataStruc, Version, fileName):
    # file-output.py
    today = datetime.date.today()

 
    f = open(outPath + "/" + fileName + "ReadMe.txt",'w')

    f.write("\n\t"+"**"*20)
    f.write("\n\t\tGE 4D Flow DICOM reader. \n\t\tDeveloped by: Ali Bakhshinejad \n\t\tali.bakhshinejad@gmail.com \n\t\t Version="+ Version)
    f.write("\n\t"+"**"*20)
    
    f.write("\nReading data for case:")
    f.write("\n\t Patient ID: " + PatientDataStruc.PatientID )
    f.write("\n\t Manufacturer Name: " + PatientDataStruc.Manufacturer )
 

    f.write('\nImage information:')
    f.write("\n\t Scan resolution: " + str(PatientDataStruc.PixelSize) )
    f.close() 
