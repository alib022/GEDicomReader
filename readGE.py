import dicom,os, glob, scipy.io, numpy, vtk, sys, datetime, argparse, timeit, math
from clint.textui import colored
from DICOMClasses import PatientData
from readGEFlow import readGEFlow, readGEcMRA, readGETOF
from GEReadPatientInfo import readPatientInfo

''' This function reads GE Flow data '''

def printReport(outPath, PatientDataStruc, Version):
    # file-output.py
    today = datetime.date.today()

 
    f = open(outPath + "/readMe.txt",'w')

    f.write("\n\t"+"**"*20)
    f.write("\n\t\tGE 4D Flow DICOM reader. \n\t\tDeveloped by: Ali Bakhshinejad \n\t\tali.bakhshinejad@gmail.com \n\t\t Version="+ Version)
    f.write("\n\t"+"**"*20)
    
    f.write("\nReading data for case:")
    f.write("\n\t Patient ID: " + PatientDataStruc.PatientID )
    f.write("\n\t Manufacturer Name: " + PatientDataStruc.Manufacturer )
 

 #   f.write('This is the report for reading GE produced DICOM files. \n In case any problems contact: ali.bakhshinejad@gmail.com \n Produced at ' + str(today))
#    f.write('\n' + '--'*20)
#    f.write('\n Patient information')
 #   f.write('\n Patient Name: ' + RefDs.PatientName)
#    f.write('\n Patient ID: ' + RefDs.PatientID)
#    f.write('\n Patient Position: ' + RefDs.PatientPosition)
#    f.write('\n'+'--'*5)
    f.write('\nImage information:')
 #   f.write('\n Image Orientation Position: ' + RefDs.ImageOrientationPosition)
    f.write("\n\t Scan resolution: " + str(PatientDataStruc.PixelSize) )
    f.close() 
    
def main():

    start_time = timeit.default_timer()

    parser = argparse.ArgumentParser(description="GE 4D Flow DICOM reader developed by Ali Bakhshinejad. contact: ali.bakhshinejad@gmail.com")

    parser.add_argument("-i", "--input", help="Path to the main folder.")
    parser.add_argument("-v", "--velocityorder", help="The order of reading velocity files, default value is [1,0,2] which reresents [y,x,z]")
    parser.add_argument("-si", "--velocitysign", help="Sign for each velocity component, default value is [1,1,-1]")
    parser.add_argument("-e", "--eddycurrent", action="store_true", help="Activating Eddy current correction function")
    parser.add_argument("-p", "--eddyplane", type=int, help="The plane order to fit on the static tissue. Currently we support 1st and second order (value: 1 or 2). Default value is 2nd order polynominal.")
    parser.add_argument("-t", "--eddythreshold", type=int, help="The threshold value to generate static tissue mask (default value is standard deviation less than 20)")
    parser.add_argument("-n", "--randomnoise", help="Threshold for random noise correction. (In percentage)")
    parser.add_argument("-ol", "--output", help="Output location")
    parser.add_argument("--vtk", action="store_true", help="save in VTK format")
    parser.add_argument("--mat", action="store_true", help="save in MAT format")
    parser.add_argument("-se", "--segmentation",  action="store_true", help="Only save magnitude file to be used for segmentation purposes.")

    parser.add_argument("--cmra", action="store_true", help="Read cMRA dataset, (No Flow Data).")
    parser.add_argument("--tof", action="store_true", help="Read Time Of Flght (TOF) database.")

    inputFlags = parser.parse_args()
    

    if inputFlags.input is None:
        print(colored.red("FatalError: Input location is missing."))
        sys.exit()
    else:
        #print(colored.green("We are looking to read data from: "))
        PatientDataStruc, Version = readPatientInfo(inputFlags.input, inputFlags.cmra, inputFlags.tof)

    if inputFlags.velocityorder is None:
        inputFlags.velocityorder = numpy.array([1,0,2])
    else:
        inputFlags.velocityorder = numpy.array(inputFlags.velocityorder)

    if inputFlags.velocitysign is None:
        inputFlags.velocitysign = numpy.array([1,1,-1])
    else:
        inputFlags.velocitysign = numpy.array(inputFlags.velocitysign)

    if inputFlags.output is None:
        print(colored.red("FatalError: output location is missing."))
        sys.exit()
    else:
        if not os.path.exists(inputFlags.output):
            os.makedirs(inputFlags.output)

    if inputFlags.eddythreshold is None:
        inputFlags.eddythreshold = 20

    if inputFlags.randomnoise is None:
        print(colored.yellow("Warning: No random noise correction will happen!"))


    if (inputFlags.eddycurrent and inputFlags.eddyplane is None):
        inputFlags.eddyplane = 2

    
    
   
    #print(args)
    if inputFlags.cmra:
        readGEcMRA(inputFlags, PatientDataStruc)

    elif inputFlags.tof:
        readGETOF(inputFlags, PatientDataStruc)
    else:
        readGEFlow(inputFlags, PatientDataStruc)

    printReport(inputFlags.output, PatientDataStruc, Version)
    # code you want to evaluate
    elapsed = timeit.default_timer() - start_time
    print(colored.yellow("Execuation time: " + str(elapsed)+ " s \nDone!"))
    
    
main()
    

