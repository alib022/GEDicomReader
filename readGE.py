import dicom,os, glob, scipy.io, numpy, vtk, sys, argparse, timeit, math, printReport, readGEFlow, readGEMRA
from clint.textui import colored
from DICOMClasses import PatientData
from GEReadPatientInfo import readPatientInfo

''' This function reads GE Flow data '''
    
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
        PatientDataStruc, Version = readPatientInfo(inputFlags.input, inputFlags.cmra, inputFlags.tof)

    if inputFlags.velocityorder is None:
        inputFlags.velocityorder = numpy.array([1,0,2])
    else:
        inputFlags.velocityorder = numpy.array(inputFlags.velocityorder)

    if inputFlags.velocitysign is None:
        inputFlags.velocitysign = numpy.array([-1,1,-1])
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

    if inputFlags.cmra:
        readGEMRA.readGEcMRA(inputFlags, PatientDataStruc)

    elif inputFlags.tof:
        readGEMRA.readGETOF(inputFlags, PatientDataStruc)
    else:
        readGEFlow.readGEFlow(inputFlags, PatientDataStruc)

    if inputFlags.segmentation:
        printReport.printReport(inputFlags.output, PatientDataStruc, Version, "seg")
    elif inputFlags.tof:
        printReport.printReport(inputFlags.output, PatientDataStruc, Version, "tof")
    elif inputFlags.cmra:
        printReport.printReport(inputFlags.output, PatientDataStruc, Version, "cMRA")
    else:
        printReport.printReport(inputFlags.output, PatientDataStruc, Version, "flow")


    # code you want to evaluate
    elapsed = timeit.default_timer() - start_time
    print(colored.yellow("Execuation time: " + str(elapsed)+ " s \nDone!"))
    
    
main()
    

