import math, pydicom, glob, DICOMClasses


def main(foldersList):
   

    PatientDataStruc = DICOMClasses.PatientData()

    for dirName in foldersList:
        filesListTEMP = glob.glob(dirName + "/*") 
                    
        ds = pydicom.read_file(filesListTEMP[0])
        GESoftwareVersion = int(ds.SoftwareVersions[0])

        if GESoftwareVersion == 27:

            if  int(ds.SeriesNumber) % 100 == 0:
                PatientDataStruc.MagPath = dirName
                seriesBase = int(ds.SeriesNumber)

    for dirName in foldersList:

        filesListTEMP = glob.glob(dirName + "/*") 
                    
        ds = pydicom.read_file(filesListTEMP[0])
        GESoftwareVersion = int(ds.SoftwareVersions[0])

        if GESoftwareVersion == 27:


            if  int(ds.SeriesNumber) % 100 == 0:
                PatientDataStruc.MagPath = dirName
                seriesBase = int(ds.SeriesNumber)

            if int(ds.SeriesNumber) == seriesBase+1:
                PatientDataStruc.FlowPathRL = dirName
                PatientDataStruc.MagVecSize = (int(ds.Rows), int(ds.Columns), int(math.ceil(len(filesListTEMP)/ int(ds.CardiacNumberOfImages))), int(ds.CardiacNumberOfImages))
                PatientDataStruc.FlowVecSize = (int(ds.Rows), int(ds.Columns), int(math.ceil(len(filesListTEMP)/ int(ds.CardiacNumberOfImages))), 3, int(ds.CardiacNumberOfImages))
                dXY = ds.PixelSpacing
                dZ = ds.SpacingBetweenSlices
                PatientDataStruc.PixelSize = (dXY[0],dXY[1],dZ)


            if int(ds.SeriesNumber) == seriesBase+2:
                PatientDataStruc.FlowPathAP = dirName

            if int(ds.SeriesNumber) == seriesBase+3:
                PatientDataStruc.FlowPathSI = dirName

        elif GESoftwareVersion == 25:
            if 100 <= int(ds.SeriesNumber) <= 199:
                    PatientDataStruc.MagPath = dirName
                
                                  
            if 200 <= int(ds.SeriesNumber) <= 299:

                PatientDataStruc.FlowPathRL = dirName                               
                PatientDataStruc.MagVecSize = (int(ds.Rows), int(ds.Columns), int(ds.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
                PatientDataStruc.FlowVecSize = (int(ds.Rows), int(ds.Columns), int(math.ceil(len(filesListTEMP)/ int(ds.CardiacNumberOfImages))), 3, int(ds.CardiacNumberOfImages))
                dXY = ds.PixelSpacing
                dZ = ds.SpacingBetweenSlices
                PatientDataStruc.PixelSize = (dXY[0],dXY[1],dZ)

            if 300 <= int(ds.SeriesNumber) <= 399:
                PatientDataStruc.FlowPathAP = dirName
                               
            if 400 <= int(ds.SeriesNumber) <= 499:
                PatientDataStruc.FlowPathSI = dirName

        PatientDataStruc.PatientID = ds.PatientID
        PatientDataStruc.Manufacturer = ds.Manufacturer

                

    return PatientDataStruc
