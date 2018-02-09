import math


def GE27ReadInfo(ds, dirName):

    if  int(ds.SeriesNumber) % 100 == 0:
        PathFlowDataMAG = dirName
        seriesBase = int(ds.SeriesNumber)
                            
    #if 200 <= int(ds.SeriesNumber) <= 299:
    if int(ds.SeriesNumber) == seriesBase+1:
        PathFlowDataRL = dirName
                            
 #      ConstDimsTemp = (int(ds.Rows), int(ds.Columns), int(ds.ImagesInAcquisition), int(ds.CardiacNumberOfImages))
        ConstDimsTemp = (int(ds.Rows), int(ds.Columns), math.ceil(len(filesListTEMP)/ int(ds.CardiacNumberOfImages)), int(ds.CardiacNumberOfImages))

        dXY = ds.PixelSpacing
        dZ = ds.SpacingBetweenSlices
        pixel_spc = (dXY[0],dXY[1],dZ)
        #print(pixel_spc)
        if int(ds.SeriesNumber) == seriesBase+2:
        #if 300 <= int(ds.SeriesNumber) <= 399:
            PathFlowDataAP = dirName
                            
                       # if 400 <= int(ds.SeriesNumber) <= 499:
        if int(ds.SeriesNumber) == seriesBase+3:
            PathFlowDataSI = dirName

        return PathFlowDataMAG, PathFlowDataRL, PathFlowDataAP, PathFlowDataSI 
