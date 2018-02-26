import vtk, numpy
from vtk.util import numpy_support

def saveVTK(magDataTemp, flowCorrected,pixel_spc, totalNodes, outPath):
        
        for timeIter in range(0,flowCorrected.shape[4]):
        
            magDataOut = numpy.reshape(magDataTemp, (totalNodes, 1), order='F')
        
            flowDataOut = numpy.reshape(flowCorrected[:,:,:,:,timeIter], (totalNodes, 3), order='F')
        
        #        #########################
        #        # Convert the VTK array to vtkImageData
        
            img_vtk = vtk.vtkImageData()
            img_vtk.SetDimensions(magDataTemp.shape)
            img_vtk.SetSpacing(pixel_spc)
            img_vtk.SetOrigin(0,0,0)
    
        
            magArray = vtk.vtkDoubleArray()
            magArray.SetNumberOfComponents(1)
            magArray.SetNumberOfTuples(img_vtk.GetNumberOfPoints())
            magArray.SetName("Magnitude")
    
            for i in range(0,totalNodes):
                magArray.SetValue(i,magDataOut[i])
        
            flowArray = vtk.vtkDoubleArray()
            flowArray.SetNumberOfComponents(3)
            flowArray.SetNumberOfTuples(img_vtk.GetNumberOfPoints())
            flowArray.SetName("Velocity")
            for i in range(0,totalNodes):
                flowArray.SetTuple3(i,flowDataOut[i,0], flowDataOut[i,1], flowDataOut[i,2])
        
            img_vtk.GetPointData().AddArray(magArray)
            img_vtk.GetPointData().AddArray(flowArray)
        
        
            writer = vtk.vtkXMLImageDataWriter()
            writer.SetFileName(outPath + '/FlowData_%d.vti' % timeIter)
            writer.SetInputData(img_vtk)
            ## This is set so we can see the data in a text editor.
            writer.SetDataModeToAscii()
            writer.Write()


def saveVTKSeg(magDataTemp, cMRA,TOF, pixel_spc, totalNodes, outPath):
        
        
            MagVTK = numpy_support.numpy_to_vtk(num_array=magDataTemp.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
                    
            img_vtk = vtk.vtkImageData()
            img_vtk.SetDimensions(magDataTemp.shape)
            img_vtk.AllocateScalars(vtk.VTK_FLOAT,1)
            img_vtk.SetSpacing(pixel_spc)
            img_vtk.SetOrigin(0,0,0)
            img_vtk.GetPointData().SetScalars(MagVTK)        
        
                
            writer = vtk.vtkXMLImageDataWriter()

            if cMRA:
                img_vtk.SetName("cMRA")
                writer.SetFileName(outPath + '/cMRAData.vti')
            elif TOF:
                img_vtk.SetName("TOF")
                writer.SetFileName(outPath + '/TOFData.vti')                
            else:
                img_vtk.SetName("Magnitude")
                writer.SetFileName(outPath + '/MagData.vti')

            writer.SetInputData(img_vtk)
            ## This is set so we can see the data in a text editor.
            writer.SetDataModeToAscii()
            writer.Write()

            return 0

