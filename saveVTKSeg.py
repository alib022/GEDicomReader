#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 16:57:45 2018

@author: ali
"""
import vtk
from vtk.util import numpy_support

def saveVTKSeg(magDataTemp, cMRA,TOF, pixel_spc, totalNodes, outPath):
        

            MagVTK = numpy_support.numpy_to_vtk(num_array=magDataTemp.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
                    
            #        #########################
            #        # Convert the VTK array to vtkImageData
            img_vtk = vtk.vtkImageData()
            img_vtk.SetDimensions(magDataTemp.shape)
            img_vtk.AllocateScalars(vtk.VTK_FLOAT,1)
            img_vtk.SetSpacing(pixel_spc)
            img_vtk.SetOrigin(0,0,0)
            img_vtk.GetPointData().SetScalars(MagVTK)        
        
   
            writer = vtk.vtkXMLImageDataWriter()

            if cMRA:
                writer.SetFileName(outPath + '/cMRAData.vti')
            elif TOF:
                writer.SetFileName(outPath + '/TOFData.vti')                
            else:
                writer.SetFileName(outPath + '/MagData.vti')

            writer.SetInputData(img_vtk)
            ## This is set so we can see the data in a text editor.
            writer.SetDataModeToAscii()
            writer.Write()

            return 0