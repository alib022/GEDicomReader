import numpy
from eddyNoise import eddyCurrentCorrection



flowData = numpy.load("FlowData.npy")
print(flowData.shape)

UOrg = flowData[:,:,:,0,:]
VOrg = flowData[:,:,:,1,:]
WOrg = flowData[:,:,:,2,:]


eddyCurrentCorrection(UOrg, VOrg, WOrg)
