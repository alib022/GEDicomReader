import numpy, sys, eddyNoise




flowData = numpy.load("FlowData.npy")
print(flowData.shape)

UOrg = flowData[:,:,:,0,:]
VOrg = flowData[:,:,:,1,:]
WOrg = flowData[:,:,:,2,:]


flowCorrected = eddyNoise.randNoise(UOrg, VOrg, WOrg, plotBool=0)


eddyNoise.eddyCurrentCorrection(flowCorrected[:,:,:,0,:], flowCorrected[:,:,:,1,:], flowCorrected[:,:,:,2,:])
