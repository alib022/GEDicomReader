import numpy, sys, eddyNoise




flowData = numpy.load("FlowData.npy")
print(flowData.shape)

UOrg = flowData[:,:,:,0,:]
VOrg = flowData[:,:,:,1,:]
WOrg = flowData[:,:,:,2,:]


flowCorrectedNoise = eddyNoise.randNoise(UOrg, VOrg, WOrg, plotBool=0)


flowCorrectedEddy = eddyNoise.eddyCurrentCorrection(flowCorrectedNoise[:,:,:,0,:], flowCorrectedNoise[:,:,:,1,:], flowCorrectedNoise[:,:,:,2,:])
