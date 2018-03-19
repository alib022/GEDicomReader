import numpy, eddyNoise




flowData = numpy.load("FlowData.npy")
magData = numpy.load("mag.npy")
print("Flow shape")
print(flowData.shape)
print("Mag shape")
print(magData.shape)

UOrg = flowData[:,:,:,0,:]
VOrg = flowData[:,:,:,1,:]
WOrg = flowData[:,:,:,2,:]


flowCorrectedNoise = eddyNoise.randNoise( UOrg, VOrg, WOrg, 30, 0)


#flowCorrectedEddy = eddyNoise.eddyCurrentCorrection(flowCorrectedNoise[:,:,:,0,:], flowCorrectedNoise[:,:,:,1,:], flowCorrectedNoise[:,:,:,2,:])
