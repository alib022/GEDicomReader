import numpy, eddyNoise, saveVTK




flowData = numpy.load("FlowData.npy")
magData = numpy.load("mag.npy")
print("Flow shape")
print(flowData.shape)
print("Mag shape")
print(magData.shape)

UOrg = flowData[:,:,:,0,:]
VOrg = flowData[:,:,:,1,:]
WOrg = flowData[:,:,:,2,:]


flowCorrectedNoise = eddyNoise.randNoise( UOrg, VOrg, WOrg, 25)


flowCorrectedEddy = eddyNoise.eddyCurrentCorrection(flowCorrectedNoise[:,:,:,0,:], flowCorrectedNoise[:,:,:,1,:], flowCorrectedNoise[:,:,:,2,:], magData,15, 5, STDPower=2)

PixelSize = numpy.array([0.70, 0.70, 0.4])
magSize = magData.shape
totalNodes = magSize[0] * magSize[1] * magSize[2]

saveVTK.saveVTK(magData, flowCorrectedEddy,  PixelSize, totalNodes, "EddyTest/")
