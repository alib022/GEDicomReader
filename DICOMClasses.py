

class PatientData:
    def __init__(self, MagPath=None, FlowPathRL=None, FlowPathAP=None, FlowPathSI=None, \
                 FlowVecSize=None, MagVecSize=None, PixelSize=None, PatientID=None, Manufacturer=None):
        self.MagPath = MagPath
        self.FlowPathRL = FlowPathRL
        self.FlowPathAP = FlowPathAP
        self.FlowPathSI = FlowPathSI
        self.FlowVecSize = FlowVecSize
        self.MagVecSize = MagVecSize
        self.PixelSize = PixelSize
        self.PatientID = PatientID
        self.Manufacturer = Manufacturer


