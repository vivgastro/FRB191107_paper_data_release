class Survey(object):

  def __init__(self, name, idealize = False):
    self.name = name
    self.cname = name
    self.set_params()

    if idealize:
      self.idealize()

  def set_params(self):
    if self.name == "SUPERB":
      self.TSAMP = 64e-6    #s
      self.CHW = 0.390625e6 #Hz
      self.NCH = 1024
      self.BW = self.NCH * self.CHW
      self.CFREQ = 1381.804688e6  #Hz
      self.SEFD = 36 #Jy  TODO  -- These numbers are for HTRU taken from Chawla et al 2017
      self.NPOL = 2

    elif self.name == "UTMOST":
      self.TSAMP = 327.68e-6  #s
      self.NCH = 320
      self.BW = 31.25e6 #Hz
      self.CHW = self.BW / self.NCH
      self.CFREQ = 835.5e6 #Hz
      self.SEFD = 180 #/ 25#Jy
      self.NPOL = 1

    elif self.name == "ASKAP_FE":
      self.TSAMP = 1.26e-3   #s
      self.CHW = 1e6   #Hz
      self.NCH = 336
      self.BW = self.CHW * self.NCH
      self.CFREQ = 1320e6   #Hz
      self.NPOL = 2
      self.SEFD = 1800. #/55  #Jy   Taken from Chawla et al 2017
    
    elif self.name == "ASKAP_ICS":
      self.TSAMP = 1.728e-3   #s
      self.CHW = 1e6   #Hz
      self.NCH = 336
      self.BW = self.CHW * self.NCH
      self.CFREQ = 1320e6   #Hz
      self.NPOL = 2
      self.SEFD = 1800 / 6. #/60  #Jy   Taken from Chawla et al 2017
    
    elif self.name == "CHIME":
      self.TSAMP = 0.983e-3 #s
      self.BW = 400e6
      self.CFREQ = 600e6
      self.CHW = 24.4e3   #Hz
      self.NCH = 16000
      self.NPOL = 2
      self.SEFD = 36.2  #Jy     Taken from Chawla et al 2017

    else:
      raise ValueError("Unknown Survey name provided : {}".format(self.name))

  def idealize(self):
    self.TSAMP = 1e-68
    self.CHW = 1e-68
    self.cname = "{}_0CHW_0TS".format(self.name)

  def narrow_chans(self):
    self.CHW = 1e-16
    self.cname = "{}_0CHW".format(self.name)
  
  def reset(self):
    self.set_params()
    self.cname = self.name
