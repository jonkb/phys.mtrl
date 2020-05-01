from region import Region

class Member:
	def __init__(self, material, xsection, length):
		self.material = material
		self.xsection = xsection
		self.length = length
	def __str__(self):
		s = self.material["name"] + " member with cross section=("
		s += str(self.xsection) + "), and length=" + str(self.length)
		return s
	@property
	def E(self):
		return self.material["E"]
	@property
	def sig_y(self):
		return self.material["sig_y"]
	@property
	def sig_u(self):
		return self.material["sig_u"]
	@property
	def xarea(self):
		return self.xsection.area
	@property
	def xIx(self):
		return self.xsection.Ix
	@property
	def xIy(self):
		return self.xsection.Iy
	@property
	def xIp(self):
		return self.xsection.Ip
	@property
	def mass(self):
		return self.xarea * self.length * self.material["rho"]

#This class is really just for reference. I'm not sure if this is the best way to do this
class Materials:
	materials = ("steel", "aluminum")
	#Used values for ASTM-A514 steel
	steel = {
		"name": "structural steel",
		"color": "#43464B",
		#Mass density
		"rho": 7850, #kg/m^3
		#Young's modulus
		"E": 200e9, #200GPa
		#yield stress
		"sig_y": 700e6, #700MPa
		#ultimate stress
		"sig_u": 830e6 #830MPa
	}
	#Used values for 6061-T6 alloy
	aluminum = {
		"name": "aluminum",
		"color": "#848789",
		#Mass density
		"rho": 2700, #kg/m^3
		#Young's modulus
		"E": 70e9, #70GPa
		#yield stress
		"sig_y": 270e6, #270MPa
		#ultimate stress
		"sig_u": 310e6 #310MPa
	}
