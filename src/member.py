#from region import Region

#prismatic uniform members
class Member:
	def __init__(self, material, xsection, length):
		self.material = material
		self.xsection = xsection
		self.length = length
		self.placed = False
		self.img_ref = None
		self.sup = [None, None]
	def __str__(self):
		s = self.material["name"] + " member with cross section=("
		s += str(self.xsection) + "), and length=" + str(self.length)
		return s
	#Define the position in the xy plane. Units are meters.
	#vh True = vertical. vh False = horizontal.
	def place(self, x, y, vh):
		self.x0 = x
		self.y0 = y
		if vh:
			self.x1 = x
			self.y1 = y+self.length
		else:
			self.x1 = x+self.length
			self.y1 = y
		self.placed = True
	def get_pos(self):
		if self.placed:
			return (self.x0, self.y0, self.x1, self.y1)
		return None
	def is_vert(self):
		if self.placed:
			if self.x1 == self.x0:
				return True
			else:
				return False
		return None
	#side: 0 or 1 (start or end)
	#direction: 0,1,2,3 --> Up, Left, Down, Right (Pointing away from member)
	def sup_dir(self, side):
		# (vert/horiz, side) --> dir
		dir_switch = {
			(True, 0) : 2,
			(True, 1) : 0,
			(False, 0) : 1,
			(False, 1) : 3
		}
		return dir_switch[(self.is_vert(), side)]
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
