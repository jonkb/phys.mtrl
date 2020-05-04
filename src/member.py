import numpy
#from region import Region
eps = 1e-15

def N_to_kN_str(N):
	k = N/1e3
	return str(round(k, 3))+"kN"
def Pa_to_MPa_str(P):
	M = P/1e6
	return str(round(M, 3))+"MPa"
def coords_str(x,y):
	return "("+str(round(x,3))+","+str(round(y,3))+")"

#prismatic uniform members
class Member:
	def __init__(self, material, xsection, length):
		self.material = material
		self.xsection = xsection
		self.length = length
		self.placed = False
		self.img_ref = None
		self.sup = [None, None]
		self.loads = []
	def __repr__(self):
		return "Member({},{},{})".format(self.material, self.xsection, self.length)
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
	def v_axis(self):
		return numpy.array([self.x1-self.x0, self.y1-self.y0])
	def uv_axis(self):
		#numpy.array(v_axis) / numpy.linalg.norm(v_axis)
		return self.v_axis() / self.length
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
	#Returns the degrees of freedom of the member #(x,y,th)
	def constraints(self):
		c = numpy.array((0,0,0))
		for s in self.sup:
			if not s == None:
				c += s.constraints()
		return c.tolist()
	def d_of_f(self):
		c = self.constraints()
		if c[0] > 1:#Two different x --> x+th
			c[0] -= 1
			c[2] += 1
		if c[1] > 1:#Two different y --> y+th
			c[1] -= 1
			c[2] += 1
		d = numpy.array((1,1,1)) - c
		return d.tolist()
	def gen_report(self, type):
		if type == 0:
			return self.axial_stress()
	def axial_stress(self):
		if max(self.d_of_f()) > 0:
			return "Underconstrained"
		if min(self.d_of_f()) < 0:
			#Overconstrained--not statically determinate
			#I can add a case for this later, with delta=PL/EA
			return "Overconstrained"
		isv = self.is_vert()
		if isv == None:
			return "Not Placed"
		for i,s in enumerate(self.sup):
			if not s==None:
				if s.constraints()[isv]: #0-->x, 1-->y
					fixed_end = i
		#Sort starting at the free end
		if fixed_end == 0:
			self.loads.sort(key=lambda p:-p.ax_dist)
			uv_2free = self.uv_axis()
		elif fixed_end == 1:
			self.loads.sort(key=lambda p:p.ax_dist)
			uv_2free = -self.uv_axis()
		#seg: [[xmin,xmax],internal P]
		p_segments = []
		p_internal = 0
		prev_d = 0
		for p in self.loads:
			#The axial component of the load
			p_ax = numpy.dot(p.get_comp(), uv_2free)*1000#kN
			#d_2free = dist to free end
			if fixed_end == 0:
				d_2free = self.length - p.ax_dist
			if fixed_end == 1:
				d_2free = p.ax_dist
			#Don't define a segment smaller than eps (Happens @ ends)
			if abs(prev_d - d_2free) > eps:
				p_segments.append([[prev_d,d_2free],p_internal])
			p_internal += p_ax
			prev_d = d_2free
		if self.length - prev_d > eps:
			p_segments.append([[prev_d,self.length],p_internal])
		#Now, take those internal forces and convert to stresses and report on it.
		s_segments = []
		xA = self.xarea
		max_s = [(0,0),0]
		min_s = [(0,0),0]
		for ps in p_segments:
			domain,p = ps
			sig = p/xA
			ss = [domain, sig]
			if sig > max_s[1]:
				max_s = ss
			if sig < min_s[1]:
				min_s = ss
			s_segments.append(ss)
		rep_text = "Measuring distance 'x' (in m) from the free end of the member,"
		for ps,ss in zip(p_segments,s_segments):
			domain, p = ps
			domain, s = ss
			rep_text += "\nfor x \u03F5 "+coords_str(domain[0],domain[1])+", "
			rep_text += "tension = "+N_to_kN_str(p)+", "
			rep_text += "\u03C3 = "+Pa_to_MPa_str(s)
		if max_s[1] > 0:
			rep_text += "\nMax Tensile \u03C3 = "+Pa_to_MPa_str(max_s[1])
			rep_text += " @ x \u03F5 "+coords_str(max_s[0][0],max_s[0][1])+", "
		if min_s[1] < 0:
			rep_text += "\nMax Compressive \u03C3 = "+Pa_to_MPa_str(-min_s[1])
			rep_text += " @ x \u03F5 "+coords_str(min_s[0][0],min_s[0][1])+", "
		return rep_text

#This class is really just for reference. I'm not sure if this is the best way to do this.
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
