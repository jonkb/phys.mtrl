import math
import numpy as np
import sympy as sym
from load import *
#from region import Region
eps = 1e-14
default_sigfig = 4
grav = 9.81 #m/s^2

#round x to n significant figures
def sigfig(x, n=default_sigfig):
	if abs(x) < eps:
		return 0
	return round(x, int(-math.log10(abs(x))+n))
#Takes an iterable and rounds those elements that are very close to int
def eps_round(A):
	for i,e in enumerate(A):
		if (abs(e)+eps/2) % 1 < eps:
			A[i] = round(e)
	return A
def sigfig_iter(A, n=default_sigfig):
	for i,e in enumerate(A):
		A[i] = sigfig(e,n)
	return A
def N_to_kN_str(N, n=default_sigfig):
	k = N/1e3
	return str(sigfig(k, n))+"kN"
def Nm_to_kNm_str(Nm, n=default_sigfig):
	km = Nm/1e3
	return str(sigfig(km, n))+"kN-m"
def Pa_to_MPa_str(P, n=default_sigfig):
	M = P/1e6
	return str(sigfig(M, n))+"MPa"
def coords_str(x,y, n=default_sigfig):
	return "("+str(sigfig(x,n))+","+str(sigfig(y,n))+")"

#prismatic uniform members
class Member:
	rep_err = ["Underconstrained", "Overconstrained", "Not Placed", "No Solution Found"]
	def __init__(self, material, xsection, length):
		self.material = material
		self.xsection = xsection
		self.length = length
		self.placed = False
		self.img_ref = None
		self.popups = 0
		self.sup = [None, None]
		self.loads = []
		self.has_weight = False
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
	def xImin(self):
		return self.xsection.Imin
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
		return np.array([self.x1-self.x0, self.y1-self.y0])
	#Unit vector along the axis
	def uv_axis(self):
		#np.array(v_axis) / np.linalg.norm(v_axis)
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
		c = np.array((0,0,0))
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
		d = np.array((1,1,1)) - c
		return d.tolist()
	#Return a distributed load representing the weight of the member
	def weight_dl(self):
		q = grav*self.material["rho"]*self.xarea
		wdl = Distr_Load(None, 0, -q, 0, 0, -q, self.length, self.is_vert())
		return wdl
	#Return a copy of acting loads. Include weight if we're counting weight.
	def my_loads(self):
		loads_c = self.loads.copy()
		if self.has_weight:
			loads_c.append(self.weight_dl())
		return loads_c
	#Returns two functions (of 'd') representing the sum of all distributed loads.
	def sum_my_dl(self):
		fx = sym.Integer(0)
		fy = sym.Integer(0)
		for p in self.my_loads():
			if isinstance(p, Distr_Load):
				(Qx, Qy) = p.to_symf()
				fx += Qx
				fy += Qy
		return (fx, fy)
	#Return shear as a function of 'd'.
	#Convention: shear that would cause clockwise rotation is positive.
	def shear_symf(self):
		(qx, qy) = self.sum_my_dl()
		(s0x, s0y, *_) = self.reactions()
		V = sym.Integer(0)
		d = sym.symbols('d')
		#I should be able to do this with dot products if I ever allow angled members
		if(self.is_vert()):
			V += -s0x
			V += sym.integrate(-qx, d)
			for p in self.my_loads():
				if not isinstance(p, Distr_Load):
					V += -p.xc * sym.Heaviside(d-p.ax_dist, .5)
		else:
			V += s0y
			V += sym.integrate(qy, d)
			for p in self.my_loads():
				if not isinstance(p, Distr_Load):
					V += p.yc * sym.Heaviside(d-p.ax_dist, .5)
		return V
	def moment_symf(self):
		d = sym.symbols('d')
		#M = integrate(V, 0, d)
		M = sym.integrate(self.shear_symf(), d)
		M -= M.subs(d, 0)
		(_, _, s0m, *_) = self.reactions()
		M -= s0m
		return M
	#Do the statics to calculate the reaction forces with two supports
	#Returns a np.array([s0x, s0y, s0m, s1x, s1y, s1m])
	def reactions(self):
		if max(self.d_of_f()) > 0:
			return self.rep_err[0]# "Underconstrained"
		if min(self.d_of_f()) < 0:
			#Overconstrained--not statically determinate
			#I can add a case for this later, with delta=PL/EA
			return self.rep_err[1]# "Overconstrained"
		isv = self.is_vert()
		if isv == None:
			return self.rep_err[2] # "Not Placed"
		#3 Statics equations (left hand only) - sum_x, sum_y, sum_m_0
		st_eq = np.array([[1,0,0,1,0,0], [0,1,0,0,1,0], [0,0,1,0,self.length,1]])
		if isv:
			st_eq[2][3] = -self.length
			st_eq[2][4] = 0
		#3 Support equations, saying which support constraints are zero
		sp_eq = np.array([[0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]])
		
		#They can't both be none, because that would have been underconstrained
		assert self.sup[0] != None or self.sup[1] != None, "141. DoF was calculated wrong."
		if self.sup[0] == None:
			sp_eq = np.array([[1,0,0,0,0,0], [0,1,0,0,0,0], [0,0,1,0,0,0]])
		elif self.sup[1] == None:
			sp_eq = np.array([[0,0,0,1,0,0], [0,0,0,0,1,0], [0,0,0,0,0,1]])
		else:
			n = 0
			sp_c = np.concatenate( (self.sup[0].constraints(), self.sup[1].constraints()) )
			for i, c_i in enumerate(sp_c):
				if c_i == 0:
					assert n < 3, "151. DoF was calculated wrong."
					sp_eq[n][i] = 1
					n += 1
		M_A = np.concatenate( (st_eq,sp_eq) )
		s_px = 0
		s_py = 0
		s_m = 0
		for p in self.my_loads():
			px,py = p.get_comp()
			s_px += px
			s_py += py
			if isv:
				s_m += -px * p.ax_dist
			else:
				s_m += py * p.ax_dist
		M_B = [-s_px, -s_py, -s_m, 0, 0, 0]
		try:
			SOL = np.linalg.solve(M_A, M_B)
		except:
			return self.rep_err[3]
		return sigfig_iter(SOL, 6) #eps_round(SOL)
	def gen_report(self, type):
		if type == 0:
			return self.axial_stress_rep()
		if type == 1:
			return self.axial_buckling_rep()
		if type == 2:
			return self.internal_rep()
	def axial_loads(self):
		if max(self.d_of_f()) > 0:
			return self.rep_err[0]# "Underconstrained"
		if min(self.d_of_f()) < 0:
			#Overconstrained--not statically determinate
			#I can add a case for this later, with delta=PL/EA
			return self.rep_err[1]# "Overconstrained"
		isv = self.is_vert()
		if isv == None:
			return self.rep_err[2] #"Not Placed"
		for i,s in enumerate(self.sup):
			if s != None:
				if s.constraints()[isv]: #0-->x, 1-->y
					fixed_end = i
		loads = self.my_loads()
		#Sort starting at the free end
		if fixed_end == 0:
			loads.sort(key=lambda p:-p.ax_dist)
			uv_2free = self.uv_axis()
		elif fixed_end == 1:
			loads.sort(key=lambda p:p.ax_dist)
			uv_2free = -self.uv_axis()
		#seg: [[xmin,xmax],internal P]
		p_segments = []
		p_internal = 0
		prev_d = 0
		for p in loads:
			#The axial component of the load
			p_ax = np.dot(p.get_comp(), uv_2free)
			#d_2free = dist to free end
			if fixed_end == 0:
				d_2free = self.length - p.ax_dist
			if fixed_end == 1:
				d_2free = p.ax_dist
			#Don't define a segment smaller than eps (Happens @ ends)
			if abs(prev_d - d_2free) > eps:
				p_segments.append( ((prev_d, d_2free), p_internal) )
			p_internal += p_ax
			prev_d = d_2free
		if self.length - prev_d > eps:
			p_segments.append( ((prev_d, self.length), p_internal) )
		return p_segments
	def axial_stress(self):
		p_segments = self.axial_loads()
		if p_segments in self.rep_err:
			return p_segments
		s_segments = []
		xA = self.xarea
		#max_s = ((0,0),0,0)
		#min_s = ((0,0),0,0)
		for ps in p_segments:
			domain,p = ps
			sig = p/xA
			ss = (domain, p, sig)
			#if sig > max_s[1]:
			#	max_s = ss
			#if sig < min_s[1]:
			#	min_s = ss
			s_segments.append(ss)
		return s_segments# (p_segments, s_segments, min_s, max_s)
	def axial_strain(self):
		s_segments = self.axial_stress()
		if s_segments in self.rep_err:
			return s_segments
		eps_segments = []
		for s in s_segments:
			domain, p, sig = s
			eps = sig / self.E
			eps_segments.append( (domain, p, sig, eps) )
		return eps_segments
	def axial_stress_rep(self):
		axeps = self.axial_strain()
		if axeps in self.rep_err:
			return axeps
		max_s = ((0,0), 0, 0, 0)
		min_s = ((0,0), 0, 0, 0)
		for seg in axeps:
			#domain, p, sig, eps = seg
			if seg[1] > max_s[1]:
				max_s = seg
			if seg[1] < min_s[1]:
				min_s = seg
		rep_text = "Measuring distance 'x' (in m) from the free end of the member,"
		for seg in axeps:
			domain, p, sig, eps = seg
			rep_text += "\nfor x \u03F5 "+coords_str(domain[0],domain[1])+", "
			rep_text += "tension = "+N_to_kN_str(p)+ ", "
			rep_text += "\u03C3 = "+Pa_to_MPa_str(sig)+ ", "
			rep_text += "\u03B5 = "+str(sigfig(eps,4))
		if max_s[1] > 0:
			rep_text += "\nMax Tensile \u03C3 = "+Pa_to_MPa_str(max_s[2])
			rep_text += " @ x \u03F5 "+coords_str(max_s[0][0],max_s[0][1])+", "
			typer = round(max_s[2] / self.material["sig_y"] * 100, 2)
			tuper = round(max_s[2] / self.material["sig_u"] * 100, 2)
			rep_text += "\nWhich is "+str(typer)+"% of the yield stress and "
			rep_text += str(tuper)+"% of the ultimate stress."
		if min_s[1] < 0:
			rep_text += "\nMax Compressive \u03C3 = "+Pa_to_MPa_str(-min_s[2])
			rep_text += " @ x \u03F5 "+coords_str(min_s[0][0],min_s[0][1])+", "
			cyper = round(-min_s[2] / self.material["sig_y"] * 100, 2)
			cuper = round(-min_s[2] / self.material["sig_u"] * 100, 2)
			rep_text += "\nWhich is "+str(cyper)+"% of the yield stress and "
			rep_text += str(cuper)+"% of the ultimate stress."
			rep_text += "\nThis evaluation does not consider buckling."
		return rep_text
	def axial_buckling_rep(self):
		p_seg = self.axial_loads()
		if p_seg in self.rep_err:
			return p_seg
		if len(p_seg) > 1:
			rep_text = "Sorry, I'm not exactly sure how to deal with "
			rep_text += "\nbuckling in a column with multiple loads."
			#I could give it a shot, though.
			return rep_text
		#I need to check the supports to find out the end conditions
		rep_text = "Sorry, this is super basic right now."
		rep_text += "\nI can only deal with fixed-free columns."
		Pmax = -p_seg[0][1]
		if Pmax < 0:
			rep_text += "\nThere is no compressive load on this member, so it will not buckle."
			return rep_text
		Pcr = math.pi**2 * self.E * self.xImin / (4 * self.length**2)
		Pper = round(Pmax / Pcr * 100, 2)
		rep_text += "\nThe load applied to this member P = " + N_to_kN_str(Pmax)
		rep_text += "\nThis is " + str(Pper) + "% of the critical load Pcr = " + N_to_kN_str(Pcr)
		if Pmax < Pcr:
			rep_text += "\nSo the member is stable"
		if Pmax >= Pcr:
			rep_text += "\nSo the member is unstable"
		return rep_text
	#Report on internal forces, sheer, and moment
	def internal_rep(self):
		(s0x, s0y, s0m, s1x, s1y, s1m) = self.reactions()
		rep_text = "Support reactions (Rx, Ry, M)"
		rep_text += "\nSupport 0: " + "("+N_to_kN_str(s0x)+","
		rep_text += N_to_kN_str(s0y)+","+Nm_to_kNm_str(s0m)+")"
		rep_text += "\nSupport 1: " + "("+N_to_kN_str(s1x)+","
		rep_text += N_to_kN_str(s1y)+","+Nm_to_kNm_str(s1m)+")"
		rep_text += "\nMeasuring 'd' (in m) from end zero (left or bottom) of the member,"
		rep_text += "\nFor internal tension, see \"axial stress report\""
		rep_text += "\nV(d) = "+str(self.shear_symf())
		rep_text += "\nM(d) = "+str(self.moment_symf())
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
