import math
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
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
	d_resolution = 100
	h_resolution = 50
	dq_resolution = 50
	hq_resolution = 25
	#W, S, E, N
	plot_margins = (.11, .19, 1, .9)
	rep_err = ["Underconstrained", "Overconstrained", "Not Placed", "No Solution Found"]
	#Names of each Evaluation Report
	eval_names = {
		0: "Mass Properties",
		1: "Axial and Shear Stress",
		2: "Shear and Moment",
		3: "Euler Buckling"
	}
	
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
		
		self.eval_reports = {
			0: self.mass_prop_rep,
			1: self.sig_tau_rep,
			2: self.VM_rep,
			3: self.axial_buckling_rep
		}
	def __repr__(self):
		return "Member({},{},{})".format(self.material, self.xsection, self.length)
	def __str__(self):
		s = self.material["name"] + " member with cross section=("
		s += str(self.xsection) + "), and length=" + str(self.length)
		return s
	def to_xml(self):
		data = """
		<mem>
			<def>
				<material>"""+str(self.material["vname"])+"</material>"
		data += self.xsection.to_xml()
		data += """
				<length>"""+str(self.length)+"""</length>
			</def>
			<place>
				<x0>"""+str(self.x0)+"""</x0>
				<y0>"""+str(self.y0)+"""</y0>
				<vh>"""+self.VH_char()+"""</vh>
			</place>
			<!-- Sup Type= 0:Fixed, 1:Pin, 2: Slot(x), 3: Slot(y) -->"""
		for i, sup in enumerate(self.sup):
			if sup == None:
				continue
			data += """
			<sup type="{}" end="{}"/>""".format(sup.stype, i)
		data += """
			<!-- Ld Type= 0:Point, 1:Distributed -->"""
		for p in self.loads:
			data += p.to_xml()
		data += """
		</mem>"""
		return data
	
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
	def volume(self):
		return self.xarea * self.length
	@property
	def mass(self):
		return self.volume * self.material["rho"]
	@property
	def weight(self):
		return self.mass * grav
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
	def VH_char(self):
		if self.is_vert():
			return "V"
		else:
			return "H"
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
		reactions = self.reactions()
		if isinstance(reactions, str):
		#if reactions in self.rep_err:
			return reactions
		(s0x, s0y, *_) = reactions
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
		#Rewriting as piecewise helps it to integrate the Heaviside
		return V.rewrite(sym.Piecewise).doit()
	def moment_symf(self):
		reactions = self.reactions()
		if isinstance(reactions, str):
		#if reactions in self.rep_err:
			return reactions
		(_, _, s0m, *_) = reactions
		d = sym.symbols('d')
		#M = integrate(V, 0, d)
		M = sym.integrate(self.shear_symf(), d)
		M -= M.subs(d, 0)
		M -= s0m
		return M.rewrite(sym.Piecewise).doit()
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
	#Do mohr's circle transformations on an element with sigma_x and tau_xy (assumes no sig_y)
	#Return two vectors, each in polar and cartesian: smvp, smvc, tmvp, tmvc
	#fm = module with appropriate functions. Default numpy.
	@staticmethod
	def mohr_trsfm(sigx, tau, fm=np):
		convert_to_neg = False
		try:
			#numpy uses "arctan"
			atan = fm.arctan
		except:
			#sympy uses "atan"
			atan = fm.atan
		th_a = atan(2*tau/sigx)/2
		#th_a = np.arctan2(2*tau,sigx)/2   --   This would be just to avoid /0 Repercussions?
		th_b = th_a + math.pi/2
		sig_a = sigx/2*(1+fm.cos(2*th_a)) + tau*fm.sin(2*th_a)
		sig_b = sigx/2*(1+fm.cos(2*th_b)) + tau*fm.sin(2*th_b)
		#I converted the ternaries to array-compatible logic
		b_gt_a = sig_b > sig_a
		b_agt_a = abs(sig_b) > abs(sig_a)
		#th_smax = th_b if abs(sig_b) > abs(sig_a) else th_a
		th_smax = th_b*b_agt_a + th_a*(1 - b_agt_a)
		#sig_max = sig_b if abs(sig_b) > abs(sig_a) else sig_a
		sig_max = sig_b*b_agt_a + sig_a*(1 - b_agt_a)
		#th_s1 = th_b if sig_b > sig_a else th_a #alg. max.
		th_s1 = th_b*b_gt_a + th_a*(1 - b_gt_a)
		th_t1 = th_s1 - math.pi/4 #Alg. max tau
		th_t2 = th_t1 + math.pi/2
		#Force the angle to be more horizontal so you get negative instead of perpindicular
		if convert_to_neg:
			t1_horiz = math.pi/2 - abs(th_t1%math.pi - math.pi/2)
			t2_horiz = math.pi/2 - abs(th_t2%math.pi - math.pi/2)
			#is t1 the more horizontal one?
			t1_is_h = t2_horiz >= t1_horiz
			th_tmax = th_t1*t1_is_h + th_t2*(1-t1_is_h)
			tau_max = -.5*sigx*fm.sin(2*th_tmax) + tau*fm.cos(2*th_tmax)
		else: #Just accept the alg. max. tau
			th_tmax = th_t1
			tau_max = fm.sqrt((sigx/2)**2 + tau**2) #alg. max
		#return (sig_max, th_smax), (tau_max, th_tmax)
		sig_max_vp = (sig_max, th_smax)
		sig_max_vc = (sig_max*fm.cos(th_smax), sig_max*fm.sin(th_smax))
		tau_max_vp = (tau_max, th_tmax)
		tau_max_vc = (tau_max*fm.cos(th_tmax), tau_max*fm.sin(th_tmax))
		return sig_max_vp, sig_max_vc, tau_max_vp, tau_max_vc
	#Return report text (& any figures) for each type of evaluation
	#Be sure to make this line up with the dictionary in Lab
	def gen_report(self, type):
		return self.eval_reports[type]()
	def mass_prop_rep(self):
		rep_text = "cross sectional area A = "+str(sigfig(self.xarea))+" m^2"
		rep_text += "\n\tsecond moment of area about horizontal axis I_x = "+str(sigfig(self.xIx))+" m^4"
		rep_text += "\n\tsecond moment of area about vertical axis I_y = "+str(sigfig(self.xIy))+" m^4"
		rep_text += "\nvolume V = "+str(sigfig(self.volume))+" m^3"
		rep_text += "\nmass m = "+str(sigfig(self.mass))+" kg"
		rep_text += "\nweight W = "+str(sigfig(self.weight))+" N"
		rep_text += "\n\nProperties of "+self.material["name"]+":"
		rep_text += "\n\tmass density \u03C1 = "+str(self.material["rho"])+" kg/m^3"
		rep_text += "\n\tYoung's modulus E = "+str(self.material["E"]/1e9)+" GPa"
		rep_text += "\n\tyield stress \u03C3_y = "+str(self.material["sig_y"]/1e6)+" MPa"
		rep_text += "\n\tultimate stress \u03C3_y = "+str(self.material["sig_u"]/1e6)+" MPa"
		return rep_text
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
	#Return internal axial loads as a function of 'd'.
	#Convention: tension is positive.
	def axial_loads_sym(self):
		uv_ax = self.uv_axis()
		(qx, qy) = self.sum_my_dl()
		reactions = self.reactions()
		if isinstance(reactions, str):
			return reactions
		(s0x, s0y, *_) = reactions
		P = sym.Float(np.dot([s0x, s0y], -uv_ax))
		d = sym.symbols('d')
		P -= sym.integrate(qx*uv_ax[0], d) + sym.integrate(qy*uv_ax[1], d)
		for p in self.my_loads():
			#The distributed loads were all counted already
			if not isinstance(p, Distr_Load):
				p_ax = np.dot(uv_ax, p.get_comp())
				P -= p_ax*sym.Heaviside(d-p.ax_dist, .5)
		return P.rewrite(sym.Piecewise).doit()
	#Return internal axial stress as a function of 'd' and 'h'.
	#Counts axial loads as well as bending stress
	def axial_stress_sym(self):
		#sig_x = Px/A - My/I
		h = sym.symbols('h')
		P = self.axial_loads_sym()
		if isinstance(P, str):
			return P
		M = self.moment_symf()
		if isinstance(M, str):
			return M
		Sig_P = P / self.xarea
		Sig_M = - M * h / self.xIx
		return Sig_P + Sig_M
	#Return internal shear stress due to bending as a function of 'd' and 'h'.
	#Returns (tau, h domain)
	def tau_sym(self):
		h = sym.symbols('h')
		V = self.shear_symf()
		if isinstance(V, str):
			return V
		tau = -V * self.xsection.Q_div_Ib(h)
		return (tau, self.xsection.Q_domain())
	def axial_stress(self):
		p_segments = self.axial_loads()
		if isinstance(p_segments, str):
		#if p_segments in self.rep_err:
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
		if isinstance(s_segments, str):
		#if s_segments in self.rep_err:
			return s_segments
		eps_segments = []
		for s in s_segments:
			domain, p, sig = s
			eps = sig / self.E
			eps_segments.append( (domain, p, sig, eps) )
		return eps_segments
	#Old, non-sym function
	def axial_stress_rep(self):
		axeps = self.axial_strain()
		if isinstance(axeps, str):
		#if axeps in self.rep_err:
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
		if isinstance(p_seg, str):
		#if p_seg in self.rep_err:
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
	#Report on shear and moment
	#Returns text and a pyplot figure
	def VM_rep(self):
		reactions = self.reactions()
		if isinstance(reactions, str):
		#if reactions in self.rep_err:
			return reactions
		(s0x, s0y, s0m, s1x, s1y, s1m) = reactions
		rep_text = "Support reactions (Rx, Ry, M)"
		rep_text += "\nSupport 0: " + "("+N_to_kN_str(s0x)+","
		rep_text += N_to_kN_str(s0y)+","+Nm_to_kNm_str(s0m)+")"
		rep_text += "\nSupport 1: " + "("+N_to_kN_str(s1x)+","
		rep_text += N_to_kN_str(s1y)+","+Nm_to_kNm_str(s1m)+")"
		rep_text += "\n'd' is measured (in m) from end \"zero\" (W or S) of the member."
		rep_text += "\nFor internal tension and stresses, see \"Axial and Shear Stress\""
		V = self.shear_symf()
		if V in self.rep_err:
			return V
		M = self.moment_symf()
		if M in self.rep_err:
			return M
		#Make plots for V and M
		d = sym.symbols("d")
		Vf = sym.lambdify(d, V/1000)
		Mf = sym.lambdify(d, M/1000)
		dax = np.linspace(0, self.length, self.d_resolution);
		Vax = Vf(dax)
		Max = Mf(dax)
		try:
			assert len(dax) == len(Vax)
			assert len(dax) == len(Max)
		except:
			#Check to see if they're constant functions
			Vc = Vf(0)
			for xi in dax:
				if Vf(xi) != Vc:
					Vc = "error"
					break
			if Vc != "error":
				Vax = Vc*np.ones(self.d_resolution)
			Mc = Mf(0)
			for xi in dax:
				if Mf(xi) != Mc:
					Mc = "error"
					break
			if Mc != "error":
				Max = Mc*np.ones(self.d_resolution)
		fig, (sp1, sp2) = plt.subplots(2, sharex=True)
		sp1.plot(dax, Vax, color="blue")
		sp1.grid(True)
		sp1.set_title("Shear V (kN)")
		sp2.plot(dax, Max, color="green")
		#Pf = sym.lambdify(d, self.axial_loads_sym()/1000)
		#Pax = Pf(dax)
		#sp2.plot(dax, Pax, color="green")
		sp2.grid(True)
		sp2.set_title("Moment M (kN-m)")
		sp2.set(xlabel="Axial Distance d (m)")
		return (rep_text, fig)
	def sig_tau_rep(self):
		reps = []
		d, h = sym.symbols("d h")
		sig = self.axial_stress_sym()
		if isinstance(sig, str):
			return sig #ERROR
		y1min, y1max = self.xsection.y1_domain()
		tau_res = self.tau_sym()
		if isinstance(tau_res, str):
			return tau_res
		tau, h_dom = tau_res
		sigf = sym.lambdify((d,h), sig/1e6, "numpy")
		tauf = sym.lambdify((d,h), tau/1e6, "numpy")
		#TEMP
		(smvp, smvc, tmvp, tmvc) = self.mohr_trsfm(sig, tau, fm=sym)
		smx, smy = smvc
		#\TEMP
		d_ls = np.linspace(0, self.length, self.d_resolution)
		h_sig_ls = np.linspace(y1min, y1max, self.h_resolution)
		h_tau_ls = np.linspace(*h_dom, self.h_resolution)
		D, H_s = np.meshgrid(d_ls, h_sig_ls)
		_, H_t = np.meshgrid(d_ls, h_tau_ls)
		SIG = sigf(D, H_s)
		sig_max = np.max(SIG)
		sig_min = np.min(SIG)
		sig_rng = max(abs(sig_min), abs(sig_max))
		s1_rep_text = "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		s1_rep_text += " and 'h' (in mm) from the neutral axis of bending,"
		if True: #OLD
			if sig_max > 0:
				s1_rep_text += "\nMax Tensile Axial Stress = " + str(sigfig(sig_max)) + " MPa"
				smax_coords = np.where(SIG == sig_max)
				smax_d = smax_coords[1][0] * self.length/(self.d_resolution-1)
				smax_h = y1min + smax_coords[0][0] * (y1max-y1min) / (self.h_resolution-1)
				s1_rep_text += " at d="+str(sigfig(smax_d))+"m, h="+str(sigfig(smax_h*1e3))+"mm"
			if sig_min < 0:
				s1_rep_text += "\nMax Compressive Axial Stress = " + str(sigfig(abs(sig_min))) + " MPa"
				smin_coords = np.where(SIG == sig_min)
				smin_d = smin_coords[1][0] * self.length/(self.d_resolution-1)
				smin_h = y1min + smin_coords[0][0] * (y1max-y1min) / (self.h_resolution-1)
				s1_rep_text += " at d="+str(sigfig(smin_d))+"m, h="+str(sigfig(smin_h*1e3))+"mm"
		TAU = tauf(D, H_t)
		tau_max = np.max(TAU)
		tau_min = np.min(TAU)
		tau_rng = max(abs(tau_min), abs(tau_max))
		t1_rep_text = "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		t1_rep_text += " and 'h' (in mm) from the neutral axis of bending,"
		if True:
			if tau_max > 0:
				t1_rep_text += "\nMax Positive Sheer Stress = " + str(sigfig(tau_max)) + " MPa"
				tmax_coords = np.where(TAU == tau_max)
				tmax_d = tmax_coords[1][0] * self.length/(self.d_resolution-1)
				tmax_h = h_dom[0] + tmax_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.h_resolution-1)
				t1_rep_text += " at d="+str(sigfig(tmax_d))+"m, h="+str(sigfig(tmax_h*1e3))+"mm"
			if tau_min < 0:
				t1_rep_text += "\nMax Negative Sheer Stress = " + str(sigfig(abs(tau_min))) + " MPa"
				tmin_coords = np.where(TAU == tau_min)
				tmin_d = tmin_coords[1][0] * self.length/(self.d_resolution-1)
				tmin_h = h_dom[0] + tmin_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.h_resolution-1)
				t1_rep_text += " at d="+str(sigfig(tmin_d))+"m, h="+str(sigfig(tmin_h*1e3))+"mm"
		
		#Plot of in-plane axial stress
		fig1, ax1 = plt.subplots(figsize=(7,3))
		ax1.set_title("Axial Stress \u03C3 (MPa)")
		ax1.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		im1 = ax1.imshow(SIG, cmap=plt.cm.RdBu, interpolation="bilinear", aspect="auto", 
			extent=[0, self.length, y1min*1e3, y1max*1e3], vmin=-sig_rng, vmax=sig_rng, origin="lower")
		fig1.colorbar(im1, ax=ax1)
		plt.subplots_adjust(*self.plot_margins, wspace=.09)
		reps.append(("In-plane \u03C3", s1_rep_text, fig1))
		#Plot of in-plane shear stress
		fig2, ax2 = plt.subplots(figsize=(7,3))
		ax2.set_title("Shear Stress \u03C4 (MPa)")
		ax2.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		im2 = ax2.imshow(TAU, cmap=plt.cm.BrBG, interpolation="bilinear", aspect="auto", 
			extent=[0, self.length, h_dom[0]*1e3, h_dom[1]*1e3], vmin=-tau_rng, vmax=tau_rng, origin="lower")
		fig2.colorbar(im2, ax=ax2)
		plt.subplots_adjust(*self.plot_margins, wspace=.09)
		reps.append(("In-plane \u03C4", t1_rep_text, fig2))
		
		#Quiver plots for max shear & tension
		d_qls = np.linspace(0, self.length, self.dq_resolution)
		h_qls = np.linspace(*h_dom, self.hq_resolution)
		Dq, Hq = np.meshgrid(d_qls, h_qls)
		Sq = sigf(Dq, Hq)
		Tq = tauf(Dq, Hq)
		(S1m, S1th), (S1x, S1y), (T1m, T1th), (T1x, T1y) = self.mohr_trsfm(Sq, Tq)
		#Swap these for gradients
		sig_max = np.max(S1m)
		sig_min = np.min(S1m)
		sig_rng = max(abs(sig_min), abs(sig_max))
		tau_max = np.max(T1m)
		tau_min = np.min(T1m)
		tau_rng = max(abs(tau_min), abs(tau_max))
		s2_rep_text = "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		s2_rep_text += " and 'h' (in mm) from the neutral axis of bending,"
		if sig_max > 0:
			s2_rep_text += "\nMax Tensile Stress = " + str(sigfig(sig_max)) + " MPa"
			smax_coords = np.where(S1m == sig_max)
			smax_d = smax_coords[1][0] * self.length/(self.dq_resolution-1)
			smax_h = h_dom[0] + smax_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.hq_resolution-1)
			smax_th = S1th[smax_coords[0][0], smax_coords[1][0]]
			s2_rep_text += " at d="+str(sigfig(smax_d))+"m, h="+str(sigfig(smax_h*1e3))+"mm"
			s2_rep_text += ", \u03B8="+str(sigfig(math.degrees(smax_th)))+"\u00B0"
		if sig_min < 0:
			s2_rep_text += "\nMax Compressive Stress = " + str(sigfig(-sig_min)) + " MPa"
			smin_coords = np.where(S1m == sig_min)
			smin_d = smin_coords[1][0] * self.length/(self.dq_resolution-1)
			smin_h = h_dom[0] + smin_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.hq_resolution-1)
			smin_th = S1th[smin_coords[0][0], smin_coords[1][0]]
			s2_rep_text += " at d="+str(sigfig(smin_d))+"m, h="+str(sigfig(smin_h*1e3))+"mm"
			s2_rep_text += ", \u03B8="+str(sigfig(math.degrees(smin_th)))+"\u00B0"
		t2_rep_text = "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		t2_rep_text += " and 'h' (in mm) from the neutral axis of bending,"
		if tau_rng > 0:
			t2_rep_text += "\nMax Sheer Stress = " + str(sigfig(tau_rng)) + " MPa"
			tmax_coords = np.where(abs(T1m) == tau_rng)
			tmax_d = tmax_coords[1][0] * self.length/(self.dq_resolution-1)
			tmax_h = h_dom[0] + tmax_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.hq_resolution-1)
			tmax_th = T1th[tmax_coords[0][0], tmax_coords[1][0]]
			t2_rep_text += " at d="+str(sigfig(tmax_d))+"m, h="+str(sigfig(tmax_h*1e3))+"mm"
			t2_rep_text += ", \u03B8="+str(sigfig(math.degrees(tmax_th)))+"\u00B0"
		
		#Plot of principal (out-of-plane) stress
		fig3, ax3 = plt.subplots(figsize=(7,3))
		ax3.set_title("Principal Stress \u03C3_p (MPa)")
		ax3.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		q1 = ax3.quiver(Dq, Hq*1e3, S1x, S1y, S1m, cmap=plt.cm.RdBu, clim=(-sig_rng, sig_rng),
			pivot="mid", headaxislength=0, headlength = 0, headwidth=1, width=.009)
		fig3.colorbar(q1, ax=ax3)
		plt.subplots_adjust(*self.plot_margins, wspace=.09)
		reps.append(("Out-of-plane \u03C3", s2_rep_text, fig3))
		#Plot of principal (out-of-plane) shear stress
		fig4, ax4 = plt.subplots(figsize=(7,3))
		ax4.set_title("Principal Shear Stress \u03C4_p (MPa)")
		ax4.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		q2 = ax4.quiver(Dq, Hq*1e3, T1x, T1y, T1m, cmap=plt.cm.BrBG, clim=(-tau_rng, tau_rng),
			pivot="mid", headaxislength=0, headlength = 0, headwidth=1, width=.009)
		fig4.colorbar(q2, ax=ax4)
		plt.subplots_adjust(*self.plot_margins, wspace=.09)
		reps.append(("Out-of-plane \u03C4", t2_rep_text, fig4))
		
		return reps

#This class is really just for reference. I'm not sure if this is the best way to do this.
class Materials:
	materials = ("steel", "aluminum", "pla", "cell_pvc", "oak")
	#Used values for ASTM-A514 steel
	steel = {
		"vname": "steel", #Variable name
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
		"vname": "aluminum",
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
	#PLA Plastic (average values from http://www.matweb.com/search/DataSheet.aspx?MatGUID=ab96a4c0655c4018a8785ac4031b9278&ckck=1)
	pla = {
		"vname": "pla",
		"name": "PLA Plastic",
		"color": "#65b565",
		#Mass density
		"rho": 1290, #1.00 - 2.47 g/cc
		#Young's modulus
		"E": 2.91e9, #0.0850 - 13.8 GPa	
		#yield stress (TENSILE)
		"sig_y": 38.0e6, #2.00 - 103 MPa	
		#ultimate stress (TENSILE)
		"sig_u": 47.2e6 #14.0 - 117 MPa	
	}
	#cellular PVC board (https://azekexteriors.com/docs/technical/azek-csi-format-spec-1-22-19v1-1.pdf)
	cell_pvc = {
		"vname": "cell_pvc",
		"name": "Cellular PVC",
		"color": "#f0f0f0",
		"rho": 550,
		"E": 992845050, #.99 GPa
		"sig_y": 0, #Unknown
		"sig_u": 15.56e6
	}
	#Oak wood (http://www.matweb.com/search/datasheet_print.aspx?matguid=3a971164050b4313930591eed2539366)
	oak = {
		"vname": "oak",
		"name": "Oak Wood",
		"color": "#d0ac7c",
		"rho": 630,
		"E": 12.56e9,
		"sig_y": 47.0e6, #Compressive, perpindicular to grain
		"sig_u": 5.5e6 #Tensile
	}
