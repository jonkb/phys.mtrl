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
	axis_resolution = 200
	y1_resolution = 50
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
		d = sym.symbols('d')
		#M = integrate(V, 0, d)
		M = sym.integrate(self.shear_symf(), d)
		M -= M.subs(d, 0)
		reactions = self.reactions()
		if isinstance(reactions, str):
		#if reactions in self.rep_err:
			return reactions
		(_, _, s0m, *_) = reactions
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
	#Return two vectors, along with their magnitude: sig_max, smv, tau_max, tmv
	@staticmethod
	def mohr_trsfm(sigx, tau):
		th_a = np.arctan(2*tau/sigx)/2
		#th_a = np.arctan2(2*tau,sigx)/2   --   This would be just to avoid /0 Repercussions?
		th_b = th_a + math.pi/2
		sig_a = sigx/2*(1+np.cos(2*th_a)) + tau*np.sin(2*th_a)
		sig_b = sigx/2*(1+np.cos(2*th_b)) + tau*np.sin(2*th_b)
		#I converted the ternaries to array-compatible logic
		b_gt_a = sig_b > sig_a
		b_agt_a = abs(sig_b) > abs(sig_a)
		#th_smax = th_b if abs(sig_b) > abs(sig_a) else th_a
		th_smax = th_b*b_agt_a + th_a*(1 - b_agt_a)
		#sig_max = sig_b if abs(sig_b) > abs(sig_a) else sig_a
		sig_max = sig_b*b_agt_a + sig_a*(1 - b_agt_a)
		#th_s1 = th_b if sig_b > sig_a else th_a #alg. max.
		th_s1 = th_b*b_gt_a + th_a*(1 - b_gt_a)
		th_tmax = th_s1 - math.pi/4
		tau_max = np.sqrt((sigx/2)**2 + tau**2)
		#return (sig_max, th_smax), (tau_max, th_tmax)
		sig_max_v = (sig_max*np.cos(th_smax), sig_max*np.sin(th_smax))
		tau_max_v = (tau_max*np.cos(th_tmax), tau_max*np.sin(th_tmax))
		return sig_max, sig_max_v, tau_max, tau_max_v
	#Return report text (& any figures) for each type of evaluation
	#Be sure to make this line up with the dictionary in Lab
	def gen_report(self, type):
		if type == 0:
			return self.axial_stress_rep()
		elif type == 1:
			return self.axial_buckling_rep()
		elif type == 2:
			return self.VM_rep()
		elif type == 3:
			return self.sig_tau_rep()
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
		Sig_P = self.axial_loads_sym() / self.xarea
		Sig_M = - self.moment_symf() * h / self.xIx
		return Sig_P + Sig_M
	#Return internal shear stress due to bending as a function of 'd' and 'h'.
	#Returns (tau, h domain)
	def tau_sym(self):
		h = sym.symbols('h')
		tau = -self.shear_symf() * self.xsection.Q_div_Ib(h)
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
		rep_text += "\nMeasuring 'd' (in m) from end zero (left or bottom) of the member,"
		rep_text += "\nFor internal tension, see \"axial stress report\""
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
		dax = np.linspace(0, self.length, self.axis_resolution);
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
				Vax = Vc*np.ones(self.axis_resolution)
			Mc = Mf(0)
			for xi in dax:
				if Mf(xi) != Mc:
					Mc = "error"
					break
			if Mc != "error":
				Max = Mc*np.ones(self.axis_resolution)
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
		rep_text = "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		rep_text += "\nand 'h' (in m) from the neutral axis of bending,"
		d, h = sym.symbols("d h")
		sig = self.axial_stress_sym()
		y1min, y1max = self.xsection.y1_domain()
		tau, h_dom = self.tau_sym()
		sigf = sym.lambdify((d,h), sig/1e6, "numpy")
		tauf = sym.lambdify((d,h), tau/1e6, "numpy")
		d_ls = np.linspace(0, self.length, self.axis_resolution)
		h_sig_ls = np.linspace(y1min, y1max, self.y1_resolution)
		h_tau_ls = np.linspace(*h_dom, self.y1_resolution)
		D, H_s = np.meshgrid(d_ls, h_sig_ls)
		_, H_t = np.meshgrid(d_ls, h_tau_ls)
		SIG = sigf(D, H_s)
		sig_max = np.max(SIG)
		sig_min = np.min(SIG)
		sig_rng = max(abs(sig_min), abs(sig_max))
		if sig_max > 0:
			rep_text += "\nMax Tensile Axial Stress = " + str(sigfig(sig_max)) + " MPa"
			smax_coords = np.where(SIG == sig_max)
			smax_d = smax_coords[1][0] * self.length/(self.axis_resolution-1)
			smax_h = y1min + smax_coords[0][0] * (y1max-y1min) / (self.y1_resolution-1)
			rep_text += " at d="+str(sigfig(smax_d))+"m, h="+str(sigfig(smax_h*1e3))+"mm"
		if sig_min < 0:
			rep_text += "\nMax Compressive Axial Stress = " + str(sigfig(abs(sig_min))) + " MPa"
			smin_coords = np.where(SIG == sig_min)
			smin_d = smin_coords[1][0] * self.length/(self.axis_resolution-1)
			smin_h = y1min + smin_coords[0][0] * (y1max-y1min) / (self.y1_resolution-1)
			rep_text += " at d="+str(sigfig(smin_d))+"m, h="+str(sigfig(smin_h*1e3))+"mm"
		TAU = tauf(D, H_t)
		tau_max = np.max(TAU)
		tau_min = np.min(TAU)
		tau_rng = max(abs(tau_min), abs(tau_max))
		if tau_max > 0:
			rep_text += "\nMax Positive Sheer Stress = " + str(sigfig(tau_max)) + " MPa"
			tmax_coords = np.where(TAU == tau_max)
			tmax_d = tmax_coords[1][0] * self.length/(self.axis_resolution-1)
			tmax_h = h_dom[0] + tmax_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.y1_resolution-1)
			rep_text += " at d="+str(sigfig(tmax_d))+"m, h="+str(sigfig(tmax_h*1e3))+"mm"
		if tau_min < 0:
			rep_text += "\nMax Negative Sheer Stress = " + str(sigfig(abs(tau_min))) + " MPa"
			tmin_coords = np.where(TAU == tau_min)
			tmin_d = tmin_coords[1][0] * self.length/(self.axis_resolution-1)
			tmin_h = h_dom[0] + tmin_coords[0][0] * (h_dom[1]-h_dom[0]) / (self.y1_resolution-1)
			rep_text += " at d="+str(sigfig(tmin_d))+"m, h="+str(sigfig(tmin_h*1e3))+"mm"
		
		fig, ((sp1, sp2), (sp3, sp4)) = plt.subplots(2, 2, sharex=True, figsize=(11,5))
		sp1.set_title("Axial Stress \u03C3 (MPa)")
		sp1.set(ylabel="Height h (mm)")
		im1 = sp1.imshow(SIG, cmap=plt.cm.RdBu, interpolation="bilinear", aspect="auto", 
			extent=[0, self.length, y1min*1e3, y1max*1e3], vmin=-sig_rng, vmax=sig_rng, origin="lower")
			#Using vmin&max normalizes zero on the colormap
		fig.colorbar(im1, ax=sp1)
		sp3.set_title("Shear Stress \u03C4 (MPa)")
		sp3.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		im2 = sp3.imshow(TAU, cmap=plt.cm.BrBG, interpolation="bilinear", aspect="auto", 
			extent=[0, self.length, h_dom[0]*1e3, h_dom[1]*1e3], vmin=-tau_rng, vmax=tau_rng, origin="lower")
		fig.colorbar(im2, ax=sp3)
		
		#Quiver plots for max shear & tension
		d_qls = np.linspace(0, self.length, 48)
		h_qls = np.linspace(*h_dom, 24)
		Dq, Hq = np.meshgrid(d_qls, h_qls)
		Sq = sigf(Dq, Hq)
		Tq = tauf(Dq, Hq)
		S1m, (S1x, S1y), T1m, (T1x, T1y) = self.mohr_trsfm(Sq, Tq)
		sig_max = np.max(S1m)
		sig_min = np.min(S1m)
		sig_rng = max(abs(sig_min), abs(sig_max))
		tau_max = np.max(T1m)
		tau_min = np.min(T1m)
		tau_rng = max(abs(tau_min), abs(tau_max))
		sp2.set_title("Principal Stress \u03C3_p (MPa)")
		q1 = sp2.quiver(Dq, Hq*1e3, S1x, S1y, S1m, cmap=plt.cm.RdBu, clim=(-sig_rng, sig_rng),
			pivot="mid", headaxislength=0, headlength = 0, headwidth=1, width=.009)
		fig.colorbar(q1, ax=sp2)
		sp4.set(xlabel="Axial Distance d (m)")
		sp4.set_title("Principal Shear Stress \u03C4_p (MPa)")
		q2 = sp4.quiver(Dq, Hq*1e3, T1x, T1y, T1m, cmap=plt.cm.BrBG, clim=(-tau_rng, tau_rng),
			pivot="mid", headaxislength=0, headlength = 0, headwidth=1, width=.009)
		fig.colorbar(q2, ax=sp4)
		
		plt.subplots_adjust(.07, .1, .99, .93, wspace=.09)
		
		return (rep_text, fig)

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
