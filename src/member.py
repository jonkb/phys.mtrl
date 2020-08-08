import math
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from load import *
import joint as jt
#from region import Region
from math_util import *



#WARNING: reactions() and axial_loads() are broken under the new system!!!



#prismatic uniform members
class Member:
	d_resolution = 100
	h_resolution = 50
	dq_resolution = 50
	hq_resolution = 25
	#W, S, E, N
	plot_margins = (.11, .19, .95, .9)
	implot_margins = (.11, .19, 1, .9)
	rep_err = ["Underconstrained", "Overconstrained", "Not Placed", "No Solution Found"]
	rep_wrng = ["Warning: because of the shape of the cross section, tau is only known \
at the neutral axis, so the max and min values presented here may be inaccurate.\n",
		"Warning: because of the shape of the cross section, the calculation for tau is not \
valid over the whole height, so the max and min values presented here may be inaccurate.\n"]
	
	#Names of each Evaluation Report
	#Match this to the dictionary in __init__
	eval_names = {
		0: "Mass Properties",
		1: "Static Equilibrium",
		2: "Axial Stress and Strain",
		3: "Shear and Moment",
		4: "Axial and Shear Stress",
		5: "Euler Buckling"
	}
	
	def __init__(self, material, xsection, length):
		self.material = material
		self.xsection = xsection
		self.length = length
		self.placed = False
		self.img_ref = None
		self.popups = 0
		self.supports = []
		self.joints = []
		self.loads = []
		self.has_weight = False
		
		self.eval_reports = {
			0: self.mass_prop_rep,
			1: self.STEQ_rep,
			2: self.axial_stress_rep,
			3: self.VM_rep,
			4: self.sig_tau_rep,
			5: self.buckling_rep
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
				<th>"""+str(self.th)+"""</th>
			</place>
			<!-- Sup Type= 0:Fixed, 1:Pin, 2: Slot, 3: Thrust -->"""
		for sup in self.supports:
			data += """
			<sup type="{}">
				<axd>{}</axd>
				<th>{}</th>
			</sup>""".format(sup.stype, sup.ax_dist, sup.th)
		data += """
			<!-- Joint Types are the same as Sup Types -->"""
		#I think m1 can be discovered instead of stored.
		#NOTE: joints must be on the axis. That is an assumption of the statics calculations.
		for jt in self.joints:
			data += """
			<jt type="{}">
				<axd>{}</axd>
				<th>{}</th>
			</jt>""".format(jt.stype, jt.axd(self), jt.th)
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
	
	#Return a string representation of the given support or joint,
	#as it relates to this member.
	def sj_str(self, sj):
		s = type(sj).__name__
		s += " joint" if isinstance(sj, jt.Joint) else " support"
		s += " ("+sj.tag+")"
		s += " @axd="+m_str(sj.axd(self))
		return s
	
	#Define the position in the xy plane. Units are meters.
	#th is in degrees
	def place(self, x, y, th):
		self.th = th
		th *= math.pi / 180
		self.x0 = x
		self.y0 = y
		v_ax = self.length*np.array((math.cos(th), math.sin(th)))
		(self.x1, self.y1) = np.array((x,y)) + v_ax
		self.placed = True
	
	def get_pos(self):
		if self.placed:
			return (self.x0, self.y0, self.x1, self.y1)
		return None
	
	#Return an np.array of the coordinates of s0 of the member.
	def get_s0(self):
		return np.array((self.x0, self.y0))
	
	def half_h(self):
		y1_dom = self.xsection.y1_domain()
		assert y1_dom[0] == -y1_dom[1]
		return abs(y1_dom[0])
	
	def v_axis(self):
		assert self.placed
		return np.array((self.x1-self.x0, self.y1-self.y0))
	
	#Unit vector along the axis
	def uv_axis(self):
		#np.array(v_axis) / np.linalg.norm(v_axis)
		return self.v_axis() / self.length
	
	#Find the intersection point of two members, if it exists.
	#If axd0 is already known, it can be supplied as axd.
	#Returns: ((xc, yc), axd0, axd1) or None
	def intersection(self, mem, axd=None):
		uv0 = np.array([self.uv_axis()]).transpose()
		uv1 = np.array([mem.uv_axis()]).transpose()
		A = np.zeros((2,2))
		A[0:2, 0:1] = uv0
		A[0:2, 1:2] = -uv1
		s00 = np.array([[self.x0], [self.y0]])
		s01 = np.array([[mem.x0], [mem.y0]])
		B = np.zeros((2,1))
		B[0:2, 0:1] = s01-s00
		#print(206, A, B)
		try: N = np.linalg.solve(A, B)
		except np.linalg.LinAlgError: return None #Parallel
		#print(209, N)
		if N[0,0] < 0 or N[0,0] > self.length: return None
		if N[1,0] < 0 or N[1,0] > mem.length: return None
		intsx = s00 + uv0*N[0,0]
		if axd is not None and not eps_eq(N[0,0], axd):
			return None
		return ((intsx[0,0], intsx[1,0]), N[0,0], N[1,0])
	
	#Checks if there is a Load, Support, or Joint attached to this member
	#at or close to the given axial distance.
	#Returns: ((xc, yc), axd, lsj) or None
	def lsj_at(self, axd, snap_dist):
		assert self.placed
		for p in self.loads:
			if isinstance(p, Distr_Load):
				(axd0, axd1) = p.limits()
				if axd > axd0-snap_dist and axd < axd1+snap_dist:
					if axd < axd0: axd = axd0
					if axd > axd1: axd = axd1
					(xc, yc) = (self.x0, self.y0) + self.uv_axis() * axd
					return ((xc, yc), axd, p)
			else:
				p_dist = p.ax_dist - axd
				if abs(p_dist) < snap_dist:
					(xc, yc) = (self.x0, self.y0) + self.uv_axis() * p.ax_dist
					return ((xc, yc), p.ax_dist, p)
		for s in self.supports:
			s_dist = s.ax_dist - axd
			if abs(s_dist) < snap_dist:
				(xc, yc) = (self.x0, self.y0) + self.uv_axis() * s.ax_dist
				return ((xc, yc), s.ax_dist, s)
		for j in self.joints:
			j_axd = j.axd(self)
			j_dist = j_axd - axd
			if abs(j_dist) < snap_dist:
				(xc, yc) = (self.x0, self.y0) + self.uv_axis() * j.ax_dist
				return ((xc, yc), j.ax_dist, j)
		return None
	
	#Return a tuple representing the constraints on the member in each direction.
	#It is in this order: (axial, perpindicular, th)
	#I'm not sure if this will work correctly with the new angles.
	#Correction: It does not work with joints!
	def constraints(self):
		c = np.array((0,0,0))
		for s in self.supports:
			sc = s.constraints()
			if eps_round(math.sqrt(sc[0]**2+sc[1]**2)) == 1:
				del_th = s.th - self.th
				#The eps_round does .999999999999999 --> 1
				c[0] += abs(eps_round(math.cos(math.radians(del_th))))
				c[1] += abs(eps_round(math.sin(math.radians(del_th))))
				print("WARNING: member.py 257")
			elif sc[0] == 1 and sc[1] == 1:
				c[0] += 1
				c[1] += 1
			c[2] += sc[2]
		print(269, c)
		return tuple(c)
	
	#Returns the degrees of freedom of the member (axial, perpindicular, th)
	def d_of_f(self):
		c = list(self.constraints())
		#2 perp. sup   ~~=   1 perp. and 1 moment
		if c[1] > 1:
			c[2] += c[1] - 1
			c[1] = 1
		#2 Axial do not create a moment, since they're on the same axis!
		d = np.array((1,1,1)) - c
		return d
	
	#Returns True if the beam is statically determinate or an error
	def is_stdet(self):
		DOF = self.d_of_f()
		if DOF.max() > 0:
			return self.rep_err[0]# "Underconstrained"
		if DOF.min() < 0:
			#Overconstrained -- not statically determinate
			#I can add a case for this later, with delta=PL/EA
			return self.rep_err[1]# "Overconstrained"
		assert np.all(self.d_of_f() == 0)
		return True
	
	#Returns True if the member can support an axial load
	# 1. Must be constrained in translation and rotation
	# 2. Must have precisely one constraint in axial direction
	def sup_axp(self):
		CON = self.constraints()
		#More compact, less descriptive:
		#if CON[0] == 1 and CON[1] > 0: #Ax = 1, Prp constraint
		#	if CON[2] > 0 or CON[1] == 2: #Moment constraint
		#		return True
		if CON[2] > 0 or CON[1] == 2:
			if CON[0] > 0 and CON[1] > 0:
				if CON[0] == 1:
					return True
				else:
					return self.rep_err[1] #Overconstrained in the axial direction
			else:
				return self.rep_err[0] #Underconstrained in Translation
		else:
			return self.rep_err[0] #Underconstrained in Moment
	
	#Return a distributed load representing the weight of the member
	def weight_dl(self):
		q = grav*self.material["rho"]*self.xarea
		wdl = Distr_Load(None, 0, -q, 0, 0, -q, self.length, self.th)
		return wdl
	
	#Return a copy of acting loads. Include weight if we're counting weight.
	def my_loads(self):
		loads_c = self.loads.copy()
		if self.has_weight:
			loads_c.append(self.weight_dl())
		return loads_c
	
	#Return my_loads() after converting Distr_Loads to equivalent point loads.
	def my_loads_pt(self):
		loads = self.my_loads()
		ptlds = []
		for ld in loads:
			if isinstance(ld, Distr_Load):
				ptlds += ld.pt_equiv()
			else:
				ptlds.append(ld)
		return ptlds
	
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
	
	#Takes my_loads() and sums on the reactions at supports and joints
	def my_P_and_R(self):
		#print(337)
		P = self.my_loads()
		reactions = self.static_eq()
		if isinstance(reactions, str):
			return reactions
		for (R, sj) in reactions:
			i = 0
			xym = sj.constraints()
			if xym[0] == 1:
				P.append(Load(None, R[i], 0, sj.axd(self)))
				i += 1
			if xym[1] == 1:
				P.append(Load(None, 0, R[i], sj.axd(self)))
				i += 1
			if (0 < abs(xym[0]) < 1) and (0 < abs(xym[1]) < 1):
				#Support at an angle
				assert eps_round(math.sqrt(xym[0]**2+xym[1]**2)) == 1
				#Fx=R*cos(th); Fy=R*sin(th) -- One unknown
				P.append(Load(None, R[i]*xym[0], R[i]*xym[1], sj.axd(self)))
				i += 1
			if xym[2]:
				P.append(Moment(None, mz=R[i], ax_dist=sj.axd(self)))
		return P
	
	#Return shear as a function of 'd'.
	#Convention: shear that would cause clockwise rotation is positive.
	def shear_symf(self):
		(qx, qy) = self.sum_my_dl()
		PR = self.my_P_and_R()
		if isinstance(PR, str):
			return PR
		V = sym.Integer(0)
		d = sym.symbols('d')
		uvax = self.uv_axis()
		uvprp = np.array((-uvax[1], uvax[0]))
		#Dot product spelled out here (b/c it's all sym)
		V += uvprp[0]*sym.integrate(qx, d) + uvprp[1]*sym.integrate(qy, d)
		for p in PR:
			#The distributed loads were all counted already
			#Moments have no influence on V
			if not isinstance(p, Distr_Load) and not isinstance(p, Moment):
				Vp = np.dot(uvprp, p.get_comp())
				V += Vp*sym.Heaviside(d-p.ax_dist, 1) #, .5)
		#Rewriting as piecewise helps it to integrate the Heaviside
		return V.rewrite(sym.Piecewise).doit()
	
	#Return moment as a function of 'd'.
	#Convention: moment that causes concave up curvature is positive.
	def moment_symf(self):
		PR = self.my_P_and_R()
		if isinstance(PR, str):
			return PR
		V = self.shear_symf()
		if isinstance(V, str):
			return V
		
		d = sym.symbols('d')
		M = sym.integrate(V, d)
		#Zero out s0 - constant of integration
		M -= M.subs(d, 0)
		
		for p in PR:
			if isinstance(p, Moment):
				M -= p.mz*sym.Heaviside(d-p.ax_dist, 1)
		
		return M.rewrite(sym.Piecewise).doit()
	
	#Returns a matrix representing the influence of the given joint or support on the 
	#static equilibrium of the member.
	#	sj: support or joint
	def steq_mat(self, sj):
		mat = np.empty((3,0))
		uv_ax = self.uv_axis()
		sj_axd = sj.ax_dist
		if isinstance(sj, jt.Joint):
			sj_axd = sj.axd(self)
		sj_con = sj.constraints() #Format: (x:0/1,y:0/1,th:0/1)
		#print(359, sj_con, self, sj)
		#To Do: See if these next two are just special cases of the third
		if sj_con[0] == 1:
			v0 = np.array([[1.], [0.], [0.]])
			v0[2, 0] = -sj_axd*uv_ax[1] #F_x influences M
			mat = np.concatenate((mat, v0), axis=1)
		if sj_con[1] == 1:
			v1 = np.array([[0.], [1.], [0.]])
			v1[2, 0] = sj_axd*uv_ax[0] #F_y influences M
			mat = np.concatenate((mat, v1), axis=1)
		if (0 < abs(sj_con[0]) < 1) and (0 < abs(sj_con[1]) < 1):
			#Support at an angle
			assert eps_round(math.sqrt(sj_con[0]**2+sj_con[1]**2)) == 1
			#Fx=R*cos(th); Fy=R*sin(th) -- One unknown
			v = np.array([[sj_con[0]], [sj_con[1]], [0.]])
			v[2, 0] = -sj_axd*uv_ax[1]*sj_con[0] + sj_axd*uv_ax[0]*sj_con[1]
			mat = np.concatenate((mat, v), axis=1)
		if sj_con[2]:
			v2 = np.array([[0.], [0.], [1.]])
			mat = np.concatenate((mat, v2), axis=1)
		return mat
	
	#Returns the matrices for the static equilibrium equations related to all supports 
	#and joints attached to the member and all other connected members.
	#	neg_joints: joints to negate (equal and opposite), since we already
	#		have the positive version.
	#	done_joints: joints that have already been done, positive and negative.
	#	returns: [([SE], m, sj), ([SE], m, sj), ...]
	def steq_mats(self, neg_joints=None, done_joints=None):
		if neg_joints is None:
			neg_joints = []
		if done_joints is None:
			done_joints = []
		mats = []
		#To Do: if this works, remove the list of [js]
		connected = [] #Format: [[m1, [j1,j2]], [m2, j3], ...]
		for s in self.supports:
			mats.append((self.steq_mat(s), self, s))
		for j in self.joints:
			jmat = self.steq_mat(j)
			if j in neg_joints:
				mats.append((-jmat, self, j))
				done_joints.append(j)
				neg_joints.remove(j)
			elif j not in done_joints:
				#First time we're seeing this joint
				mats.append((jmat, self, j))
				neg_joints.append(j)
				other_m = j.other_mem(self)
				#If other_m is already in the list of connected members, then
				#append j to js. Otherwise, add other_m and j.
				for (m, js) in connected:
					if other_m is m:
						js.append(j)
						break
				else:
					connected.append([other_m, [j]])
		for (m, js) in connected:
			#Recursive call
			mats.extend(m.steq_mats(neg_joints, done_joints))
		return mats
	
	#Run the calculations for finding reaction forces with static equilibrium
	#The system must be statically determinate
	#Returns: [(reactions, sj), (reactions, sj), ...]
	def static_eq(self):
		SE_mats = self.steq_mats()
		#print(469, len(SE_mats))
		mems = []
		sjs = []
		sj_cols = []
		#i = 0
		for SEm in SE_mats:
			#print(i); i+=1;
			m = SEm[1]
			#print(473, len(mems), m)
			sj = SEm[2]
			if m not in mems:
				mems.append(m)
			if sj not in sjs:
				sjs.append(sj)
				sj_cols.append(SEm[0].shape[1])
		eqs = 3*len(mems) #sum_x, sum_y, sum_th
		unknowns = sum(sj_cols)
		if eqs > unknowns:
			print(511, eqs, unknowns, SE_mats)
			return self.rep_err[0]
		if eqs < unknowns:
			print(514, eqs, unknowns, SE_mats)
			return self.rep_err[1]
		#print(518, eqs, unknowns)
		STEQ = np.zeros((eqs, unknowns))
		for SEm in SE_mats:
			SE = SEm[0]
			m = SEm[1]
			sj = SEm[2]
			row = mems.index(m)*3
			sji = sjs.index(sj)
			col = sum(sj_cols[:sji])
			rows, cols = SE.shape
			STEQ[row:row+rows, col:col+cols] = SE
		#print(440, SE_mats)
		#print(401, sj_cols)
		#print(442, STEQ)
		
		#The B in Ax = B. sum(R) = -sum(P)
		M_B = np.zeros((eqs, 1))
		for m in mems:
			row = mems.index(m)*3
			s_px = 0
			s_py = 0
			s_m = 0
			uv_ax = m.uv_axis()
			for p in m.my_loads_pt():
				px,py = p.get_comp()
				s_px += px
				s_py += py
				s_m += py*p.ax_dist*uv_ax[0] - px*p.ax_dist*uv_ax[1]
			M_B[row:row+3, 0:1] = np.array([[-s_px], [-s_py], [-s_m]])
		#print(418, M_B)
		try:
			SOL = np.linalg.solve(STEQ, M_B)
		except np.linalg.LinAlgError:
			print(552, STEQ, M_B)
			return self.rep_err[3]
		
		#print(458, sigfig(SOL))
		result = []
		rows_scanned = 0
		for i, sj in enumerate(sjs):
			R = SOL[rows_scanned:rows_scanned+sj_cols[i], 0]
			rows_scanned += sj_cols[i]
			if sj in self.supports or sj in self.joints:
				result.append((R, sj))
		return result
	
	#Do mohr's circle transformations on an element with sigma_x and tau_xy (assumes no sig_y)
	#Return two vectors, each in polar and cartesian: smvp, smvc, tmvp, tmvc
	#fm = module with appropriate functions. Default numpy.
	#not quite working yet for symbolic expressions
	@staticmethod
	def mohr_trsfm(sigx, tau, fm=np):
		#Force the angle to be more horizontal so you get negative instead of perpindicular
		convert_to_neg = False
		#numpy uses "arctan" but sympy and math use "atan"
		try: atan = fm.arctan
		except AttributeError: atan = fm.atan
		
		#print(294, tau/sigx, type(tau/sigx))
		#Does not work with sym expr
		th_a = np.array(atan(2*tau/sigx)/2, dtype="float64")
		#print(295, th_a, type(th_a))
		#Replace oo with pi/4 (which is (pi/2)/2)
		th_a = np.nan_to_num(th_a, nan=np.pi/4, posinf=np.pi/4, neginf=-np.pi/4)
		#print(299, th_a, type(th_a))
		th_b = th_a + math.pi/2
		sig_a = sigx/2*(1+fm.cos(2*th_a)) + tau*fm.sin(2*th_a)
		sig_b = sigx/2*(1+fm.cos(2*th_b)) + tau*fm.sin(2*th_b)
		#print(300, sigx, tau, th_a, th_b, sig_a, sig_b)
		#I converted the ternaries to array-compatible logic
		#Must be int because sympy isn't dealing with T/F multiplication
		b_gt_a = np.array(sig_b > sig_a, dtype="bool")
		b_agt_a = np.array(abs(sig_b) > abs(sig_a), dtype="bool")
		#th_smax = th_b if abs(sig_b) > abs(sig_a) else th_a
		th_smax = th_b*b_agt_a + th_a*(1 - b_agt_a)
		#sig_max = sig_b if abs(sig_b) > abs(sig_a) else sig_a
		sig_max = sig_b*b_agt_a + sig_a*(1 - b_agt_a)
		#th_s1 = th_b if sig_b > sig_a else th_a #alg. max.
		th_s1 = th_b*b_gt_a + th_a*(1 - b_gt_a)
		th_t1 = th_s1 - math.pi/4 #Alg. max tau
		th_t2 = th_t1 + math.pi/2
		if convert_to_neg:
			t1_horiz = math.pi/2 - abs(th_t1%math.pi - math.pi/2)
			t2_horiz = math.pi/2 - abs(th_t2%math.pi - math.pi/2)
			#is t1 the more horizontal one?
			t1_is_h = t2_horiz >= t1_horiz
			th_tmax = th_t1*t1_is_h + th_t2*(1-t1_is_h)
			tau_max = -.5*sigx*fm.sin(2*th_tmax) + tau*fm.cos(2*th_tmax)
		else: #Just accept the (positive) alg. max. tau
			th_tmax = th_t1
			tau_max = fm.sqrt((sigx/2)**2 + tau**2) #alg. max
		#return (sig_max, th_smax), (tau_max, th_tmax)
		sig_max_vp = (sig_max, th_smax)
		sig_max_vc = (sig_max*fm.cos(th_smax), sig_max*fm.sin(th_smax))
		tau_max_vp = (tau_max, th_tmax)
		tau_max_vc = (tau_max*fm.cos(th_tmax), tau_max*fm.sin(th_tmax))
		#print(331, sig_max_vp, sig_max_vc, tau_max_vp, tau_max_vc)
		return sig_max_vp, sig_max_vc, tau_max_vp, tau_max_vc
	
	#Return report text (& any figures) for each type of evaluation
	#Be sure to make this line up with the dictionary in Lab
	def gen_report(self, type):
		return self.eval_reports[type]()
	
	#Static Equilibrium Report
	#Lists the reaction forces at each support and joint along the member
	def STEQ_rep(self):
		STEQ = self.static_eq()
		if isinstance(STEQ, str):
			return STEQ
		rep_text = "Reactions at supports and joints: "
		for (R, sj) in STEQ:
			#The order is important here, since (Joint << Support)
			rep_text += "\n" + self.sj_str(sj) + ": "
			xym = sj.constraints()
			i = 0
			if xym[0] == 1:
				rep_text += "\n\tR_x=" + N_to_kN_str(R[i])
				i += 1
			if xym[1] == 1:
				rep_text += "\n\tR_y=" + N_to_kN_str(R[i])
				i += 1
			if (0 < abs(xym[0]) < 1) and (0 < abs(xym[1]) < 1):
				#Support at an angle
				assert eps_round(math.sqrt(xym[0]**2+xym[1]**2)) == 1
				rep_text += "\n\tR=" + N_to_kN_str(R[i])
				th = math.degrees(math.atan2(xym[1], xym[0]))
				rep_text += " @"+str(sigfig(th))+"\u00B0"
				i += 1
			if xym[2]:
				rep_text += "\n\tR_M=" + Nm_to_kNm_str(R[i])
		return rep_text
	
	#Mass Properties Report
	def mass_prop_rep(self):
		rep_text = "cross sectional area A = "+str(sigfig(self.xarea))+" m^2"
		rep_text += "\n    | second moment of area about horizontal axis I_x = "+str(sigfig(self.xIx))+" m^4"
		rep_text += "\n    | second moment of area about vertical axis I_y = "+str(sigfig(self.xIy))+" m^4"
		rep_text += "\nvolume V = "+str(sigfig(self.volume))+" m^3"
		rep_text += "\nmass m = "+str(sigfig(self.mass))+" kg"
		rep_text += "\nweight W = "+str(sigfig(self.weight))+" N"
		rep_text += "\n\nProperties of "+self.material["name"]+":"
		rep_text += "\n    | mass density \u03C1 = "+str(self.material["rho"])+" kg/m^3"
		rep_text += "\n    | Young's modulus E = "+str(self.material["E"]/1e9)+" GPa"
		rep_text += "\n    | yield stress \u03C3_y = "+str(self.material["sig_y"]/1e6)+" MPa"
		rep_text += "\n    | ultimate stress \u03C3_y = "+str(self.material["sig_u"]/1e6)+" MPa"
		return rep_text
	
	#Return internal axial loads as a function of 'd'.
	#Convention: tension is positive.
	def axial_loads_sym(self):
		d = sym.symbols('d')
		uv_ax = self.uv_axis()
		(qx, qy) = self.sum_my_dl()
		P_ax = -( sym.integrate(qx*uv_ax[0], d) + sym.integrate(qy*uv_ax[1], d) )
		PR = self.my_P_and_R()
		if isinstance(PR, str):
			return PR
		for p in PR:
			#The distributed loads were all counted already
			#Moments have no influence on P_ax
			if not isinstance(p, Distr_Load) and not isinstance(p, Moment):
				p_ax = np.dot(uv_ax, p.get_comp())
				P_ax -= p_ax*sym.Heaviside(d-p.ax_dist, 1) #, .5)
		return P_ax.rewrite(sym.Piecewise).doit()
	
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
	
	#Report of only the stress and strain due to axial loads
	def axial_stress_rep(self):
		P = self.axial_loads_sym()
		if isinstance(P, str):
			return P
		d = sym.symbols("d")
		(Pmax, dPmax), (Pmin, dPmin) = max1d(P, d, (0, self.length))
		
		rep_text = "'d' is measured (in m) from end \"zero\" of the member."
		if Pmax > 0:
			rep_text += "\nMax tension is at "+m_str(dPmax)
			rep_text += "\nMax tension = "+N_to_kN_str(Pmax)
			Smax = Pmax/self.xarea
			rep_text += "\nMax tensile stress \u03C3 = "+Pa_to_MPa_str(Smax)
			typer = round(Smax / self.material["sig_y"] * 100, 3)
			tuper = round(Smax / self.material["sig_u"] * 100, 3)
			rep_text += "\nWhich is "+str(typer)+"% of the yield stress and "
			rep_text += str(tuper)+"% of the ultimate stress."
			rep_text += "\nMax tensile strain \u03B5 = "+str(sigfig(Smax/self.E))
		if Pmin < 0:
			rep_text += "\nMax compression is at "+m_str(dPmin)
			rep_text += "\nMax compression = "+N_to_kN_str(-Pmin)
			Smin = Pmin/self.xarea
			rep_text += "\nMax compressive stress \u03C3 = "+Pa_to_MPa_str(-Smin)
			cyper = round(-Smin / self.material["sig_y"] * 100, 3)
			cuper = round(-Smin / self.material["sig_u"] * 100, 3)
			rep_text += "\nWhich is "+str(cyper)+"% of the yield stress and "
			rep_text += str(cuper)+"% of the ultimate stress."
			rep_text += "\nMax compressive strain \u03B5 = "+str(sigfig(-Smin/self.E))
			rep_text += "\nThis evaluation does not consider buckling."
		
		#Make plot for P(d)
		Pf = sym.lambdify(d, P / 1e3)
		if P.is_number:
			#print(738, P)
			#print(739, float(P))
			Pf = lambda d: 0*d + float(P)
		
		dax = np.linspace(0, self.length, self.d_resolution)
		Pax = Pf(dax)
		
		fig, sp = plt.subplots(figsize=(7,3))
		sp.plot(dax, Pax, color="red")
		sp.grid(True)
		sp.set_title("Axial Load P (kN)")
		sp.set(xlabel="Axial Distance d (m)")
		plt.subplots_adjust(*self.plot_margins, wspace=.09)
		return (rep_text, fig)
	
	#Euler Buckling Report
	#To Do: Fix it so it deals with st. overconstrained columns like Fixed-Slot(ax)
	def buckling_rep(self):
		STEQ = self.static_eq()
		if isinstance(STEQ, str):
			return STEQ
		P = self.axial_loads_sym()
		if isinstance(P, str):
			return P
		Pcr = math.pi**2 * self.E * self.xImin / self.length**2
		d = sym.symbols("d")
		
		STEQ.sort(key=lambda r_sj: r_sj[1].axd(self))
		#Make a list of intervals along the member between each support or joint
		intervals = []
		last_axd = 0
		last_sj = None
		for r_sj in STEQ:
			sj = r_sj[1]
			sj_axd = sj.axd(self)
			if last_axd < sj_axd:
				intervals.append(((last_axd, last_sj), (sj_axd, sj)))
			last_axd = sj_axd
			last_sj = sj
		if last_axd < self.length:
			intervals.append(((last_axd, last_sj), (self.length, None)))
		rep_text = "Measuring distance 'x' (in m) from end 'zero' of the member,"
		for sdint in intervals:
			(axd0, sj0), (axd1, sj1) = sdint
			rep_text += "\nfor x \u03F5 "+coords_str(axd0, axd1)
			P_sub = f_restrict(P, d, (axd0, axd1), lims=False)
			#print(782, P_sub)
			#Check that P_sub is constant
			#(2 steps since sym.is_constant() is flawed with Piecewise)
			if len(pw_sdls(P_sub, d)) <= 2 and P_sub.is_constant():
				d_test = (axd0 + axd1) / 2 #avg of the limits --> in the middle
				Pc = -float(P_sub.subs(d, d_test))
				if Pc <= 0:
					rep_text += ", there is no compressive load."
				else:
					rep_text += ", the compressive load Pc="+N_to_kN_str(Pc)
					c0 = (0,0,0) if sj0 is None else sj0.constraints()
					c1 = (0,0,0) if sj1 is None else sj1.constraints()
					#print(798, c0, c1)
					if isinstance(sj0, jt.Joint) or isinstance(sj1, jt.Joint):
						rep_text += "\n\tWarning: this Euler Buckling analysis may be inaccurate for"
						rep_text += " complex structures, since it was designed for simple columns."
					#How many perpindicular constraints? How many moment constraints?
					c = [0,0,0]
					if eps_round(math.sqrt(c0[0]**2+c0[1]**2)) == 1:
						del_th = sj0.th - self.th
						print(811, del_th, sj0.th, self.th, c)
						c[0] += abs(eps_round(math.cos(math.radians(del_th))))
						c[1] += abs(eps_round(math.sin(math.radians(del_th))))
					elif c0[0] == 1 and c0[1] == 1:
						c[0] += 1
						c[1] += 1
					c[2] += c0[2]
					if eps_round(math.sqrt(c1[0]**2+c1[1]**2)) == 1:
						del_th = sj1.th - self.th
						print(820, del_th, sj1.th, self.th, c)
						c[0] += abs(eps_round(math.cos(math.radians(del_th))))
						c[1] += abs(eps_round(math.sin(math.radians(del_th))))
					elif c1[0] == 1 and c1[1] == 1:
						c[0] += 1
						c[1] += 1
					c[2] += c1[2]
					assert c[0] == 1, "Something confusing is happening. c=({}, {}, {})".format(
						c[0], c[1], c[2])
					Pcr_i = Pcr
					if c[1] == 1 and c[2] == 1: #Fixed-Free or Pin-M or Sy-Tx or Ty-Sx
						Pcr_i /= 4
					elif c[1] == 2 and c[2] == 1: #Fixed-Sx or Pin-Tx
						Pcr_i *= 2.046 #From book - table on columns
					elif c[1] == 2 and c[2] == 2: #Fixed-Tx
						Pcr_i *= 4
					elif c[1] == 2 and c[2] == 0: #Pin-Sx
						pass #Pcr_i *= 1
					elif c[1] == 1 and cp[2] == 2: #Fixed-M or Ty-Tx
						pass #Pcr_i *= 1
					else:
						Pcr_i = None
					if Pcr_i is None:
						rep_text += "\n\tI don't understand the supports for this column."
					else:
						Pper = round(Pc / Pcr_i * 100, 3)
						rep_text += "\n\tThis is " + str(Pper) + "% of the critical load Pcr="
						rep_text += N_to_kN_str(Pcr_i)
						if Pc < Pcr_i:
							rep_text += "\n\tSo the member is stable in this interval."
						if Pc >= Pcr_i:
							rep_text += "\n\tSo the member is unstable in this interval."
			else:
				rep_text += ", the axial load is not constant."
				rep_text += " Sorry, I don't know how to deal with buckling when"
				rep_text += " the axial load is not constant."
				#I could give it a shot later, but I think it includes messing with the
				#differential equation and solving it generically instead of using the 4[/10] cases.
		return rep_text
	
	#Report on shear and moment
	#Returns text and a pyplot figure
	def VM_rep(self):
		rep_text = "'d' is measured (in m) from end \"zero\" (W or S) of the member."
		rep_text += "\nFor internal tension and stresses, see \"Axial and Shear Stress\""
		#print(797)
		V = self.shear_symf()
		if V in self.rep_err:
			return V
		#print(801)
		M = self.moment_symf()
		if M in self.rep_err:
			return M
		#Make plots for V and M
		d = sym.symbols("d")
		Vf = sym.lambdify(d, V/1000)
		Mf = sym.lambdify(d, M/1000)
		dax = np.linspace(0, self.length, self.d_resolution)
		Vax = Vf(dax)
		Max = Mf(dax)
		try:
			assert len(dax) == len(Vax)
			assert len(dax) == len(Max)
		except (AssertionError, TypeError) as e:
			print(860, e)
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
		#print(586, ": sig (Pa) =", sig)
		sig = sig/1e6 #MPa
		y1min, y1max = self.xsection.y1_domain()
		tau_res = self.tau_sym()
		if isinstance(tau_res, str):
			return tau_res
		tau, h_dom = tau_res
		tau = tau/1e6
		y1_rng = y1max - y1min
		h_rng = h_dom[1] - h_dom[0]
		sigf = sym.lambdify((d,h), sig, "numpy")
		tauf = sym.lambdify((d,h), tau, "numpy")
		#Fix this for zero / constant functions. See VM_rep assertions
		#This way it returns an array of the same size if arrays are passed
		if sig.is_number:
			sigf = lambda d,h: 0*d + float(sig)
		if tau.is_number:
			tauf = lambda d,h: 0*d + float(tau)
		
		#TEMP
		#print(608, "y1 domain:", [y1min, y1max])
		#print(609, "h_dom:", h_dom)
		#(smvp, smvc, tmvp, tmvc) = self.mohr_trsfm(sig, tau, fm=sym)
		#smx, smy = smvc
		#/TEMP
		d_ls = np.linspace(0, self.length, self.d_resolution)
		h_sig_ls = np.linspace(y1min, y1max, self.h_resolution)
		h_tau_ls = np.linspace(*h_dom, self.h_resolution)
		D, H_s = np.meshgrid(d_ls, h_sig_ls)
		_, H_t = np.meshgrid(d_ls, h_tau_ls)
		
		#Sigma - report 1
		SIG = sigf(D, H_s)
		#print(606, ": sig =", sig)
		(sig_max, smax_d, smax_h), (sig_min, smin_d, smin_h) = max2d(sig, d, h, 
			(0, self.length), (y1min, y1max), xres=self.d_resolution*5, yres=self.h_resolution*5)
		print(609, ": smax,min", sig_max, sig_min)
		sig_rng = max(abs(sig_min), abs(sig_max))
		s1_rep_text = "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		s1_rep_text += " and 'h' (in mm) from the neutral axis of bending,"
		if sig_max > 0:
			s1_rep_text += "\nMax Tensile Axial Stress = " + str(sigfig(sig_max)) + " MPa"
			s1_rep_text += " at d="+str(sigfig(smax_d))+"m, h="+str(sigfig(smax_h*1e3))+"mm"
		if sig_min < 0:
			s1_rep_text += "\nMax Compressive Axial Stress = " + str(sigfig(-sig_min)) + " MPa"
			s1_rep_text += " at d="+str(sigfig(smin_d))+"m, h="+str(sigfig(smin_h*1e3))+"mm"
		#Plot of in-plane axial stress
		fig1, ax1 = plt.subplots(figsize=(7,3))
		ax1.set_title("Axial Stress \u03C3 (MPa)")
		ax1.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		#print(651, SIG)
		#print(652, [0, self.length, y1min*1e3, y1max*1e3])
		im1 = ax1.imshow(SIG, cmap=plt.cm.RdBu, interpolation="bilinear", aspect="auto", 
			extent=[0, self.length, y1min*1e3, y1max*1e3], vmin=-sig_rng, vmax=sig_rng, origin="lower")
		fig1.colorbar(im1, ax=ax1)
		plt.subplots_adjust(*self.implot_margins, wspace=.09)
		reps.append(("In-plane \u03C3", s1_rep_text, fig1))
		
		#TAU - report 2
		#Plot of in-plane shear stress TAU
		fig2, ax2 = plt.subplots(figsize=(7,3))
		ax2.set_title("Shear Stress \u03C4 (MPa)")
		t1_rep_text = ""
		if h_rng == 0:
			t1_rep_text += self.rep_wrng[0]
			(tau_max, tmax_d), (tau_min, tmin_d) = max1d(tau.subs(h, h_dom[0]),
				d, (0, self.length))
			tmin_h = 0
			tmax_h = 0
			ax2.set(xlabel="Axial Distance d (m)", ylabel="Shear at Neutral Axis (MPa)")
			TAU_nax = tauf(d_ls, h_dom[0])#TAU along the neutral axis
			ax2.plot(d_ls, TAU_nax)
			plt.subplots_adjust(*self.plot_margins, wspace=.09)
		else:
			if h_rng < y1_rng:
				t1_rep_text += self.rep_wrng[1]
			(tau_max, tmax_d, tmax_h), (tau_min, tmin_d, tmin_h) = max2d(tau, d, h, 
				(0, self.length), h_dom, xres=self.d_resolution*5, yres=self.h_resolution*5)
			tau_rng = max(abs(tau_min), abs(tau_max))
			ax2.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
			TAU = tauf(D, H_t)
			im2 = ax2.imshow(TAU, cmap=plt.cm.BrBG, interpolation="bilinear", aspect="auto", origin="lower",
				extent=[0, self.length, h_dom[0]*1e3, h_dom[1]*1e3], vmin=-tau_rng, vmax=tau_rng)
			fig2.colorbar(im2, ax=ax2)
			plt.subplots_adjust(*self.implot_margins, wspace=.09)
		t1_rep_text += "Measuring 'd' (in m) from end zero (left or bottom) of the member"
		t1_rep_text += " and 'h' (in mm) from the neutral axis of bending,"
		if tau_max > 0:
			t1_rep_text += "\nMax Positive Sheer Stress = " + str(sigfig(tau_max)) + " MPa"
			t1_rep_text += " at d="+str(sigfig(tmax_d))+"m, h="+str(sigfig(tmax_h*1e3))+"mm"
		if tau_min < 0:
			t1_rep_text += "\nMax Negative Sheer Stress = " + str(sigfig(abs(tau_min))) + " MPa"
			t1_rep_text += " at d="+str(sigfig(tmin_d))+"m, h="+str(sigfig(tmin_h*1e3))+"mm"
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
		
		s2_rep_text = ""
		t2_rep_text = ""
		if h_rng == 0:
			s2_rep_text += self.rep_wrng[0]
			t2_rep_text += self.rep_wrng[0]
			#Use a different kind of plot?
			#Though I'm not sure what I'd like it to look like.
		elif h_rng < y1_rng:
			s2_rep_text += self.rep_wrng[1]
			t2_rep_text += self.rep_wrng[1]
		s2_rep_text += "Measuring 'd' (in m) from end zero (left or bottom) of the member"
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
		t2_rep_text += "Measuring 'd' (in m) from end zero (left or bottom) of the member"
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
		plt.subplots_adjust(*self.implot_margins, wspace=.09)
		reps.append(("Out-of-plane \u03C3", s2_rep_text, fig3))
		#Plot of principal (out-of-plane) shear stress
		fig4, ax4 = plt.subplots(figsize=(7,3))
		ax4.set_title("Principal Shear Stress \u03C4_p (MPa)")
		ax4.set(xlabel="Axial Distance d (m)", ylabel="Height h (mm)")
		q2 = ax4.quiver(Dq, Hq*1e3, T1x, T1y, T1m, cmap=plt.cm.BrBG, clim=(-tau_rng, tau_rng),
			pivot="mid", headaxislength=0, headlength = 0, headwidth=1, width=.009)
		fig4.colorbar(q2, ax=ax4)
		plt.subplots_adjust(*self.implot_margins, wspace=.09)
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



def steq_test():
	import support as sup
	matl = None
	xsec = None
	m0 = Member(matl, xsec, 1.0)
	m0.place(1, 1, 0)
	m1 = Member(matl, xsec, 1.0)
	m1.place(1.25,1.25, -45)
	s0 = sup.Fixed(None, 0, 0)
	m0.supports.append(s0)
	
#steq_test()
