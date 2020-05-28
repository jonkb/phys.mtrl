import math
import numpy as np
import sympy as sym

class Load:
	ah_width = 6#px
	ball_rad = ah_width*.65
	load_types = {
		0: "Point",
		1: "Distributed"
	}
	def __init__(self, tag, xc, yc, ax_dist=0):
		self.tag = tag
		self.xc = xc
		self.yc = yc
		self.ax_dist = ax_dist
	def get_comp(self):
		return (self.xc,self.yc)
	#Draw a point load. All args in pixels. The tip is at the point specified by (px,py)
	def draw(self, lab, px, py):
		if self.xc == 0 and self.yc == 0:
			lab.canv.create_oval(px - self.ball_rad, py - self.ball_rad, px + self.ball_rad, 
				py + self.ball_rad, fill="red", outline="red", tags=self.tag)
			return
		pxc, pyc = lab.kN_to_px(self.xc/1000, self.yc/1000)
		lab.canv.create_line(px-pxc,py-pyc,px,py, width=2, fill="red", tags=self.tag)
		L = math.sqrt(pxc**2+pyc**2)
		ux = pxc/L
		uy = pyc/L
		Trx = px - self.ah_width*3*ux + self.ah_width/2*uy
		Try = py - self.ah_width*3*uy - self.ah_width/2*ux
		Tlx = px - self.ah_width*3*ux - self.ah_width/2*uy
		Tly = py - self.ah_width*3*uy + self.ah_width/2*ux
		T_points = [px,py, Trx,Try, Tlx,Tly]
		lab.canv.create_polygon(T_points, fill="red", outline="red", tags=self.tag)# width=1,

#Quadrelateral distributed loads, defined by 2 q vectors @ axd0 & axd1
class Distr_Load(Load):
	#Dist btw arrows in px
	ad_px = 25
	def __init__(self, tag, xc0, yc0, axd0, xc1, yc1, axd1, isv):
		super().__init__(tag, xc0, yc0)
		self.xc0 = xc0
		self.yc0 = yc0
		self.axd0 = axd0
		self.xc1 = xc1
		self.yc1 = yc1
		self.axd1 = axd1
		self.isv = isv
		self.a_dist = self.ad_px / 360
	def pt_equiv(self):
		L = self.axd1 - self.axd0
		Px = self.xc0*L + (self.xc1-self.xc0)*L/2
		Py = self.yc0*L + (self.yc1-self.yc0)*L/2
		if self.isv:
			if Px == 0:
				axd = L/2 # axd Doesn't matter, since the force is alligned with mem axis
			else:
				axd = ( self.xc0*L**2/2 + (self.xc1-self.xc0)*L**2/3 ) / Px
		else:
			if Py == 0:
				axd = L/2
			else:
				axd = ( self.yc0*L**2/2 + (self.yc1-self.yc0)*L**2/3 ) / Py
		ptl = Load(self.tag, Px, Py, self.axd0 + axd)
		return ptl
	def get_comp(self):
		ptl = self.pt_equiv()
		return ptl.get_comp()
	@property
	def ax_dist(self):
		ptl = self.pt_equiv()
		return ptl.ax_dist
	#Set axd doesn't mean anything here
	@ax_dist.setter
	def ax_dist(self, axd):
		pass
	#Include both ends, but keep the spacing relatively standard
	def ax_rng(self, ax0, ax1):
		L = ax1 - ax0
		n = abs(round(L / self.a_dist))
		if n == 0:
			n = 1
		step = L / n
		return np.arange(ax0, ax1+step/2, step)
	#(px,py) is for the position in px of q0
	def draw(self, lab, px, py):
		self.a_dist = self.ad_px / lab.px_per_m
		Ld = self.axd1 - self.axd0 #(Neg. if backwards)
		for axd in self.ax_rng(self.axd0, self.axd1):
			lx = self.xc0 + (self.xc1 - self.xc0) * (axd-self.axd0)/Ld
			ly = self.yc0 + (self.yc1 - self.yc0) * (axd-self.axd0)/Ld
			pxc, pyc = lab.kN_to_px(lx/1000, ly/1000)
			Ll = math.sqrt(pxc**2+pyc**2)
			if Ll == 0:
				lab.canv.create_oval(px - self.ball_rad, py - self.ball_rad, px + self.ball_rad, 
					py + self.ball_rad, fill="red", outline="red", tags=self.tag)
				continue
			if self.isv:
				ax = px
				ay = py - lab.px_per_m*(axd-self.axd0)
			else:
				ax = px + lab.px_per_m*(axd-self.axd0)
				ay = py
			lab.canv.create_line(ax-pxc,ay-pyc,ax,ay, width=2, fill="red", tags=self.tag)
			ux = pxc/Ll #Unit vector
			uy = pyc/Ll
			ah_l = self.ah_width*3
			Trx = ax - ah_l*ux + self.ah_width/2*uy
			Try = ay - ah_l*uy - self.ah_width/2*ux
			Tlx = ax - ah_l*ux - self.ah_width/2*uy
			Tly = ay - ah_l*uy + self.ah_width/2*ux
			T_points = [ax,ay, Trx,Try, Tlx,Tly]
			lab.canv.create_polygon(T_points, fill="red", outline="red", tags=self.tag)
		q0xc, q0yc = lab.kN_to_px(self.xc0/1000, self.yc0/1000)
		q1xc, q1yc = lab.kN_to_px(self.xc1/1000, self.yc1/1000)
		q0ax = px - q0xc
		q0ay = py - q0yc
		if self.isv:
			q1ax = px - q1xc
			q1ay = py - Ld*lab.px_per_m - q1yc
		else:
			q1ax = px + Ld*lab.px_per_m - q1xc
			q1ay = py - q1yc
		lab.canv.create_line(q0ax, q0ay, q1ax, q1ay, width=2, fill="red", tags=self.tag)
	#Return two sympy functions for the load at a given axial d (Qx(d), Qy(d))
	def to_symf(self):
		d = sym.symbols("d")
		Ld = self.axd1 - self.axd0
		chr_pls = (sym.Heaviside(d-self.axd0, .5)-sym.Heaviside(d-self.axd1, .5))
		px = self.xc0 + (d-self.axd0)/Ld*(self.xc1-self.xc0)
		py = self.yc0 + (d-self.axd0)/Ld*(self.yc1-self.yc0)
		return (px*chr_pls, py*chr_pls)
	
	
	
