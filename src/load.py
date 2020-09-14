import math
import numpy as np
import sympy as sym

class Load:
	ah_width = 6#px
	ball_rad = ah_width*.65
	load_types = ("Point", "Distributed", "Moment")
	#Load Type for this class (0, 1, 2)
	ltype = 0
	def __init__(self, tag, xc, yc, ax_dist=0):
		self.tag = tag
		self.xc = xc
		self.yc = yc
		self.ax_dist = ax_dist
	def to_xml(self):
		data = """
			<ld type="{}">
				<xc>{}</xc>
				<yc>{}</yc>
				<axd>{}</axd>
			</ld>""".format(self.ltype, self.xc, self.yc, self.ax_dist)
		return data
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
		T_points = self.arrow_head(px, py, pxc/L, pyc/L)
		lab.canv.create_polygon(T_points, fill="red", outline="red", tags=self.tag)# width=1,
	
	#Given a point on the canvas and a unit vector, draw an arrow head pointing to (px,py)
	#All units in pixels on the canvas, not points on the grid.
	def arrow_head(self, px, py, ux, uy):
		ah_l = self.ah_width*3
		Trx = px - ah_l*ux + self.ah_width/2*uy
		Try = py - ah_l*uy - self.ah_width/2*ux
		Tlx = px - ah_l*ux - self.ah_width/2*uy
		Tly = py - ah_l*uy + self.ah_width/2*ux
		return [px,py, Trx,Try, Tlx,Tly]
	

#Quadrelateral distributed loads, defined by 2 q vectors @ axd0 & axd1
class Distr_Load(Load):
	#Dist btw arrows in px
	ad_px = 25
	ltype = 1
	#Note: th is the angle (deg) of the axis of the member
	def __init__(self, tag, xc0, yc0, axd0, xc1, yc1, axd1, th):
		super().__init__(tag, xc0, yc0)
		self.xc0 = xc0
		self.yc0 = yc0
		self.axd0 = axd0
		self.xc1 = xc1
		self.yc1 = yc1
		self.axd1 = axd1
		self.th = th
		self.a_dist = self.ad_px / 360
	def to_xml(self):
		data = """
			<ld type="{}">
				<xc0>{}</xc0>
				<yc0>{}</yc0>
				<axd0>{}</axd0>
				<xc1>{}</xc1>
				<yc1>{}</yc1>
				<axd1>{}</axd1>
			</ld>""".format(self.ltype, self.xc0, self.yc0, 
			self.axd0, self.xc1, self.yc1, self.axd1)
		return data
	
	#Converts this distributed load into an equivalent list of point loads.
	#It needs to be multiple point loads because distributed loads may cause moments.
	def pt_equiv(self):
		ptlds = []
		L = self.axd1 - self.axd0
		#TO DO: Check to see if this all works when it's backwards so L<0
		assert L != 0
		if L < 0:
			print("WARNING: load.py 80: Not sure if this works backwards yet.")
		
		#break x-component into loads
		if self.xc0 != 0 or self.xc1 != 0:
			if np.sign(self.xc0) == -np.sign(self.xc1):
				#Must be 2 loads, since the sign flipped. Otherwise the moment will be lost.
				#axdi: axd of inflection point
				axdi = self.axd0 + L*self.xc0/(self.xc0-self.xc1)
				p0 = self.xc0*(axdi-self.axd0)/2
				axdp0 = self.axd0 + (axdi-self.axd0)/3
				ptlds.append(Load(self.tag, p0, 0, axdp0))
				p1 = self.xc1*(self.axd1-axdi)/2
				axdp1 = self.axd1 - (self.axd1-axdi)/3
				ptlds.append(Load(self.tag, p1, 0, axdp1))
			else:
				p = (self.xc0+self.xc1)/2*L
				axdp = self.axd0 + (self.xc0/2*L**2 + (self.xc1-self.xc0)/3*L**2)/p
				ptlds.append(Load(self.tag, p, 0, axdp))
		
		#break y-component into loads
		if self.yc0 != 0 or self.yc1 != 0:
			if np.sign(self.yc0) == -np.sign(self.yc1):
				#print(102, self.axd0, self.axd1, L, self.yc0, self.yc1)
				axdi = self.axd0 + L*self.yc0/(self.yc0-self.yc1)
				p0 = self.yc0*(axdi-self.axd0)/2
				axdp0 = self.axd0 + (axdi-self.axd0)/3
				ptlds.append(Load(self.tag, 0, p0, axdp0))
				p1 = self.yc1*(self.axd1-axdi)/2
				axdp1 = self.axd1 - (self.axd1-axdi)/3
				ptlds.append(Load(self.tag, 0, p1, axdp1))
			else:
				p = (self.yc0+self.yc1)/2*L
				axdp = self.axd0 + (self.yc0/2*L**2 + (self.yc1-self.yc0)/3*L**2)/p
				ptlds.append(Load(self.tag, 0, p, axdp))
				
		return ptlds
	
	def get_comp(self):
		print("WARNING load.py 117")
		return "WARNING: Not sure what you mean"
		ptl = self.pt_equiv()
		return ptl.get_comp()
	@property
	def ax_dist(self):
		print(125, "WARNING, I dont' think this works.")
		ptl = self.pt_equiv()
		return ptl.ax_dist
	#Set axd doesn't mean anything here
	@ax_dist.setter
	def ax_dist(self, axd):
		pass
	#Return axd0 and axd1, in increasing order
	def limits(self):
		return (min(self.axd0, self.axd1), max(self.axd0, self.axd1))
	#Include both ends, but keep the spacing relatively standard
	def ax_rng(self, ax0, ax1):
		L = ax1 - ax0
		n = abs(round(L / self.a_dist))
		if n == 0:
			n = 1
		step = L / n
		return np.arange(ax0, ax1+step/2, step)
	#(px,py) is for the position in pixels of q0
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
			#ax, ay is the location of the arrowhead tip
			ax = px + lab.px_per_m*(axd-self.axd0)*math.cos(math.radians(self.th))
			ay = py - lab.px_per_m*(axd-self.axd0)*math.sin(math.radians(self.th))
			#print(152, ax, ay, axd)
			lab.canv.create_line(ax-pxc,ay-pyc,ax,ay, width=2, fill="red", tags=self.tag)
			ah_l = self.ah_width*3
			ux = pxc/Ll #Unit vector along arrow
			uy = pyc/Ll
			T_points = self.arrow_head(ax, ay, ux, uy)
			lab.canv.create_polygon(T_points, fill="red", outline="red", tags=self.tag)
		q0xc, q0yc = lab.kN_to_px(self.xc0/1000, self.yc0/1000)
		q1xc, q1yc = lab.kN_to_px(self.xc1/1000, self.yc1/1000)
		q0ax = px - q0xc
		q0ay = py - q0yc
		q1ax = px - q1xc + Ld*lab.px_per_m*math.cos(math.radians(self.th))
		q1ay = py - q1yc - Ld*lab.px_per_m*math.sin(math.radians(self.th))
		lab.canv.create_line(q0ax, q0ay, q1ax, q1ay, width=2, fill="red", tags=self.tag)
	#Return two sympy functions for the load at a given axial d (Qx(d), Qy(d))
	def to_symf(self):
		d = sym.symbols("d")
		Ld = self.axd1 - self.axd0
		chr_pls = (sym.Heaviside(d-self.axd0, 1)-sym.Heaviside(d-self.axd1, 1))
		px = self.xc0 + (d-self.axd0)/Ld*(self.xc1-self.xc0)
		py = self.yc0 + (d-self.axd0)/Ld*(self.yc1-self.yc0)
		return (px*chr_pls, py*chr_pls)

#Applied moments
class Moment(Load):
	ltype = 2
	#These next two control (together with lab.px_per_kN) the scaling of Z-Moments on the canvas.
	#18kN --> 1m linear. 18kN-m --> 40px radius (by default)
	large_kNm = 18
	M_minscale = 0.5
	arrow_offset = Load.ah_width*2.5#*3 is the arrow length
	def __init__(self, tag, mx=0, my=0, mz=0, ax_dist=0):
		#NOTE: this used to be (tag, 0, 0, ax_dist)
		#Hopefully I didn't break anything
		super().__init__(tag, mx, my, ax_dist)
		self.mx = mx
		self.my = my
		self.mz = mz
	def to_xml(self):
		data = """
			<ld type="{}">
				<xc>{}</xc>
				<yc>{}</yc>
				<zc>{}</zc>
				<axd>{}</axd>
			</ld>""".format(self.ltype, self.mx, self.my, self.mz, self.ax_dist)
		return data
	def draw_z(self, lab, px, py):
		if self.mz == 0:
			return
		r = (self.M_minscale + abs(self.mz/1000)/self.large_kNm) * lab.px_per_kN
		if self.mz > 0:
			lab.canv.create_arc(px-r, py-r, px+r, py+r, start=0, extent=270, style="arc",
				width=2, outline="red", tags=self.tag)
			T_points = self.arrow_head(px+self.arrow_offset, py+r, 1, 0)
		else:
			lab.canv.create_arc(px-r, py-r, px+r, py+r, start=90, extent=270, style="arc",
				width=2, outline="red", tags=self.tag)
			T_points = self.arrow_head(px+self.arrow_offset, py-r, 1, 0)
		lab.canv.create_polygon(T_points, fill="red", outline="red", tags=self.tag)
	def draw_xy(self, lab, px, py):
		Load.draw(self, lab, px, py)
		#Add second arrowhead
		Mmag_xy = math.sqrt(self.mx**2 + self.my**2)
		#Unit vector
		ux = self.mx / Mmag_xy
		uy = - self.my / Mmag_xy
		#Arrow head location
		h2x = px - ux * self.arrow_offset
		h2y = py - uy * self.arrow_offset
		T2_pts = self.arrow_head(h2x, h2y, ux, uy)
		lab.canv.create_polygon(T2_pts, fill="red", outline="red", tags=self.tag)
	def draw(self, *args):
		if self.mz != 0:
			self.draw_z(*args)
		if self.mx != 0 or self.my != 0:
			self.draw_xy(*args)
	
	#To Do: draw function (depends on if is in plane(mx,y) or out)
	#To Do: interface for placing applied Moments
