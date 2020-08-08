import math
import math_util as m_u


class Support:
	sup_types = {
		0: "Fixed",
		1: "Pin",
		2: "Slot",
		3: "Thrust",
	}
	img_w = 20
	def __init__(self, tag, ax_dist=0, th=0):
		self.tag = tag #tag on the canvas
		self.ax_dist = ax_dist
		self.th = th
	
	def __str__(self):
		return "{} support ({}) @axd={}".format(type(self).__name__, 
			self.tag, m_u.m_str(self.ax_dist))
	
	#The *args is used by Joints; they care which member is asking.
	def axd(self, *args):
		return self.ax_dist
	
	#Returns a binary tuple (x,y,th) for which directions the joint constrains
	def constraints(self):
		return (0,0,0)
	
	#Draws self on the given canvas
	def draw(self, canv, x, y):
		canv.create_oval(x-10,y-10,x+10,y+10)

class Fixed(Support):
	stype = 0
	def __init__(self, tag, ax_dist=0, th=0):
		super().__init__(tag, ax_dist, th)
	
	#def __str__(self):
	#	return "Fixed Support ({}) @axd={}".format(self.tag, self.ax_dist)
	
	def constraints(self):
		return (1,1,1)
	
	def draw(self, canv, x, y):
		o = (x,y)
		L0_pts = m_u.rot_pts(self.th, (x, y+2*self.img_w), 
			(x, y-2*self.img_w), origin=o).flatten()
		canv.create_line(*L0_pts, width=2, tags=self.tag)
		L1_pts = m_u.rot_pts(self.th, (x-self.img_w, y+2*self.img_w),
			(x, y+self.img_w), origin=o).flatten()
		canv.create_line(*L1_pts, tags=self.tag)
		L2_pts = m_u.rot_pts(self.th, (x-self.img_w, y+self.img_w), 
			(x, y), origin=o).flatten()
		canv.create_line(*L2_pts, tags=self.tag)
		L3_pts = m_u.rot_pts(self.th, (x-self.img_w, y), 
			(x, y-self.img_w), origin=o).flatten()
		canv.create_line(*L3_pts, tags=self.tag)
		L4_pts = m_u.rot_pts(self.th, (x-self.img_w, y-self.img_w), 
			(x, y-2*self.img_w), origin=o).flatten()
		canv.create_line(*L4_pts, tags=self.tag)

class Pin(Support):
	stype = 1
	def __init__(self, tag, ax_dist=0, th=0):
		super().__init__(tag, ax_dist, th)
	
	#def __str__(self):
	#	return "Pin Support ({}) @axd={}".format(self.tag, self.ax_dist)
	
	def constraints(self):
		return (1,1,0)
	
	def draw(self, canv, x, y):
		pts = [(x,y), (x-self.img_w,y-self.img_w/2), (x-self.img_w,y+self.img_w/2)]
		rotated = m_u.rot_pts(self.th, *pts, origin=(x,y)).flatten()
		canv.create_polygon(*rotated, width=2, fill="white", outline="black", tags=self.tag)

class Slot(Support):
	R = Support.img_w/2
	gap = 4
	def __init__(self, tag, ax_dist=0, th=0):
		super().__init__(tag, ax_dist, th)
	
	#def __str__(self):
	#	return "Slot Support ({}) @axd={}".format(self.tag, self.ax_dist)
	
	@property
	def stype(self):
		return 2
	
	def constraints(self):
		#abs() ? Or something else?
		return (m_u.eps_round(math.cos(math.radians(self.th))), 
			m_u.eps_round(math.sin(math.radians(self.th))), 0)
	
	def draw_track(self, canv, x, y):
		Rsinth = self.R*math.sin(math.radians(self.th))
		Rcosth = self.R*math.cos(math.radians(self.th))
		so0x = x - 2*Rsinth
		so0y = y - 2*Rcosth
		so1x = x + 2*Rsinth
		so1y = y + 2*Rcosth
		canv.create_arc(so0x-self.R, so0y-self.R, so0x+self.R, so0y+self.R, extent=180, 
			start=self.th, width=2, style="arc", tags=self.tag)
		canv.create_arc(so1x-self.R, so1y-self.R, so1x+self.R, so1y+self.R, extent=180, 
			start=180+self.th, width=2, style="arc", tags=self.tag)
		canv.create_line(so0x+Rcosth, so0y-Rsinth, so1x+Rcosth, so1y-Rsinth, 
			width=2, tags=self.tag)
		canv.create_line(so0x-Rcosth, so0y+Rsinth, so1x-Rcosth, so1y+Rsinth, 
			width=2, tags=self.tag)
	
	def draw(self, canv, x, y):
		self.draw_track(canv, x, y)
		#Pin
		canv.create_oval(x-self.R+self.gap, y-self.R+self.gap, x+self.R-self.gap, 
			y+self.R-self.gap, fill="black", tags=self.tag)

class Thrust(Slot):
	def __init__(self, tag, ax_dist=0, th=0):
		super().__init__(tag, ax_dist, th)
	
	#def __str__(self):
	#	return "Thrust Support ({}) @axd={}".format(self.tag, self.ax_dist)
	
	@property
	def stype(self):
		return 3
	
	def constraints(self):
		scon = super().constraints()
		return (scon[0], scon[1], 1)
	
	def draw(self, canv, x, y):
		self.draw_track(canv, x, y)
		#Pins
		Rsinth = self.R*math.sin(math.radians(self.th))
		Rcosth = self.R*math.cos(math.radians(self.th))
		p0x = x - Rsinth
		p0y = y - Rcosth
		p1x = x + Rsinth
		p1y = y + Rcosth
		Rpin = self.R - self.gap
		canv.create_oval(p0x-Rpin, p0y-Rpin, p0x+Rpin, p0y+Rpin, fill="black", tags=self.tag)
		canv.create_oval(p1x-Rpin, p1y-Rpin, p1x+Rpin, p1y+Rpin, fill="black", tags=self.tag)
