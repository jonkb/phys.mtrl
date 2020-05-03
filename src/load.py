import math

class Load:
	ah_width = 6
	load_types = {
		0: "Point",
		1: "Distributed"
	}
	def __init__(self, type, tag, xc, yc):
		self.type = type
		self.tag = tag
		self.xc = xc
		self.yc = yc
		self.ax_dist = 0 #To be set
	def draw(self, canv, px, py, px_per_kN):
		if self.type == 0:
			self.draw_pt(canv, px, py, px_per_kN)
	def get_comp(self):
		return (self.xc,self.yc)
	#Draw a point load. All args in pixels. The tip is at the point specified by (px,py)
	def draw_pt(self, canv, px, py, px_per_kN):
		pxc, pyc = self.kN_to_px(px_per_kN)
		canv.create_line(px-pxc,py-pyc,px,py, width=2, fill="red", tags=self.tag)
		L = math.sqrt(pxc**2+pyc**2)
		ux = pxc/L
		uy = pyc/L
		Trx = px - self.ah_width*3*ux + self.ah_width/2*uy
		Try = py - self.ah_width*3*uy - self.ah_width/2*ux
		Tlx = px - self.ah_width*3*ux - self.ah_width/2*uy
		Tly = py - self.ah_width*3*uy + self.ah_width/2*ux
		T_points = [px,py, Trx,Try, Tlx,Tly]
		canv.create_polygon(T_points, fill="red", outline="red", tags=self.tag)# width=1,
	def kN_to_px(self, px_per_kN):
		Pxp = self.xc*px_per_kN
		Pyp = -self.yc*px_per_kN
		return (Pxp, Pyp)
	