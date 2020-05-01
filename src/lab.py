import tkinter as tk
import math
#My files
from member import *
from region import *


class Lab:
	#Constants
	px_per_m = 200#scale
	dash_len = int(px_per_m/6)
	snap_dist = 15
	
	def __init__(self, main_frm, add_mem_bar, add_sup_bar):
		self.c_wd = 800
		self.c_ht = 500
		#Remember that when c_wd = 800 pixels, this corrosponds to x=[0,799]
		add_mem_bar.add_btn.config(command=self.toggle_mem_mode)
		self.add_mem_bar = add_mem_bar
		self.floating_mem = None
		self.mem_mode = False
		add_sup_bar.add_btn.config(command=self.toggle_sup_mode)
		self.add_sup_bar = add_sup_bar
		self.floating_sup = None
		self.sup_mode = False
		self.members = []
		
		self.canv = tk.Canvas(main_frm, width=self.c_wd, height=self.c_ht, bg='white')
		self.canv.config(borderwidth=0, highlightthickness=0)
		self.canv.pack(side = tk.BOTTOM, fill=tk.BOTH, expand=1)
		self.canv.bind("<Configure>", self.resize)
		self.canv.bind("<Motion>", self.mouse_moved)
		self.canv.bind("<Button-1>", self.mouse_click)
		self.bgbox = self.canv.create_rectangle(0,0,self.c_wd-1,self.c_ht-1)
		self.bggrid = []
		#Make 1.0m^2 grid
		for linex in range(self.px_per_m, self.c_wd, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				linex, self.c_ht-1, linex, 0, fill="gray", dash=(self.dash_len,)))
		for liney in range(self.px_per_m, self.c_ht, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				0, self.c_ht-1-liney, self.c_wd-1, self.c_ht-1-liney, fill="gray", dash=(self.dash_len,)))
	
	#convert x,y coordinates in meters to their place on the canvas in pixels
	def coords_to_px(self, xc, yc):
		xp = xc*self.px_per_m
		yp = self.c_ht - yc*self.px_per_m
		return (xp,yp)
	#convert a location on the canvas in pixels to x,y coordinates in meters
	def px_to_coords(self, xp, yp):
		xc = xp / self.px_per_m
		yc = (self.c_ht - yp) / self.px_per_m
		return (xc,yc)
	#Snap x,y (in px) to the ends of existing members
	def snap_to_mem(self, xp, yp):
		xc, yc = self.px_to_coords(xp,yp)
		#Direction 0,1,2,3 --> Up, Left, Down, Right -- Pointing away from member
		dir = 0
		closest_r = self.snap_dist / self.px_per_m
		for m in self.members:
			r0 = math.sqrt((m.x0-xc)**2 + (m.y0-yc)**2)
			if r0 < closest_r:
				xp,yp = self.coords_to_px(m.x0, m.y0)
				closest_r = r0
				if m.get_vh():
					dir = 2
				else:
					dir = 1
			r1 = math.sqrt((m.x1-xc)**2 + (m.y1-yc)**2)
			if r1 < closest_r:
				xp,yp = self.coords_to_px(m.x1, m.y1)
				closest_r = r1
				if m.get_vh():
					dir = 0
				else:
					dir = 3
		#xc, yc = self.px_to_coords(x,y)
		closest_r_px = closest_r * self.px_per_m
		return ((xp, yp),closest_r_px,dir)
	def resize(self, event):
		for m in self.members:
			m.oldx, m.oldy = self.coords_to_px(m.x0, m.y0)
		self.c_wd, self.c_ht = event.width, event.height
		self.canv.coords(self.bgbox, 0, 0, self.c_wd-1,self.c_ht-1)
		for m in self.members:
			newx, newy = self.coords_to_px(m.x0, m.y0)
			self.canv.move(m.img_ref, newx-m.oldx, newy-m.oldy)
		
		self.canv.delete(*self.bggrid)
		#Make 1.0m^2 grid
		for linex in range(self.px_per_m, self.c_wd, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				linex, self.c_ht-1, linex, 0, fill="gray", dash=(self.dash_len,)))
		for liney in range(self.px_per_m, self.c_ht, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				0, self.c_ht-1-liney, self.c_wd-1, self.c_ht-1-liney, fill="gray", dash=(self.dash_len,)))
	def mouse_moved(self, event):
		x = event.x
		y = event.y
		if self.mem_mode:
			if self.floating_mem == None:
				self.floating_x = x
				self.floating_y = y
				L = float(self.add_mem_bar.get_L())
				self.floating_mem = self.canv.create_rectangle(*self.rect_mem_coords(x, y, L*self.px_per_m))
			else:
				self.canv.move(self.floating_mem, x-self.floating_x, y-self.floating_y)
				self.floating_x = x
				self.floating_y = y
		elif self.sup_mode:
			if self.floating_sup == None:
				self.floating_x = x
				self.floating_y = y
				r = 15
				self.floating_sup = self.canv.create_oval(x-r,y-r,x+r,y+r)
			else:
				(x,y),*_ = self.snap_to_mem(x,y)
				self.canv.move(self.floating_sup, x-self.floating_x, y-self.floating_y)
				self.floating_x = x
				self.floating_y = y
	def mouse_click(self, event):
		if self.mem_mode:
			self.add_member(event)
		#print(canv.find_all())
		if self.sup_mode:
			self.add_support(event)
	def add_member(self, event):
		x = event.x
		y = event.y
		try:
			L = float(self.add_mem_bar.get_L())
		except:#Flash the L field red here?
			return
		matl = getattr(Materials,self.add_mem_bar.get_matl())
		rgb = matl["color"]
		xsec = self.add_mem_bar.get_xsec()
		if xsec == "circle":
			try:
				r = float(self.add_mem_bar.get_xparams()[0])/1000 #mm->m
			except:
				return
			xsec = Circle(r)
		if xsec == "rectangle":
			try:
				b = float(self.add_mem_bar.get_xparams()[0])/1000 #mm->m
				h = float(self.add_mem_bar.get_xparams()[1])/1000 
			except:
				return
			xsec = Rectangle(b, h)
		if xsec == "I-beam":
			try:
				d = float(self.add_mem_bar.get_xparams()[0])/1000 #mm->m
				w = float(self.add_mem_bar.get_xparams()[1])/1000 
				tf = float(self.add_mem_bar.get_xparams()[2])/1000 
				tw = float(self.add_mem_bar.get_xparams()[3])/1000 
			except:
				return
			xsec = W_F_I(d, w, tf, tw)
		m = Member(matl, xsec, L)
		xc,yc = self.px_to_coords(x,y)
		m.place(xc, yc, self.add_mem_bar.get_vh())
		m.img_ref = self.canv.create_rectangle(*self.rect_mem_coords(x, y, L*self.px_per_m), fill=rgb)
		self.members.append(m)
		#print("Added new " + str(m))
		self.set_mem_mode(False)
	def add_support(self, event):
		x = event.x
		y = event.y
		(xp,yp), closest_r, dir = self.snap_to_mem(x,y)
		if closest_r < self.snap_dist:
			self.set_sup_mode(False)
	def rect_mem_coords(self, x,y,L_px):
		hh = self.add_mem_bar.half_h()
		if hh=="NaN":
			return "NaN"
		hh *= self.px_per_m
		if(self.add_mem_bar.get_vh()):
			x1 = x - hh
			y1 = y
			x2 = x + hh
			y2 = y - L_px
		else:
			x1 = x
			y1 = y - hh
			x2 = x + L_px
			y2 = y + hh
		return (x1, y1, x2, y2)
	def set_mem_mode(self, on_off):
		if on_off:
			try:
				L = float(self.add_mem_bar.get_L())
			except:#Flash the L field red here?
				return
			self.mem_mode = True
			if self.sup_mode:
				self.set_sup_mode(False)
		else:
			self.mem_mode = False
			self.canv.delete(self.floating_mem)
			self.floating_mem = None
	def set_sup_mode(self, on_off):
		if on_off:
			self.sup_mode = True
			if self.mem_mode:
				self.set_mem_mode(False)
		else:
			self.sup_mode = False
			self.canv.delete(self.floating_sup)
			self.floating_sup = None
	def toggle_mem_mode(self):
		self.set_mem_mode(False) if self.mem_mode else self.set_mem_mode(True)
	def toggle_sup_mode(self):
		self.set_sup_mode(False) if self.sup_mode else self.set_sup_mode(True)
