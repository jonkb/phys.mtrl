import tkinter as tk
import math
#My files
from member import *
from region import *
from support import Support


class Lab:
	#Constants
	px_per_m = 180#scale
	dash_len = int(px_per_m/6)
	snap_dist = 15
	
	def __init__(self, main_frm, add_mem_bar, add_sup_bar):
		self.c_wd = 800
		#Remember that when c_wd = 800 pixels, this corrosponds to x=[0,799]
		self.c_ht = 500
		#mode: 0=None, 1=add_mem, 2=add_sup
		self.add_mode = 0#
		self.floating = None
		add_mem_bar.add_btn.config(command=self.toggle_mem_mode)
		self.add_mem_bar = add_mem_bar
		add_sup_bar.add_btn.config(command=self.toggle_sup_mode)
		self.add_sup_bar = add_sup_bar
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
		side = 0
		closest_r = self.snap_dist / self.px_per_m
		closest = None
		for m in self.members:
			r0 = math.sqrt((m.x0-xc)**2 + (m.y0-yc)**2)
			if r0 < closest_r:
				xp,yp = self.coords_to_px(m.x0, m.y0)
				closest_r = r0
				closest = m
				side = 0
			r1 = math.sqrt((m.x1-xc)**2 + (m.y1-yc)**2)
			if r1 < closest_r:
				xp,yp = self.coords_to_px(m.x1, m.y1)
				closest_r = r1
				closest = m
				side = 1
		#xc, yc = self.px_to_coords(x,y)
		closest_r_px = closest_r * self.px_per_m
		return ((xp, yp), closest_r_px, side, closest)
	def resize(self, event):
		#This can actually be way simplified by calculating a delta y just once.
		for m in self.members:
			m.oldx, m.oldy = self.coords_to_px(m.x0, m.y0)
		self.c_wd, self.c_ht = event.width, event.height
		self.canv.coords(self.bgbox, 0, 0, self.c_wd-1,self.c_ht-1)
		for m in self.members:
			newx, newy = self.coords_to_px(m.x0, m.y0)
			self.canv.move(m.img_ref, newx-m.oldx, newy-m.oldy)
			for s in m.sup:
				if not s == None:
					self.canv.move(s.tag, newx-m.oldx, newy-m.oldy)
		
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
		if self.add_mode == 1:
			if self.floating == None:
				self.floating_x = x
				self.floating_y = y
				L = float(self.add_mem_bar.get_L())
				self.floating = self.canv.create_rectangle(*self.rect_mem_coords(x, y, L*self.px_per_m))
			else:
				self.canv.move(self.floating, x-self.floating_x, y-self.floating_y)
				self.floating_x = x
				self.floating_y = y
		elif self.add_mode == 2:
			if self.floating == None:
				self.floating_x = x
				self.floating_y = y
				r = 15
				self.floating = self.canv.create_oval(x-r,y-r,x+r,y+r)
			else:
				(x,y),*_ = self.snap_to_mem(x,y)
				self.canv.move(self.floating, x-self.floating_x, y-self.floating_y)
				self.floating_x = x
				self.floating_y = y
	def mouse_click(self, event):
		if self.add_mode == 1:
			self.add_member(event)
		elif self.add_mode == 2:
			self.add_support(event)
		print(self.canv.find_all())
	def add_member(self, event):
		if not self.add_mem_bar.has_float_vals():
			return
		x = event.x
		y = event.y
		L = float(self.add_mem_bar.get_L())
		matl = getattr(Materials,self.add_mem_bar.get_matl())
		rgb = matl["color"]
		xsec = self.add_mem_bar.get_xsec()
		if xsec == "circle":
			r = float(self.add_mem_bar.get_xparams()[0])/1000 #mm->m
			xsec = Circle(r)
		elif xsec == "rectangle":
			b = float(self.add_mem_bar.get_xparams()[0])/1000 #mm->m
			h = float(self.add_mem_bar.get_xparams()[1])/1000 
			xsec = Rectangle(b, h)
		elif xsec == "I-beam":
			d = float(self.add_mem_bar.get_xparams()[0])/1000 #mm->m
			w = float(self.add_mem_bar.get_xparams()[1])/1000 
			tf = float(self.add_mem_bar.get_xparams()[2])/1000 
			tw = float(self.add_mem_bar.get_xparams()[3])/1000 
			xsec = W_F_I(d, w, tf, tw)
		m = Member(matl, xsec, L)
		xc,yc = self.px_to_coords(x,y)
		m.place(xc, yc, self.add_mem_bar.get_vh())
		m.img_ref = self.canv.create_rectangle(*self.rect_mem_coords(x, y, L*self.px_per_m), fill=rgb)
		self.members.append(m)
		#print("Added new " + str(m))
		self.set_add_mode(0)
	def add_support(self, event):
		x = event.x
		y = event.y
		(xp,yp), closest_r, side, mem = self.snap_to_mem(x,y)
		if closest_r < self.snap_dist:
			stype = self.add_sup_bar.get_sup_type()
			stag = 'sup_'+str(mem.img_ref)+'s'+str(side)
			self.canv.delete(stag)#Delete the support img if it exists
			sup = Support(stype, stag)
			sup.draw(self.canv, xp, yp, mem.sup_dir(side))
			mem.sup[side] = sup
			self.set_add_mode(0)
	def rect_mem_coords(self, x, y, L_px):
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
	#Mode: 0=None, 1=mem, 2=sup
	def set_add_mode(self, mode):
		if mode == self.add_mode:
			return
		if mode == 1:
			if not self.add_mem_bar.has_float_vals():
				return
		self.canv.delete(self.floating)
		self.floating = None
		self.add_mode = mode
	def toggle_mem_mode(self):
		self.set_add_mode(0) if self.add_mode==1 else self.set_add_mode(1)
	def toggle_sup_mode(self):
		self.set_add_mode(0) if self.add_mode==2 else self.set_add_mode(2)
