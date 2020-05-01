import tkinter as tk
#My files
from member import *
from region import *


class Lab:
	#Constants
	px_per_m = 200#scale
	dash_len = int(px_per_m/6)
	
	def __init__(self, main_frm, add_mem_bar):
		self.c_wd = 800
		self.c_ht = 500
		#Remember that when c_wd = 800 pixels, this corrosponds to x=[0,799]
		self.add_mem_bar = add_mem_bar
		self.floating_mem = None
		
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
				linex, 0, linex, self.c_ht-1, fill="gray", dash=(self.dash_len,)))
		for liney in range(self.px_per_m, self.c_ht, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				0, self.c_ht-1-liney, self.c_wd-1, self.c_ht-1-liney, fill="gray", dash=(self.dash_len,)))
	
	def resize(self, event):
		self.c_wd, self.c_ht = event.width, event.height
		self.canv.coords(self.bgbox, 0, 0, self.c_wd-1,self.c_ht-1)
		#? Add something to move everything else, so it's all anchored to the bottom left
		self.canv.delete(*self.bggrid)
		#Make 1.0m^2 grid
		for linex in range(self.px_per_m, self.c_wd, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				linex, 0, linex, self.c_ht-1, fill="gray", dash=(self.dash_len,)))
		for liney in range(self.px_per_m, self.c_ht, self.px_per_m):
			self.bggrid.append(self.canv.create_line(
				0, self.c_ht-1-liney, self.c_wd-1, self.c_ht-1-liney, fill="gray", dash=(self.dash_len,)))
	def mouse_moved(self, event):
		if self.add_mem_bar.add_mode:
			x = event.x
			y = event.y
			if self.floating_mem == None:
				self.floating_x = x
				self.floating_y = y
				L = float(self.add_mem_bar.add_L_entry.get())
				self.floating_mem = self.canv.create_rectangle(*self.rect_coords(x, y, L*self.px_per_m))
			else:
				self.canv.move(self.floating_mem, x-self.floating_x, y-self.floating_y)
				self.floating_x = x
				self.floating_y = y
	def mouse_click(self, event):
		if self.add_mem_bar.add_mode:
			self.add_member(event)
		#print(canv.find_all())
	def add_member(self, event):
		x = event.x
		y = event.y
		try:
			L = float(self.add_mem_bar.add_L_entry.get())
		except:#Flash the L field red here?
			return
		matl = getattr(Materials,self.add_mem_bar.add_matl.get())
		rgb = matl["color"]
		xsec = self.add_mem_bar.add_xsec.get()
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
		m = Member(matl, xsec, 1)#Store this!!!
		print("Added new " + str(m))
		self.canv.create_rectangle(*self.rect_coords(x, y, L*self.px_per_m), fill=rgb)
		self.add_mem_bar.toggle_add()
		self.canv.delete(self.floating_mem)
		self.floating_mem = None
	def rect_coords(self, x,y,L_px):
		hh = self.add_mem_bar.half_h()
		if hh=="NaN":
			return "NaN"
		hh *= self.px_per_m
		if(self.add_mem_bar.add_vh.get()):
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



