import tkinter as tk
import math
import numpy #For dot and cross products
#My files
from member import *
from region import *
from support import Support
from load import Load


class Lab:
	#Constants
	px_per_m = 160#scale
	#how long is a 1kN force arrow
	px_per_kN = 3.2#3.2 gives 160px=50kN
	dash_len = int(px_per_m/6)
	snap_dist = 15
	
	def __init__(self, main_frm, add_mem_bar, add_sup_bar, add_load_bar):
		self.c_wd = 800
		#Remember that when c_wd = 800 pixels, this corrosponds to x=[0,799]
		self.c_ht = 500
		#mode: 0=None, 1=add_mem, 2=add_sup, 3=add_load
		self.sel_mem_cb = None
		self.sel_txt = None
		self.add_mode = 0
		self.floating = None
		add_mem_bar.add_btn.config(command=self.toggle_mem_mode)
		self.add_mem_bar = add_mem_bar
		add_sup_bar.add_btn.config(command=self.toggle_sup_mode)
		self.add_sup_bar = add_sup_bar
		add_load_bar.add_btn.config(command=self.toggle_load_mode)
		self.add_load_bar = add_load_bar
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
	def snap_to_mem_ends(self, xp, yp):
		xc, yc = self.px_to_coords(xp,yp)
		side = 0
		closest_r = self.snap_dist / self.px_per_m
		closest = None
		for m in self.members:
			r0 = math.sqrt((m.x0-xc)**2 + (m.y0-yc)**2)
			if r0 < closest_r:
				closest_r = r0
				closest = m
				xp,yp = self.coords_to_px(m.x0, m.y0)
				side = 0
			r1 = math.sqrt((m.x1-xc)**2 + (m.y1-yc)**2)
			if r1 < closest_r:
				closest_r = r1
				closest = m
				xp,yp = self.coords_to_px(m.x1, m.y1)
				side = 1
		#xc, yc = self.px_to_coords(x,y)
		#closest_r_px = closest_r * self.px_per_m
		return ((xp, yp), closest, side)#((xp, yp), closest_r_px, side, closest)
	def snap_to_mem_axis(self, xp, yp):
		xc, yc = self.px_to_coords(xp,yp)
		closest_r = self.snap_dist / self.px_per_m
		closest = None
		v_comp = 0 #store the distance along the member (the comp.)
		for m in self.members:
			v_axis = m.v_axis()
			uv_axis = m.uv_axis()
			v_s0_P = (xc-m.x0, yc-m.y0)
			v_s1_P = (xc-m.x1, yc-m.y1)
			if numpy.dot(v_axis, v_s0_P) < 0: #P is before s0
				r = math.sqrt((m.x0-xc)**2 + (m.y0-yc)**2)
				if r < closest_r:
					closest_r = r
					closest = m
					xp,yp = self.coords_to_px(m.x0, m.y0)
			elif numpy.dot(v_axis, v_s1_P) > 0: #P is after s1
				r = math.sqrt((m.x1-xc)**2 + (m.y1-yc)**2)
				if r < closest_r:
					closest_r = r
					closest = m
					xp,yp = self.coords_to_px(m.x1, m.y1)
					v_comp = m.length
			else:
				r = abs(float(numpy.cross(uv_axis,v_s0_P)))
				if r < closest_r:
					closest_r = r
					closest = m
					s0 = numpy.array((m.x0,m.y0))
					v_comp = float(numpy.dot(v_s0_P,uv_axis))
					pf = s0 + v_comp*uv_axis
					xp,yp = self.coords_to_px(*pf)
		return ((xp,yp), closest, v_comp)
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
				if s != None:
					self.canv.move(s.tag, newx-m.oldx, newy-m.oldy)
			for l in m.loads:
				self.canv.move(l.tag, newx-m.oldx, newy-m.oldy)
		
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
				if not self.add_mem_bar.has_float_vals():
					return
				self.floating_x = x
				self.floating_y = y
				L = float(self.add_mem_bar.get_L())
				self.floating = self.canv.create_rectangle(*self.rect_mem_coords(x, y, L*self.px_per_m))
				return
		elif self.add_mode == 2:
			if self.floating == None:
				self.floating_x = x
				self.floating_y = y
				r = self.snap_dist
				self.floating = self.canv.create_oval(x-r,y-r,x+r,y+r)
				return
			else:
				(x,y),*_ = self.snap_to_mem_ends(x,y)
		elif self.add_mode == 3:
			if self.floating == None:
				if not self.add_load_bar.has_float_vals():
					return
				self.floating_x = x
				self.floating_y = y
				fl_tag = "fl_load"
				load = Load(0, fl_tag, *self.add_load_bar.get_P())
				load.draw(self.canv, x, y, self.px_per_kN)
				self.floating = fl_tag
				return
			else:
				(x,y),*_ = self.snap_to_mem_axis(x,y)
		else:
			return
		self.canv.move(self.floating, x-self.floating_x, y-self.floating_y)
		self.floating_x = x
		self.floating_y = y
	def mouse_click(self, event):
		#print(self.canv.find_all())
		if self.add_mode == 1:
			self.add_member(event)
		elif self.add_mode == 2:
			self.add_support(event)
		elif self.add_mode == 3:
			self.add_load(event)
		elif self.sel_mem_cb != None:
			self.sel_mem(event)
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
		(xp,yp), mem, side = self.snap_to_mem_ends(x,y)
		if mem != None:
			stype = self.add_sup_bar.get_sup_type()
			stag = 'sup_'+str(mem.img_ref)+'_s'+str(side)
			self.canv.delete(stag)#Delete the support img if it exists
			sup = Support(stype, stag)
			sup.draw(self.canv, xp, yp, mem.sup_dir(side))
			mem.sup[side] = sup
			self.set_add_mode(0)
	def add_load(self, event):
		x = event.x
		y = event.y
		(xp,yp), mem, ax_dist = self.snap_to_mem_axis(x,y)
		#(xc,yc) = self.px_to_coords(xp,yp)
		if mem != None:
			ltag = 'load_'+str(mem.img_ref)+'_d'+str(ax_dist)
			#self.canv.delete(ltag)
			ltype = 0#self.add_load_bar.something
			load = Load(ltype, ltag, *self.add_load_bar.get_P())
			load.ax_dist = ax_dist
			load.draw(self.canv, xp, yp, self.px_per_kN)
			mem.loads.append(load)
			self.set_add_mode(0)
			#popup_report(mem, 0)
			#mem.axial_stress()
	def sel_mem(self, event):
		x = event.x
		y = event.y
		(xp,yp), mem, _ = self.snap_to_mem_axis(x,y)
		if mem != None:
			#print("251: "+str(mem))
			cb = self.sel_mem_cb
			self.sel_mem_cb = None
			self.canv.delete(self.sel_txt)
			cb(mem)
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
		if mode == 3:
			if not self.add_load_bar.has_float_vals():
				return
		self.canv.delete(self.floating)
		self.floating = None
		self.add_mode = mode
	#callback is a function that accepts a reference to a member
	def select_mem(self, callback):
		self.canv.delete(self.sel_txt)
		self.sel_txt = self.canv.create_text(self.c_wd - 75, self.c_ht - 12, text="SELECT A MEMBER")
		self.sel_mem_cb = callback
	def del_mem(self, mem):
		for l in mem.loads:
			self.canv.delete(l.tag)
		for s in mem.sup:
			if s != None:
				self.canv.delete(s.tag)
		self.canv.delete(mem.img_ref)
		self.members.remove(mem)
	def clear_all(self):
		mem_copy = self.members.copy()
		for m in mem_copy:
			self.del_mem(m)
	def toggle_mem_mode(self):
		self.set_add_mode(0) if self.add_mode==1 else self.set_add_mode(1)
	def toggle_sup_mode(self):
		self.set_add_mode(0) if self.add_mode==2 else self.set_add_mode(2)
	def toggle_load_mode(self):
		self.set_add_mode(0) if self.add_mode==3 else self.set_add_mode(3)
	def eval_axial(self):
		self.select_mem(lambda m: self.popup_report(m,0))
	def del_mode(self):
		self.select_mem(lambda m: self.del_mem(m))
	def popup_report(self, mem, type):
		name = {
			0: "Axial Stress"
		}
		rep_text = mem.gen_report(type)
		popup = tk.Tk()
		popup.title(name[type]+" Report")
		popup.iconbitmap("../img/phys.ico")
		mem_lbl = tk.Label(popup, text=str(mem))
		mem_lbl.pack()
		rep_lbl = tk.Label(popup, text=rep_text)
		rep_lbl.pack()
		#Flash blue
		self.canv.itemconfig(mem.img_ref, outline="blue")
		def ol_black():
			self.canv.itemconfig(mem.img_ref, outline="black")
			popup.destroy()
		popup.protocol("WM_DELETE_WINDOW", ol_black)
		popup.mainloop()
