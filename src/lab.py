import tkinter as tk
import math
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#My files
from member import *
from region import *
import support
import joint
from load import *
from tk_wig import Txt_wig, Tk_rt
import pmfs


class Lab:
	#Constants
	version = "0.2.0"
	snap_dist = 15 #px
	dash_len = 30
	
	def __init__(self, main_frm, add_mem_bar, add_sup_bar, add_load_bar):
		self.c_wd = 756
		self.c_ht = 468
		
		add_mem_bar.add_btn.config(command=self.toggle_mem_mode)
		self.add_mem_bar = add_mem_bar
		add_sup_bar.add_btn.config(command=self.toggle_sup_mode)
		self.add_sup_bar = add_sup_bar
		add_load_bar.add_btn.config(command=self.toggle_load_mode)
		self.add_load_bar = add_load_bar
		
		self.canv = tk.Canvas(main_frm, width=self.c_wd, height=self.c_ht, bg='white')
		self.canv.config(borderwidth=0, highlightthickness=0)
		self.canv.pack(side = tk.BOTTOM, fill=tk.BOTH, expand=1)
		self.canv.bind("<Configure>", self.resize)
		self.canv.bind("<Motion>", self.mouse_moved)
		self.canv.bind("<Button-1>", self.mouse_click)
		self.bgbox = self.canv.create_rectangle(0,0,self.c_wd-1,self.c_ht-1)
		self.bggrid = [] #Switch this to a tag
		
		self.members = []
		self.popups = []
		self.reset_lab()
		
		#TEMP
		#self.main_frm = main_frm
	
	def to_xml(self):
		data = """
<lab version=\""""+self.version+"""\">
	<options>
		<c_wd>"""+str(self.c_wd)+"""</c_wd>
		<c_ht>"""+str(self.c_ht)+"""</c_ht>
		<px_per_m>"""+str(self.px_per_m)+"""</px_per_m>
		<subdivision>"""+str(self.subdivision)+"""</subdivision>
		<px_per_kN>"""+str(self.px_per_kN)+"""</px_per_kN>
	</options>
	<members>"""
		for mem in self.members:
			data += mem.to_xml()
		data += """
	</members>
</lab>"""
		return data
	
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
	
	#convert a vector in kN to px
	def kN_to_px(self, xc, yc):
		Pxp = xc*self.px_per_kN
		Pyp = -yc*self.px_per_kN
		return (Pxp, Pyp)
	
	def wtls(self):
		try:
			return self.mem_wtls.get()
		except AttributeError:# Was not set
			return True
	
	def smem(self):
		try:
			return self.show_mem.get()
		except AttributeError:# Was not set
			return True
	
	def ssup(self):
		try:
			return self.show_sup.get()
		except AttributeError:# Was not set
			return True
	
	def sld(self):
		try:
			return self.show_ld.get()
		except AttributeError:# Was not set
			return True
	
	def save(self):
		pmfs.save(self)
	
	def open(self):
		pmfs.open(self)
	
	def reset_lab(self):
		self.c_wd = 756
		#Remember that when c_wd = 800 pixels, this corrosponds to x=[0,799]
		self.c_ht = 468
		self.px_per_m = 360#scale
		self.subdivision = 5 #For grid lines smaller than 1m
		#how long is a 1kN force arrow
		self.px_per_kN = 20.0#4.0 gives 200px=50kN
		#A counter used to create unique tk tags
		self.tag_n = 0
		#Callback function to execute after selecting a member
		self.sel_mem_cb = None
		#The tk tag for the "SELECT A MEMBER" label
		self.sel_txt = None
		#mode: 0=None, 1=add_mem, 2=add_sup, 3=add_load
		self.add_mode = 0
		#The tk tag for a floating image following the cursor
		self.floating = None
		#Temp 1/2 of distributed load
		self.dlq0_temp = None
		self.cleanup()
		self.clear_all()
		self.redraw()
		self.members = []
		self.popups = []
	
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
	
	def snap_to_mem_axis(self, xp, yp, exclude=None):
		xc, yc = self.px_to_coords(xp,yp)
		closest_r = self.snap_dist / self.px_per_m
		closest = None
		#store axd along the member (the component of the projection)
		closest_axd = 0
		for m in self.members:
			if m is exclude:
				continue
			v_axis = m.v_axis()
			uv_axis = m.uv_axis()
			v_s0_P = (xc-m.x0, yc-m.y0)
			v_s1_P = (xc-m.x1, yc-m.y1)
			axd = np.dot(v_s0_P,uv_axis)
			if axd < 0: #P is before s0
				r = math.sqrt((m.x0-xc)**2 + (m.y0-yc)**2)
				if r < closest_r:
					closest_r = r
					closest = m
					xp,yp = self.coords_to_px(m.x0, m.y0)
					closest_axd = 0
			elif axd - m.length > 0: #P is after s1
				r = math.sqrt((m.x1-xc)**2 + (m.y1-yc)**2)
				if r < closest_r:
					closest_r = r
					closest = m
					xp,yp = self.coords_to_px(m.x1, m.y1)
					closest_axd = m.length
			else:
				r = abs(float(np.cross(uv_axis,v_s0_P)))
				if r < closest_r:
					closest_r = r
					closest = m
					s0 = np.array((m.x0,m.y0))
					closest_axd = axd
					pf = s0 + closest_axd*uv_axis
					xp,yp = self.coords_to_px(*pf)
		return ((xp,yp), closest, closest_axd)
	
	def snap_to_mem(self, xp, yp):
		xc, yc = self.px_to_coords(xp,yp)
		closest_r = self.snap_dist / self.px_per_m
		closest = None
		for m in self.members:
			v_axis = m.v_axis()
			uv_axis = m.uv_axis()
			uv_prp = np.array((-uv_axis[1], uv_axis[0]))
			P = np.array((xc, yc))
			v_s0_P = P - (m.x0, m.y0)
			axdP = np.dot(uv_axis, v_s0_P)
			axdP_s1 = axdP - m.length
			prpdP = np.dot(uv_prp, v_s0_P)
			hhdP = abs(prpdP) - m.half_h()
			O_ax = 0
			O_prp = 0
			if axdP < 0: #P is before s0
				O_ax = -axdP
			elif axdP_s1 > 0: #P is after s1
				O_ax = -axdP_s1
			if hhdP > 0: #P is further away from the axis than halfh
				O_prp = -hhdP * np.sign(prpdP)
			r = math.sqrt(O_ax**2 + O_prp**2)
			if r < closest_r:
				closest_r = r
				closest = m
				P += O_ax*uv_axis + O_prp*uv_prp
				xp,yp = self.coords_to_px(*P)
		return ((xp,yp), closest)
		
	def redraw(self, w=None, h=None, set_size=False):
		for m in self.members:
			m.oldx, m.oldy = self.coords_to_px(m.x0, m.y0)
		if w != None and h != None:
			self.c_wd, self.c_ht = w, h
			if set_size:
				#This next line works only if the window hasn't been manually resized yet.
				self.canv.config(width=w, height=h)
		#This can actually be way simplified by calculating a delta y just once.
		#Also I should make it change all the images if the scale changes.
		self.canv.coords(self.bgbox, 0, 0, self.c_wd-1,self.c_ht-1)
		for m in self.members:
			newx, newy = self.coords_to_px(m.x0, m.y0)
			self.canv.move(m.img_ref, newx-m.oldx, newy-m.oldy)
			for s in m.supports:
				self.canv.move(s.tag, newx-m.oldx, newy-m.oldy)
			for l in m.loads:
				self.canv.move(l.tag, newx-m.oldx, newy-m.oldy)
		self.redraw_grid()
	
	def redraw_grid(self):
		#Make it okay with no grid lines
		self.canv.delete(*self.bggrid)
		self.bggrid = []
		subdist = self.px_per_m / self.subdivision
		xn = int(self.c_wd / subdist)
		yn = int(self.c_ht / subdist)
		subdash = 4 #int(self.dash_len/2)
		#Make 1.0m^2 grid with subdivisions
		for i in range(1, xn+1):
			lnx = int(round(subdist*i))
			if i % self.subdivision == 0:
				self.bggrid.append(self.canv.create_line(
					lnx, self.c_ht-1, lnx, 0, fill="gray", dash=self.dash_len))
			else:
				self.bggrid.append(self.canv.create_line(
					lnx, self.c_ht-1, lnx, 0, fill="light gray", dash=subdash))
		for i in range(1, yn+1):
			lny = int(round(subdist*i))
			if i % self.subdivision == 0:
				self.bggrid.append(self.canv.create_line(0, self.c_ht-1-lny, 
					self.c_wd-1, self.c_ht-1-lny, fill="gray", dash=self.dash_len))
			else:
				self.bggrid.append(self.canv.create_line(0, self.c_ht-1-lny, 
					self.c_wd-1, self.c_ht-1-lny, fill="light gray", dash=subdash))
		for l in self.bggrid:
			self.canv.tag_lower(l)
	
	def resize(self, event):
		self.redraw(event.width, event.height)
	
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
				if self.add_sup_bar.is_jt():
					(x,y),*_ = self.snap_to_mem_axis(x,y)
				else: #CHANGE LATER ?
					(x,y),*_ = self.snap_to_mem_ends(x,y)
		elif self.add_mode == 3:
			if self.floating == None:
				if not self.add_load_bar.has_float_vals():
					return
				self.floating_x = x
				self.floating_y = y
				fl_tag = "fl_load"
				if self.add_load_bar.is_ds():
					if self.dlq0_temp == None:
						(xc, yc), *_ = self.add_load_bar.get_P()
					else:
						*_, (xc, yc) = self.add_load_bar.get_P()
				else:
					xc, yc = self.add_load_bar.get_P()
				load = Load(fl_tag, xc, yc)
				load.draw(self, x, y)
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
		#print("self.canv.config(width=300, height=300)")
		#self.canv.config(width=300, height=300)
		if self.add_mode == 1:
			self.add_member(event)
		elif self.add_mode == 2:
			if self.add_sup_bar.is_jt():
				self.add_joint(event)
			else:
				self.add_support(event)
		elif self.add_mode == 3:
			if self.add_load_bar.is_ds():
				self.add_distr_load(event)
			else:
				self.add_load(event)
		elif self.sel_mem_cb != None:
			self.sel_mem(event)
		#FOR TESTING
		else:
			self.sel_mem_cb = lambda m: print(m.static_eq())
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
		elif xsec == "annulus":
			ro = float(self.add_mem_bar.get_xparams()[0])/1000
			ri = float(self.add_mem_bar.get_xparams()[1])/1000
			xsec = Annulus(ro, ri)
		m = Member(matl, xsec, L)
		if not self.wtls():
			m.has_weight = True
		self.place_member(m, xp=x, yp=y, vh=self.add_mem_bar.get_vh())
		self.set_add_mode(0)
	
	def place_member(self, mem, xp=None, yp=None, xc=None, yc=None, vh=True):
		if vh == "V":
			vh = True
		elif vh == "H":
			vh = False
		if not self.wtls():
			mem.has_weight = True
		if not (xp is None) and not (yp is None):
			xc, yc = self.px_to_coords(xp, yp)
		elif not (xc is None) and not (yc is None):
			xp, yp = self.coords_to_px(xc, yc)
		else:
			return "Error: Must provide xp&yp or xc&yc"
		rgb = mem.material["color"]
		mem.place(xc, yc, vh)
		rcoords = self.rect_mem_coords(xp, yp, 
			mem.length*self.px_per_m, hh=mem.half_h(), vh=vh)
		mem.img_ref = self.canv.create_rectangle(*rcoords, fill=rgb)
		self.members.append(mem)
		#print("Added new " + str(m))
	
	def add_support(self, event):
		x = event.x
		y = event.y
		(xp,yp), mem, side = self.snap_to_mem_ends(x,y)
		if mem != None:
			stype = self.add_sup_bar.get_sup_type()
			self.place_support(mem, side, stype, xp, yp)
			self.set_add_mode(0)
	
	def place_support(self, mem, side, stype, xp=None, yp=None):
		if xp is None or yp is None:
			x0, y0, x1, y1 = mem.get_pos()
			if side == 0:
				xp, yp = self.coords_to_px(x0, y0)
			if side == 1:
				xp, yp = self.coords_to_px(x1, y1)
		axd = side*mem.length
		stag = "sup"+str(self.tag_n)+'_m'+str(mem.img_ref)
		self.tag_n += 1
		if stype == 0:
			sup = support.Fixed(stag, ax_dist=axd)
			sup.draw(self.canv, xp, yp, mem.sup_dir(side))
		elif stype == 1:
			sup = support.Pin(stag, ax_dist=axd)
			sup.draw(self.canv, xp, yp, mem.sup_dir(side))
		elif stype == 2:
			sup = support.Slot(stag, False, ax_dist=axd)
			sup.draw(self.canv, xp, yp)
		elif stype == 3:
			sup = support.Slot(stag, True, ax_dist=axd)
			sup.draw(self.canv, xp, yp)
		elif stype == 4:
			sup = support.Thrust(stag, False, ax_dist=axd)
			sup.draw(self.canv, xp, yp)
		elif stype == 5:
			sup = support.Thrust(stag, True, ax_dist=axd)
			sup.draw(self.canv, xp, yp)
		mem.supports.append(sup)
	
	def add_joint(self, event):
		x = event.x
		y = event.y
		_, mem, axd = self.snap_to_mem_axis(x,y)
		if mem != None:
			_, mem1, axd1 = self.snap_to_mem_axis(x,y, exclude=mem)
			if not (mem1 is None):
				#Include another case here for if they are in the same axis.
				#This could be doable after adding the snap grid.
				#Then m.intersection() will probably return None.
				intsx = mem.intersection(mem1)
				if not (intsx is None):
					((xc, yc), axd0, axd1) = intsx
					(xp, yp) = self.coords_to_px(xc, yc)
					jtype = self.add_sup_bar.get_sup_type()
					self.place_joint(mem, mem1, jtype, xp, yp, axd0=axd0, axd1=axd1)
					self.set_add_mode(0)
	
	def place_joint(self, m0, m1, jtype, xp=None, yp=None, axd0=None, axd1=None):
		if xp is None or yp is None:
			intsx = m0.intersection(m1)
			if intsx is None:
				return "Error: No intersection"
			((xc, yc), axd0, axd1) = intsx
			(xp, yp) = self.coords_to_px(xc, yc)
		jtag = "jt"+str(self.tag_n)+'_m0'+str(m0.img_ref)+'_m1'+str(m1.img_ref)
		self.tag_n += 1
		if jtype == 0:
			jt = joint.Fixed(jtag, m0, m1, axd0, axd1)
			jt.draw(self.canv, xp, yp)
		elif jtype == 1:
			jt = joint.Pin(jtag, m0, m1, axd0, axd1)
			jt.draw(self.canv, xp, yp)
		"""
		elif jtype == 2:
			jt = joint.Slot(jtag, False, m0, m1, axd0, axd1)
			jt.draw(self.canv, xp, yp)
		elif jtype == 3:
			jt = joint.Slot(jtag, True, m0, m1, axd0, axd1)
			jt.draw(self.canv, xp, yp)
		elif jtype == 4:
			jt = joint.Thrust(jtag, False, m0, m1, axd0, axd1)
			jt.draw(self.canv, xp, yp)
		elif jtype == 5:
			jt = joint.Thrust(jtag, True, m0, m1, axd0, axd1)
			jt.draw(self.canv, xp, yp)
		"""
		m0.joints.append(jt)
		m1.joints.append(jt)
		
	
	def add_load(self, event):
		x = event.x
		y = event.y
		(xp,yp), mem, axd = self.snap_to_mem_axis(x,y)
		#(xc,yc) = self.px_to_coords(xp,yp)
		if mem != None:
			Px, Py = self.add_load_bar.get_P()
			self.place_load(mem, Px, Py, axd, xp, yp)
			self.set_add_mode(0)
	
	def place_load(self, mem, Px, Py, axd, xp=None, yp=None):
		if xp is None or yp is None:
			x0, y0, *_ = mem.get_pos()
			xc, yc = np.array((x0, y0)) + axd*mem.uv_axis()
			xp, yp = self.coords_to_px(xc, yc)
		ltag = "ld_"+str(self.tag_n)
		self.tag_n += 1
		#ltag = 'load_'+str(mem.img_ref)+'_d'+str(ax_dist)
		#self.canv.delete(ltag)
		load = Load(ltag, Px, Py, axd)
		load.draw(self, xp, yp)
		mem.loads.append(load)
	
	def add_distr_load(self, event):
		x = event.x
		y = event.y
		(xp,yp), mem, ax_dist = self.snap_to_mem_axis(x,y)
		if mem != None:
			if self.dlq0_temp == None:
				(q0x, q0y), _ = self.add_load_bar.get_P()
				self.dlq0_temp = Load("dlq0_temp", q0x, q0y, ax_dist)
				self.dlq0_temp.draw(self, xp, yp)
				self.canv.delete(self.floating)
				self.floating = None
			else:
				(q0x, q0y), (q1x, q1y) = self.add_load_bar.get_P()
				axd0 = self.dlq0_temp.ax_dist
				if mem.is_vert():
					yp += (ax_dist - axd0)*self.px_per_m
				else:
					xp -= (ax_dist - axd0)*self.px_per_m
				self.place_distr_load(mem, q0x, q0y, axd0, q1x, q1y, ax_dist, xp, yp)
				self.canv.delete(self.dlq0_temp.tag)
				self.dlq0_temp = None
				self.set_add_mode(0)
	
	def place_distr_load(self, mem, q0x, q0y, axd0, q1x, q1y, axd1, xp=None, yp=None):
		if xp is None or yp is None:
			x0, y0, *_ = mem.get_pos()
			xc, yc = np.array((x0, y0)) + axd0*mem.uv_axis()
			xp, yp = self.coords_to_px(xc, yc)
		isv = mem.is_vert()
		dltag = "dl_"+str(self.tag_n)
		self.tag_n += 1
		#+str(mem.img_ref)+"_d"+str(self.dlq0_temp.ax_dist)+"-"+str(ax_dist)
		dl = Distr_Load(dltag, q0x, q0y, axd0, q1x, q1y, axd1, isv)
		mem.loads.append(dl)
		dl.draw(self, xp, yp)
	
	def sel_mem(self, event):
		x = event.x
		y = event.y
		(xp,yp), mem = self.snap_to_mem(x,y)
		if mem != None:
			#print("251: "+str(mem))
			cb = self.sel_mem_cb
			self.sel_mem_cb = None
			self.canv.delete(self.sel_txt)
			cb(mem)
	
	def rect_mem_coords(self, x, y, L_px, hh=None, vh=None):
		if hh is None:
			hh = self.add_mem_bar.half_h()
		if hh == "NaN":
			return "NaN"
		if vh is None:
			vh = self.add_mem_bar.get_vh()
		hh *= self.px_per_m
		if vh:
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
	def set_sel_mem_cb(self, callback):
		self.canv.delete(self.sel_txt)
		self.sel_txt = self.canv.create_text(self.c_wd - 75, self.c_ht - 12, text="SELECT A MEMBER")
		self.sel_mem_cb = callback
	
	def del_mem(self, mem):
		for l in mem.loads:
			self.canv.delete(l.tag)
		for s in mem.supports:
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
		#print(354, self.add_mode)
		self.set_add_mode(0) if self.add_mode==3 else self.set_add_mode(3)
	
	def toggle_wtls(self):
		for m in self.members:
			m.has_weight = not self.wtls()
	
	def toggle_smem(self):
		if self.smem():
			self.add_mem_bar.tb_frm.pack(side=tk.TOP, fill=tk.X)
		else:
			self.add_mem_bar.tb_frm.pack_forget()
	
	def toggle_ssup(self):
		if self.ssup():
			self.add_sup_bar.tb_frm.pack(side=tk.TOP, fill=tk.X)
		else:
			self.add_sup_bar.tb_frm.pack_forget()
	
	def toggle_sld(self):
		if self.sld():
			self.add_load_bar.tb_frm.pack(side=tk.TOP, fill=tk.X)
		else:
			self.add_load_bar.tb_frm.pack_forget()
	
	def show_allt(self):
		if not self.smem():
			self.add_mem_bar.tb_frm.pack(side=tk.TOP, fill=tk.X)
			self.show_mem.set(True)
		if not self.ssup():
			self.add_sup_bar.tb_frm.pack(side=tk.TOP, fill=tk.X)
			self.show_sup.set(True)
		if not self.sld():
			self.add_load_bar.tb_frm.pack(side=tk.TOP, fill=tk.X)
			self.show_ld.set(True)
	
	def hide_allt(self):
		if self.smem():
			self.add_mem_bar.tb_frm.pack_forget()
			self.show_mem.set(False)
		if self.ssup():
			self.add_sup_bar.tb_frm.pack_forget()
			self.show_sup.set(False)
		if self.sld():
			self.add_load_bar.tb_frm.pack_forget()
			self.show_ld.set(False)
	
	def eval_report(self, rtype):
		self.set_sel_mem_cb(lambda m: self.popup_report(m,rtype))
	
	def del_mode(self):
		self.set_sel_mem_cb(lambda m: self.del_mem(m))
	
	def popup_report(self, mem, rtype):
		popup = Tk_rt(mem.eval_names[rtype]+" Report")
		loading_lbl = tk.Label(popup, text="LOADING", padx=48, pady=24)
		loading_lbl.pack()
		def add_report():
			report = mem.gen_report(rtype)
			loading_lbl.destroy()
			mem_lbl = tk.Label(popup, text=str(mem))
			mem_lbl.pack()
			if isinstance(report, str):
				#This can be either a text-only report or an error message
				rep_lbl = Txt_wig(popup, report)
				rep_lbl.packslf()
			elif isinstance(report, tuple):
				rep_text, fig = report
				Txt_wig(popup, rep_text).packslf()
				figcanv = FigureCanvasTkAgg(fig, popup)
				figcanv.get_tk_widget().pack(fill=tk.BOTH, expand=1)
			elif isinstance(report, list):
				tk.Label(popup, text="Choose which report to show:").pack()
				rep_names = []
				rep_widgets = []
				for rep in report:
					rep_name, rep_text, fig = rep
					rep_names.append(rep_name)
					rep_lbl = Txt_wig(popup, rep_text)
					figcanv = FigureCanvasTkAgg(fig, popup)
					rep_canv = figcanv.get_tk_widget()
					#rep_canv.pack(fill=tk.BOTH, expand=1)
					rep_widgets.append((rep_lbl, rep_canv))
				#Pulldown to choose which report to show
				rep_option = tk.StringVar(popup)
				rep_option.set(rep_names[0])
				matl_option = tk.OptionMenu(popup, rep_option, *rep_names)
				matl_option.config(width=16)
				matl_option.pack()
				popup.current_rep = 0
				rep_widgets[popup.current_rep][0].packslf()
				rep_widgets[popup.current_rep][1].pack(fill=tk.BOTH, expand=1)
				def switch_rep(*args):
					rep_widgets[popup.current_rep][0].pack_forget()
					rep_widgets[popup.current_rep][1].pack_forget()
					popup.current_rep = rep_names.index(rep_option.get())
					rep_widgets[popup.current_rep][0].packslf()
					rep_widgets[popup.current_rep][1].pack(fill=tk.BOTH, expand=1)
				rep_option.trace("w", switch_rep)
			else:
				tk.Label(popup, text="Error: unknown report data", justify=tk.LEFT).pack()
		
		#Flash blue
		self.canv.itemconfig(mem.img_ref, outline="blue")
		mem.popups += 1
		def del_pop():
			mem.popups -= 1
			#This way, if you have multiple reports on the same member, it stays blue when closing one
			if mem.popups <= 0:
				self.canv.itemconfig(mem.img_ref, outline="black")
			self.popups.remove(popup)
			popup.destroy()
		popup.protocol("WM_DELETE_WINDOW", del_pop)
		self.popups.append(popup)
		popup.after(200, add_report)
		popup.mainloop()
	
	#Open window to edit px_per_m and px_per_kN
	def edit_scale(self):
		popup = Tk_rt("Edit Lab Scale")
		popup.grid_columnconfigure(0, minsize=120)
		popup.grid_columnconfigure(1, minsize=120)
		popup.grid_rowconfigure(4, pad=16)
		
		ppm_lbl = tk.Label(popup, text="px per m: ")
		ppm_lbl.grid(row=0, column=0, sticky=tk.E)
		ppm_e = tk.Entry(popup, width=8)
		ppm_e.grid(row=0, column=1, sticky=tk.W)
		ppm_e.insert(0, self.px_per_m)
		ppk_lbl = tk.Label(popup, text="px per kN: ")
		ppk_lbl.grid(row=1, column=0, sticky=tk.E)
		ppk_e = tk.Entry(popup, width=8)
		ppk_e.grid(row=1, column=1, sticky=tk.W)
		ppk_e.insert(0, self.px_per_kN)
		sd_lbl = tk.Label(popup, text="grid lines per m: ")
		sd_lbl.grid(row=2, column=0, sticky=tk.E)
		sd_e = tk.Entry(popup, width=8)
		sd_e.grid(row=2, column=1, sticky=tk.W)
		sd_e.insert(0, self.subdivision)
		
		note = "Note: For now, since resizing with members on the canvas does not work "
		note += "properly, hitting save will also remove any members from the screen."
		note_lbl = tk.Label(popup, text=note, wraplength=300)
		note_lbl.grid(row=3, column=0, columnspan=2)
		
		def set_scale():
			try:
				ppm = int(ppm_e.get())
				self.px_per_m = ppm
			except ValueError:
				ppm_e.delete(0, tk.END)
				ppm_e.insert(0, self.px_per_m)
			try:
				ppk = float(ppk_e.get())
				self.px_per_kN = ppk
			except ValueError:
				ppk_e.delete(0, tk.END)
				ppk_e.insert(0, self.px_per_kN)
			try:
				#Limit to Z or N
				sd = int(sd_e.get())
				self.subdivision = sd
			except ValueError:
				sd_e.delete(0, tk.END)
				sd_e.insert(0, self.subdivision)
			self.clear_all()
			self.redraw()
		save_btn = tk.Button(popup, text="Save", command=set_scale)
		save_btn.grid(row=4, column=0, columnspan=2, ipadx= 16, ipady=2)
		def del_pop():
			self.popups.remove(popup)
			popup.destroy()
		popup.protocol("WM_DELETE_WINDOW", del_pop)
		self.popups.append(popup)
		popup.mainloop()
	
	def cleanup(self):
		for w in self.popups:
			w.destroy()
		plt.close("all")
