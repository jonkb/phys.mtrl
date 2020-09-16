import tkinter as tk
import math

from load import Load
from member import *
from region import *
from support import Support


num_e_wid = 8 #Standard width for a numerical tk.Entry widget.
option_wid = 10 #Standard width for tk.OptionMenu pulldown. 

#Toolbar to add a new member
class Add_mem:
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		self.tb_frm.grid_rowconfigure(1,weight=1)
		
		#First label
		tb_lbl = tk.Label(self.tb_frm, text="Add New\nMember")
		tb_lbl.grid(row=0, column=0, rowspan=2)
		
		#Choose material label
		matl_lbl = tk.Label(self.tb_frm, text="Material:")
		matl_lbl.grid(row=0, column=1)
		#Choose material pulldown
		self.matl = tk.StringVar(self.tb_frm)
		self.matl.set(Materials.materials[0])
		matl_option = tk.OptionMenu(self.tb_frm, self.matl, *Materials.materials)
		matl_option.config(width=option_wid)
		matl_option.grid(row=1, column=1)
		
		#Choose xsection label
		xsec_lbl = tk.Label(self.tb_frm, text="Cross Section Type:")
		xsec_lbl.grid(row=0, column=2)
		#Choose xsection pulldown
		self.xsec = tk.StringVar(self.tb_frm)
		self.xsec.set(Region.regions[0])
		self.xsec.trace("w", self.update_xparam)
		xsec_option = tk.OptionMenu(self.tb_frm, self.xsec, *Region.regions)
		xsec_option.config(width=option_wid)
		xsec_option.grid(row=1, column=2)
		
		#Frame that adjusts itself to the chosen xsection to have the needed parameters
		self.xparam_frm = tk.Frame(self.tb_frm)
		self.xparam_frm.config(borderwidth=2, relief=tk.SUNKEN)
		self.xparam_frm.grid(row=0, column=3, rowspan=2, sticky=tk.N+tk.S)
		self.xparam_frm.grid_rowconfigure(1,weight=1)
		self.xparam_entries = []
		self.update_xparam()
		
		#Length label
		L_lbl = tk.Label(self.tb_frm, text="Length (m):")
		L_lbl.grid(row=0, column=4)
		#Length entry
		self.L_entry = tk.Entry(self.tb_frm)
		self.L_entry.config(width=num_e_wid)
		self.L_entry.grid(row=1, column=4)
		
		#Angle label
		th_lbl = tk.Label(self.tb_frm, text="Angle (\u00B0):")
		th_lbl.grid(row=0, column=5)
		#Angle entry
		self.th_entry = tk.Entry(self.tb_frm)
		self.th_entry.config(width=num_e_wid)
		self.th_entry.grid(row=1, column=5)
		self.th_entry.insert(0, "0")
		
		#Button to add the new member
		self.add_btn = tk.Button(self.tb_frm, text="Add")
		#self.add_btn.config(command=self.toggle_add)
		self.add_btn.grid(row=0, column=6, padx=2, pady=2, ipadx=8, rowspan=2, sticky=tk.N+tk.S)
	
	def update_xparam(self, *args):
		region = self.xsec.get()
		self.xparam_entries.clear()
		for widget in self.xparam_frm.winfo_children():
			widget.destroy()
		if(region == "circle"):
			rad_lbl = tk.Label(self.xparam_frm, text="Radius (mm):")
			rad_lbl.grid(row=0, column=0)
			rad_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			rad_entry.grid(row=1, column=0)
			self.xparam_entries.append(rad_entry)
		if(region == "rectangle"):
			b_lbl = tk.Label(self.xparam_frm, text="Base (mm):")
			b_lbl.grid(row=0, column=0)
			b_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			b_entry.grid(row=1, column=0)
			self.xparam_entries.append(b_entry)
			h_lbl = tk.Label(self.xparam_frm, text="Height (mm):")
			h_lbl.grid(row=0, column=1)
			h_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			h_entry.grid(row=1, column=1)
			self.xparam_entries.append(h_entry)
		if(region == "I-beam"):
			d_lbl = tk.Label(self.xparam_frm, text="Depth (mm):")
			d_lbl.grid(row=0, column=0)
			d_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			d_entry.grid(row=1, column=0)
			self.xparam_entries.append(d_entry)
			w_lbl = tk.Label(self.xparam_frm, text="Width (mm):")
			w_lbl.grid(row=0, column=1)
			w_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			w_entry.grid(row=1, column=1)
			self.xparam_entries.append(w_entry)
			tf_lbl = tk.Label(self.xparam_frm, text="Flange t (mm):")
			tf_lbl.grid(row=0, column=2)
			tf_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			tf_entry.grid(row=1, column=2)
			self.xparam_entries.append(tf_entry)
			tw_lbl = tk.Label(self.xparam_frm, text="Web t (mm):")
			tw_lbl.grid(row=0, column=3)
			tw_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			tw_entry.grid(row=1, column=3)
			self.xparam_entries.append(tw_entry)
		if(region == "annulus"):
			ro_lbl = tk.Label(self.xparam_frm, text="Outer Radius (mm):")
			ro_lbl.grid(row=0, column=0)
			ro_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			ro_entry.grid(row=1, column=0)
			self.xparam_entries.append(ro_entry)
			ri_lbl = tk.Label(self.xparam_frm, text="Inner Radius (mm):")
			ri_lbl.grid(row=0, column=1)
			ri_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			ri_entry.grid(row=1, column=1)
			self.xparam_entries.append(ri_entry)
	
	#Return half of the height of the beam being added (in m)
	def half_h(self):
		xsec = self.xsec.get()
		return Region.half_h(xsec, self.get_xparams())
	
	#Return the parameters for the cross section
	def get_xparams(self):
		params = []
		for param in self.xparam_entries:
			params.append(float(param.get()))
		return params
	
	def get_L(self):
		return float(self.L_entry.get())
	
	def get_matl(self):
		return self.matl.get()
	
	def get_xsec(self):
		return self.xsec.get()
	
	def get_th(self):
		return float(self.th_entry.get())
	
	#Return true if all fields have numbers. Also returns false if any dim<=0
	def has_float_vals(self):
		try:
			if float(self.L_entry.get()) <= 0:
				return False
			for v in self.xparam_entries:
				if float(v.get()) <= 0:
					return False
			float(self.th_entry.get())
		except ValueError:
			#Flash fields red here?
			return False
		else:
			return True

#Add support toolbar
class Add_sup:
	add_sp_txt = "Add New\nSupport"
	add_jt_txt = "Add New\nJoint"
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		
		#Add Support or Joint
		self.sp_jt = tk.IntVar(self.tb_frm)
		self.sp_jt.set(0)
		sp_btn = tk.Radiobutton(self.tb_frm, variable=self.sp_jt, value=0, command=self.setsp)
		sp_btn.config(indicatoron=0, text="Support")
		sp_btn.grid(row=0, column=0, sticky=tk.W+tk.E)
		jt_btn = tk.Radiobutton(self.tb_frm, variable=self.sp_jt, value=1, command=self.setjt)
		jt_btn.config(indicatoron=0, text="Joint")
		jt_btn.grid(row=1, column=0, sticky=tk.W+tk.E)
		
		#First label
		self.tb_lbl = tk.Label(self.tb_frm, text=self.add_sp_txt)
		self.tb_lbl.grid(row=0, column=1, rowspan=2)
		
		next_col = 2
		#Radio buttons for support type
		self.sup_type = tk.IntVar(self.tb_frm)
		self.sup_type.set(0)
		for val, txt in Support.sup_types.items():
			s_btn = tk.Radiobutton(self.tb_frm, variable=self.sup_type, value=val)
			s_btn.config(indicatoron=0, text=txt, width=8)
			s_btn.grid(row=0, column=next_col, rowspan=2)
			next_col += 1
		
		#Angle label
		th_lbl = tk.Label(self.tb_frm, text="Angle (\u00B0):")
		th_lbl.grid(row=0, column=next_col)
		#Angle entry
		self.th_entry = tk.Entry(self.tb_frm)
		self.th_entry.config(width=num_e_wid)
		self.th_entry.grid(row=1, column=next_col)
		self.th_entry.insert(0, "0")
		next_col += 1
		#TO DO: Add option to snap to axis angle of member
		# - Angle: [ --Auto-- ] : darkened out - not accepting input
		
		#Button to add the new support
		self.add_btn = tk.Button(self.tb_frm, text="Add")
		#self.add_btn.config(command=self.toggle_add)
		self.add_btn.grid(row=0, column=next_col, rowspan=2, padx=2, pady=2, ipadx=8, sticky=tk.N+tk.S)
	
	def get_sup_type(self):
		return self.sup_type.get()
	
	def get_th(self):
		if self.th_entry["state"] == "disabled":
			return "[auto]"
		try: 
			th = float(self.th_entry.get())
		except ValueError: 
			return None
		return th
	
	def is_jt(self):
		return self.sp_jt.get()
	
	def setsp(self):
		self.tb_lbl.config(text=self.add_sp_txt)
	
	def setjt(self):
		self.tb_lbl.config(text=self.add_jt_txt)
	
	#Turn on or off auto angle mode
	#	onoff - boolean to say if it should be turned on or off
	def auto_th(self, onoff):
		if onoff:
			self.th_entry.delete(0, "end")
			self.th_entry.insert(0, "[auto]")
			self.th_entry.config(state="disabled")
		else:
			self.th_entry.config(state="normal")
			self.th_entry.delete(0, "end")
			self.th_entry.insert(0, "0")

#Add load toolbar
class Add_load:
	#Text for Point loads.
	xctext = "x-comp. (kN):"
	yctext = "y-comp. (kN):"
	rtext = "load (kN):"
	thtext = "angle (deg):"
	#Text for Distr. loads. Append q0 or q1 to the start for these four
	qxtext = "x-comp. (kN/m):" 
	qytext = "y-comp. (kN/m):"
	qrtext = "mag (kN/m):"
	qthtext = "angle (deg):"
	#Text for Moments.
	mxtext = "x-comp. (kN-m):"
	mytext = "y-comp. (kN-m):"
	mztext = "z-comp. (kN-m):"
	allow_xy_moments = False
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		
		#First label
		tb_lbl = tk.Label(self.tb_frm, text="Add New\nLoad")
		tb_lbl.grid(row=0, column=0, rowspan=2)
		
		#Choose load type (point, distributed, or moment) label
		ldtype_lbl = tk.Label(self.tb_frm, text="Load Type:")
		ldtype_lbl.grid(row=0, column=1)
		#Choose load type pulldown
		self.ldtype = tk.StringVar(self.tb_frm)
		self.ldtype.set(Load.load_types[0])
		self.ldtype.trace("w", self.set_ldtype) #Trace "w" --> run function on select
		ldtype_option = tk.OptionMenu(self.tb_frm, self.ldtype, *Load.load_types)
		ldtype_option.config(width=option_wid)
		ldtype_option.grid(row=1, column=1)
		
		#Choose components or r,theta
		self.c_p = tk.IntVar(self.tb_frm)
		self.c_p.set(0)
		self.c_btn = tk.Radiobutton(self.tb_frm, variable=self.c_p, value=0, command=self.setcomp)
		self.c_btn.config(indicatoron=0, text="Components")
		self.c_btn.grid(row=0, column=2, sticky=tk.W+tk.E)
		self.p_btn = tk.Radiobutton(self.tb_frm, variable=self.c_p, value=1, command=self.setpol)
		self.p_btn.config(indicatoron=0, text="Polar")
		self.p_btn.grid(row=1, column=2, sticky=tk.W+tk.E)
		
		#Frame that adjusts itself to have the needed fields
		self.comp_frm = tk.Frame(self.tb_frm)
		self.comp_frm.config(borderwidth=2, relief=tk.SUNKEN)
		self.comp_frm.grid(row=0, column=3, rowspan=2, sticky=tk.N+tk.S)
		self.comp_frm.grid_rowconfigure(1,weight=1)
		self.setpt()
		
		#Button to add the new load
		self.add_btn = tk.Button(self.tb_frm, text="Add")
		self.add_btn.grid(row=0, column=5, padx=2, pady=2, ipadx=8, rowspan=2, sticky=tk.N+tk.S)
	
	#Returns 0, 1, or 2, corrosponding to "Point", "Distributed", or "Moment"
	def get_ldtype(self):
		return Load.load_types.index(self.ldtype.get())
	#OLD
	def is_ds(self):
		print(308, "GET RID OF THIS")
		return self.get_ldtype() == 1
	def setcomp(self):
		if self.get_ldtype() == 0: #Pt. load comp. labels
			self.Pc1_lbl.config(text=self.xctext)
			self.Pc2_lbl.config(text=self.yctext)
		elif self.get_ldtype() == 1: #Distr. load comp. labels
			self.Pc1_lbl.config(text="q0 "+self.qxtext)
			self.Pc2_lbl.config(text="q0 "+self.qytext)
			self.Pc3_lbl.config(text="q1 "+self.qxtext)
			self.Pc4_lbl.config(text="q1 "+self.qytext)
	def setpol(self):
		if self.get_ldtype() == 0: #Pt. load polar labels
			self.Pc1_lbl.config(text=self.rtext)
			self.Pc2_lbl.config(text=self.thtext)
		if self.get_ldtype() == 1: #Distr. load polar labels
			self.Pc1_lbl.config(text="q0 "+self.rtext)
			self.Pc2_lbl.config(text="q0 "+self.thtext)
			self.Pc3_lbl.config(text="q1 "+self.qrtext)
			self.Pc4_lbl.config(text="q1 "+self.qthtext)
	#This function is called whenever Load Type is selected
	def set_ldtype(self, *args):
		ldtype = self.ldtype.get()
		if ldtype == Load.load_types[0]: #Point
			self.setpt()
		elif ldtype == Load.load_types[1]: #Distr
			self.setds()
		elif ldtype == Load.load_types[2]: #Moment
			self.setmt()
		else:
			#error
			print(331, "Error, unknown load type selected")
	
	def setpt(self):
		#Restore those radio buttons if they were hidden
		self.c_btn.grid()
		self.p_btn.grid()
		
		for widget in self.comp_frm.winfo_children():
			widget.destroy()
		#Load comp. 1 label
		self.Pc1_lbl = tk.Label(self.comp_frm, text=self.xctext)
		self.Pc1_lbl.grid(row=0, column=0)
		#Load comp. 1  entry
		self.Pc1_entry = tk.Entry(self.comp_frm)
		self.Pc1_entry.config(width=num_e_wid)
		self.Pc1_entry.grid(row=1, column=0)
		#Load comp. 2  label
		self.Pc2_lbl = tk.Label(self.comp_frm, text=self.yctext)
		self.Pc2_lbl.grid(row=0, column=1)
		#Load comp. 2  entry
		self.Pc2_entry = tk.Entry(self.comp_frm)
		self.Pc2_entry.config(width=num_e_wid)
		self.Pc2_entry.grid(row=1, column=1)
		#If in polar mode, fix the labels.
		if self.c_p.get():
			self.setpol()
	
	def setds(self):
		#Restore those radio buttons if they were hidden
		self.c_btn.grid()
		self.p_btn.grid()
		
		for widget in self.comp_frm.winfo_children():
			widget.destroy()
		#Q 0 comp. 1 label
		self.Pc1_lbl = tk.Label(self.comp_frm, text="q0 "+self.qxtext)
		self.Pc1_lbl.grid(row=0, column=0)
		#Q 0 comp. 1  entry
		self.Pc1_entry = tk.Entry(self.comp_frm)
		self.Pc1_entry.config(width=num_e_wid)
		self.Pc1_entry.grid(row=1, column=0)
		#Q 0 comp. 2  label
		self.Pc2_lbl = tk.Label(self.comp_frm, text="q0 "+self.qytext)
		self.Pc2_lbl.grid(row=0, column=1)
		#Q 0 comp. 2  entry
		self.Pc2_entry = tk.Entry(self.comp_frm)
		self.Pc2_entry.config(width=num_e_wid)
		self.Pc2_entry.grid(row=1, column=1)
		
		#Q 1 comp. 1 label
		self.Pc3_lbl = tk.Label(self.comp_frm, text="q1 "+self.qxtext)
		self.Pc3_lbl.grid(row=0, column=2)
		#Q 1 comp. 1  entry
		self.Pc3_entry = tk.Entry(self.comp_frm)
		self.Pc3_entry.config(width=num_e_wid)
		self.Pc3_entry.grid(row=1, column=2)
		#Q 1 comp. 2  label
		self.Pc4_lbl = tk.Label(self.comp_frm, text="q1 "+self.qytext)
		self.Pc4_lbl.grid(row=0, column=3)
		#Q 1 comp. 2  entry
		self.Pc4_entry = tk.Entry(self.comp_frm)
		self.Pc4_entry.config(width=num_e_wid)
		self.Pc4_entry.grid(row=1, column=3)
		#If in polar mode, fix the labels.
		if self.c_p.get():
			self.setpol()
	
	#Set up the interface for adding applied moments.
	def setmt(self):
		self.c_btn.grid_remove()
		self.p_btn.grid_remove()
		
		for widget in self.comp_frm.winfo_children():
			widget.destroy()
		
		#Turn on or off x and y moments
		if self.allow_xy_moments:
			#Load comp. 1 label
			self.Pc1_lbl = tk.Label(self.comp_frm, text=self.mxtext)
			self.Pc1_lbl.grid(row=0, column=0)
			#Load comp. 1  entry
			self.Pc1_entry = tk.Entry(self.comp_frm)
			self.Pc1_entry.config(width=num_e_wid)
			self.Pc1_entry.grid(row=1, column=0)
			#Load comp. 2  label
			self.Pc2_lbl = tk.Label(self.comp_frm, text=self.mytext)
			self.Pc2_lbl.grid(row=0, column=1)
			#Load comp. 2  entry
			self.Pc2_entry = tk.Entry(self.comp_frm)
			self.Pc2_entry.config(width=num_e_wid)
			self.Pc2_entry.grid(row=1, column=1)
		else:
			#Still have the tk.Entry-s, but not interactively
			self.Pc1_entry = tk.Entry(self.comp_frm)
			self.Pc1_entry.insert(tk.END, "0")
			self.Pc2_entry = tk.Entry(self.comp_frm)
			self.Pc2_entry.insert(tk.END, "0")
		#Load comp. 3  label
		self.Pc3_lbl = tk.Label(self.comp_frm, text=self.mztext)
		self.Pc3_lbl.grid(row=0, column=2)
		#Load comp. 3  entry
		self.Pc3_entry = tk.Entry(self.comp_frm)
		self.Pc3_entry.config(width=num_e_wid)
		self.Pc3_entry.grid(row=1, column=2)
	
	#Returns the components of the load in N
	def get_P(self):
		try:
			Pc1 = float(self.Pc1_entry.get())
			Pc2 = float(self.Pc2_entry.get())
		except ValueError:
			return ("NaN","NaN")
		ldt = self.get_ldtype()
		if ldt == 0:
			if self.c_p.get() == 0: #"Components"
				return (Pc1*1000, Pc2*1000)
			if self.c_p.get() == 1: #"Polar"
				return self.p_to_c(Pc1*1000, Pc2)
		elif ldt == 1:
			try:
				Pc3 = float(self.Pc3_entry.get())
				Pc4 = float(self.Pc4_entry.get())
			except ValueError:
				return ("NaN","NaN")
			if self.c_p.get() == 0:
				return ( (Pc1*1000, Pc2*1000), (Pc3*1000, Pc4*1000) )
			if self.c_p.get() == 1:
				return ( self.p_to_c(Pc1*1000, Pc2), self.p_to_c(Pc3*1000, Pc4))
		elif ldt == 2:
			try:
				Pc3 = float(self.Pc3_entry.get())
			except ValueError:
				return ("NaN","NaN", "NaN")
			return (Pc1*1000, Pc2*1000, Pc3*1000)
	
	def mag(self):
		ldt = self.get_ldtype()
		if ldt == 0:
			x,y = self.get_P()
			if self.c_p.get() == 0:
				return math.sqrt(x**2+y**2)
			if self.c_p.get() == 1:
				return x
		elif ldt == 1:
			(x0,y0), (x1,y1) = self.get_P()
			if self.c_p.get() == 0:
				return ( math.sqrt(x0**2+y0**2), math.sqrt(x1**2+y1**2) )
			if self.c_p.get() == 1:
				return ( x0, x1 )
		elif ldt == 2:
			x,y,z = self.get_P()
			return math.sqrt(x**2 + y**2 + z**2)
	#Polar to Coords. th in deg
	@staticmethod
	def p_to_c(r, th):
		x = r*math.cos(math.radians(th))
		y = r*math.sin(math.radians(th))
		return (x,y)
	#Return true if all fields have numbers. Also returns false if mag==0
	def has_float_vals(self):
		try:
			if self.mag() == 0:
				return False
			float(self.Pc1_entry.get())
			float(self.Pc2_entry.get())
			if self.get_ldtype() == 2:
				float(self.Pc3_entry.get())
			elif self.get_ldtype() == 1:
				float(self.Pc3_entry.get())
				float(self.Pc4_entry.get())
		except ValueError:
			return False
		#print(371)
		return True

