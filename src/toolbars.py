import tkinter as tk
import math

from member import *
from region import *
from support import Support


num_e_wid = 8

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
		matl_option.config(width=12)
		matl_option.grid(row=1, column=1)

		#Choose material label
		xsec_lbl = tk.Label(self.tb_frm, text="Cross Section Type:")
		xsec_lbl.grid(row=0, column=2)
		#Choose xsection pulldown
		self.xsec = tk.StringVar(self.tb_frm)
		self.xsec.set(Region.regions[0])
		self.xsec.trace("w", self.update_xparam)
		xsec_option = tk.OptionMenu(self.tb_frm, self.xsec, *Region.regions)
		xsec_option.config(width=12)
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

		#Radio buttons for vert/horiz
		self.vh = tk.IntVar(self.tb_frm)
		self.vh.set(1)
		v_btn = tk.Radiobutton(self.tb_frm, variable=self.vh, value=1)
		v_btn.config(indicatoron=0, text="Vertical")
		v_btn.grid(row=0, column=5, sticky=tk.W+tk.E)
		h_btn = tk.Radiobutton(self.tb_frm, variable=self.vh, value=0)
		h_btn.config(indicatoron=0, text="Horizontal")
		h_btn.grid(row=1, column=5, sticky=tk.W+tk.E)

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
		if xsec == "circle" or xsec == "annulus":
			try:
				r = float(self.get_xparams()[0])/1000 #mm
			except:
				return "NaN"
			return r
		if xsec == "rectangle":
			try:
				h = float(self.get_xparams()[1])/1000
			except:
				return "NaN"
			return h/2
		if xsec == "I-beam":
			try:
				d = float(self.get_xparams()[0])/1000
			except:
				return "NaN"
			return d/2
	def get_xparams(self):
		params = []
		for param in self.xparam_entries:
			params.append(param.get())
		return params
	def get_L(self):
		return self.L_entry.get()
	def get_matl(self):
		return self.matl.get()
	def get_xsec(self):
		return self.xsec.get()
	def get_vh(self):
		return self.vh.get()
	#Return true if all fields have numbers. Also returns false if any dim<=0
	def has_float_vals(self):
		try:
			if float(self.get_L()) <= 0:
				return false
			for v in self.get_xparams():
				if float(v) <= 0:
					return false
		except:#Flash fields red here?
			return False
		return True

#Add support toolbar
class Add_sup:
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		
		#First label
		tb_lbl = tk.Label(self.tb_frm, text="Add New\nSupport")
		tb_lbl.grid(row=0, column=0)
		
		next_col = 1
		#Radio buttons for support type
		self.sup_type = tk.IntVar(self.tb_frm)
		self.sup_type.set(0)
		for val, txt in Support.sup_types.items():
			s_btn = tk.Radiobutton(self.tb_frm, variable=self.sup_type, value=val)
			s_btn.config(indicatoron=0, text=txt, width=8)
			s_btn.grid(row=0, column=next_col)
			next_col += 1
		
		#Button to add the new support
		self.add_btn = tk.Button(self.tb_frm, text="Add")
		#self.add_btn.config(command=self.toggle_add)
		self.add_btn.grid(row=0, column=next_col, padx=2, pady=2, ipadx=8, sticky=tk.N+tk.S)
	
	def get_sup_type(self):
		return self.sup_type.get()

#Add load toolbar
class Add_load:
	xctext = "x-comp. (kN):"
	yctext = "y-comp. (kN):"
	rtext = "load (kN):"
	thtext = "angle (deg):"
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		
		#First label
		tb_lbl = tk.Label(self.tb_frm, text="Add New\nLoad")
		tb_lbl.grid(row=0, column=0, rowspan=2)
		
		#Choose components or r,theta
		self.c_p = tk.IntVar(self.tb_frm)
		self.c_p.set(0)
		c_btn = tk.Radiobutton(self.tb_frm, variable=self.c_p, value=0, command=self.setcomp)
		c_btn.config(indicatoron=0, text="Components")
		c_btn.grid(row=0, column=1, sticky=tk.W+tk.E)
		p_btn = tk.Radiobutton(self.tb_frm, variable=self.c_p, value=1, command=self.setpol)
		p_btn.config(indicatoron=0, text="Polar")
		p_btn.grid(row=1, column=1, sticky=tk.W+tk.E)
		
		#Load comp. 1 label
		self.Pc1_lbl = tk.Label(self.tb_frm, text=self.xctext)
		self.Pc1_lbl.grid(row=0, column=2)
		#Load comp. 1  entry
		self.Pc1_entry = tk.Entry(self.tb_frm)
		self.Pc1_entry.config(width=num_e_wid)
		self.Pc1_entry.grid(row=1, column=2)
		
		#Load comp. 2  label
		self.Pc2_lbl = tk.Label(self.tb_frm, text=self.yctext)
		self.Pc2_lbl.grid(row=0, column=3)
		#Load comp. 2  entry
		self.Pc2_entry = tk.Entry(self.tb_frm)
		self.Pc2_entry.config(width=num_e_wid)
		self.Pc2_entry.grid(row=1, column=3)
		
		#Add option for distributed loads
		
		#Button to add the new load
		self.add_btn = tk.Button(self.tb_frm, text="Add")
		self.add_btn.grid(row=0, column=4, padx=2, pady=2, ipadx=8, rowspan=2, sticky=tk.N+tk.S)
		
		
	def setcomp(self):
		self.Pc1_lbl.config(text=self.xctext)
		self.Pc2_lbl.config(text=self.yctext)
	def setpol(self):
		self.Pc1_lbl.config(text=self.rtext)
		self.Pc2_lbl.config(text=self.thtext)
	
	#Returns the components of the load
	def get_P(self):
		try:
			Pc1 = float(self.Pc1_entry.get())
			Pc2 = float(self.Pc2_entry.get())
		except:
			return ("NaN","NaN")
		if self.c_p.get() == 0: #"Components"
			return (Pc1, Pc2)
		if self.c_p.get() == 1: #"Polar"
			Px = Pc1*math.cos(math.radians(Pc2))
			Py = Pc1*math.sin(math.radians(Pc2))
			return (Px, Py)
	def mag(self):
		x,y = self.get_P()
		if self.c_p.get() == 0:
			return math.sqrt(x**2+y**2)
		if self.c_p.get() == 1:
			return x
	#Return true if all fields have numbers. Also returns false if mag==0
	def has_float_vals(self):
		try:
			if self.mag() == 0:
				return False
			float(self.Pc1_entry.get())
			float(self.Pc2_entry.get())
		except:
			return False
		return True

