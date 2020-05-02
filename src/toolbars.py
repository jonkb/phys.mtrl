import tkinter as tk

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
	#Return half of the height of the beam being added (in m)
	def half_h(self):
		xsec = self.xsec.get()
		if xsec == "circle":
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

#Add support toolbar
class Add_sup:
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		
		#First label
		tb_lbl = tk.Label(self.tb_frm, text="Add New\nSupport")
		tb_lbl.grid(row=0, column=0, rowspan=2)
		
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



#END



