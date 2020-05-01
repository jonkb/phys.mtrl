import tkinter as tk

from member import *
from region import *

num_e_wid = 8

class Add_mem:
	def __init__(self, main_frm):
		self.add_mode = False
		
		#Toolbar to add a new member
		self.add_mem_frm = tk.Frame(main_frm)
		self.add_mem_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.add_mem_frm.pack(side=tk.TOP, fill=tk.X)
		self.add_mem_frm.grid_rowconfigure(1,weight=1)

		#Choose material label
		self.matl_lbl = tk.Label(self.add_mem_frm, text="Material:")
		self.matl_lbl.grid(row=0, column=0)
		#Choose material pulldown
		self.add_matl = tk.StringVar(self.add_mem_frm)
		self.add_matl.set(Materials.materials[0])
		self.add_matl_option = tk.OptionMenu(self.add_mem_frm, self.add_matl, *Materials.materials)
		self.add_matl_option.config(width=12)
		self.add_matl_option.grid(row=1, column=0)

		#Choose material label
		self.xsec_lbl = tk.Label(self.add_mem_frm, text="Cross Section Type:")
		self.xsec_lbl.grid(row=0, column=1)
		#Choose xsection pulldown
		self.add_xsec = tk.StringVar(self.add_mem_frm)
		self.add_xsec.set(Region.regions[0])
		self.add_xsec.trace("w", self.update_xparam)
		self.add_xsec_option = tk.OptionMenu(self.add_mem_frm, self.add_xsec, *Region.regions)
		self.add_xsec_option.config(width=12)
		self.add_xsec_option.grid(row=1, column=1)

		#Frame that adjusts itself to the chosen xsection to have the needed parameters
		self.xparam_frm = tk.Frame(self.add_mem_frm)
		self.xparam_frm.config(borderwidth=2, relief=tk.SUNKEN)
		self.xparam_frm.grid(row=0, column=2, rowspan=2, sticky=tk.N+tk.S)
		self.xparam_frm.grid_rowconfigure(1,weight=1)
		self.xparam_entries = []
		self.update_xparam()

		#Length label
		self.L_lbl = tk.Label(self.add_mem_frm, text="Length (m):")
		self.L_lbl.grid(row=0, column=3)
		#Length entry
		self.add_L_entry = tk.Entry(self.add_mem_frm)
		self.add_L_entry.config(width=num_e_wid)
		self.add_L_entry.grid(row=1, column=3)

		#Radio buttons for vert/horiz
		self.add_vh = tk.IntVar(self.add_mem_frm)
		self.add_vh.set(1)
		self.add_v_btn = tk.Radiobutton(self.add_mem_frm, variable=self.add_vh, value=1)
		self.add_v_btn.config(indicatoron=0, text="Vertical")
		self.add_v_btn.grid(row=0, column=4, sticky=tk.W+tk.E)
		self.add_h_btn = tk.Radiobutton(self.add_mem_frm, variable=self.add_vh, value=0)
		self.add_h_btn.config(indicatoron=0, text="Horizontal")
		self.add_h_btn.grid(row=1, column=4, sticky=tk.W+tk.E)

		#Button to add the new member
		self.add_btn = tk.Button(self.add_mem_frm, text="Add")
		self.add_btn.config(command=self.toggle_add)
		self.add_btn.grid(row=0, column=5, padx=2, pady=2, ipadx=8, rowspan=2, sticky=tk.N+tk.S)

	def update_xparam(self, *args):
		region = self.add_xsec.get()
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
			d_lbl.grid(row=0, column=1)
			d_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			d_entry.grid(row=1, column=1)
			self.xparam_entries.append(d_entry)
			w_lbl = tk.Label(self.xparam_frm, text="Width (mm):")
			w_lbl.grid(row=0, column=0)
			w_entry = tk.Entry(self.xparam_frm, width=num_e_wid)
			w_entry.grid(row=1, column=0)
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
	def toggle_add(self):
		try:
			L = float(self.add_L_entry.get())
		except:#Flash the L field red here?
			return
		self.add_mode = not self.add_mode
	#Return half of the height of the beam being added (in m)
	def half_h(self):
		xsec = self.add_xsec.get()
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
