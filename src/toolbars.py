import tkinter as tk
from tkinter import font
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
		return Region.half_h(xsec, self.get_xparams())
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
		except ValueError:#Flash fields red here?
			return False
		else:
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
	qxtext = "x-comp. (kN/m)" # Append q0 or q1 to the start for these four
	qytext = "y-comp. (kN/m)"
	qrtext = "mag (kN/m)"
	qthtext = "angle (deg)"
	def __init__(self, main_frm):
		self.tb_frm = tk.Frame(main_frm)
		self.tb_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
		self.tb_frm.pack(side=tk.TOP, fill=tk.X)
		
		#First label
		tb_lbl = tk.Label(self.tb_frm, text="Add New\nLoad")
		tb_lbl.grid(row=0, column=0, rowspan=2)
		
		#Choose point or distributed
		self.pt_ds = tk.IntVar(self.tb_frm)
		self.pt_ds.set(0)
		pt_btn = tk.Radiobutton(self.tb_frm, variable=self.pt_ds, value=0, command=self.setpt)
		pt_btn.config(indicatoron=0, text="Point")
		pt_btn.grid(row=0, column=1, sticky=tk.W+tk.E)
		ds_btn = tk.Radiobutton(self.tb_frm, variable=self.pt_ds, value=1, command=self.setds)
		ds_btn.config(indicatoron=0, text="Distributed")
		ds_btn.grid(row=1, column=1, sticky=tk.W+tk.E)
		
		#Choose components or r,theta
		self.c_p = tk.IntVar(self.tb_frm)
		self.c_p.set(0)
		c_btn = tk.Radiobutton(self.tb_frm, variable=self.c_p, value=0, command=self.setcomp)
		c_btn.config(indicatoron=0, text="Components")
		c_btn.grid(row=0, column=2, sticky=tk.W+tk.E)
		p_btn = tk.Radiobutton(self.tb_frm, variable=self.c_p, value=1, command=self.setpol)
		p_btn.config(indicatoron=0, text="Polar")
		p_btn.grid(row=1, column=2, sticky=tk.W+tk.E)
		
		#Frame that adjusts itself to have the needed fields
		self.comp_frm = tk.Frame(self.tb_frm)
		self.comp_frm.config(borderwidth=2, relief=tk.SUNKEN)
		self.comp_frm.grid(row=0, column=3, rowspan=2, sticky=tk.N+tk.S)
		self.comp_frm.grid_rowconfigure(1,weight=1)
		self.setpt()
		
		#Button to add the new load
		self.add_btn = tk.Button(self.tb_frm, text="Add")
		self.add_btn.grid(row=0, column=5, padx=2, pady=2, ipadx=8, rowspan=2, sticky=tk.N+tk.S)
		
	def is_ds(self):
		return self.pt_ds.get()
	def setcomp(self):
		self.Pc1_lbl.config(text=self.xctext)
		self.Pc2_lbl.config(text=self.yctext)
		if self.pt_ds.get() == 1:
			self.Pc3_lbl.config(text="q0 "+self.qxtext)
			self.Pc4_lbl.config(text="q0 "+self.qytext)
	def setpol(self):
		self.Pc1_lbl.config(text=self.rtext)
		self.Pc2_lbl.config(text=self.thtext)
		if self.pt_ds.get() == 1:
			self.Pc3_lbl.config(text="q1 "+self.qrtext)
			self.Pc4_lbl.config(text="q1 "+self.qthtext)
	def setpt(self):
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
	def setds(self):
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
	
	#Returns the components of the load in N
	def get_P(self):
		try:
			Pc1 = float(self.Pc1_entry.get())
			Pc2 = float(self.Pc2_entry.get())
		except ValueError:
			return ("NaN","NaN")
		if self.pt_ds.get() == 0:
			if self.c_p.get() == 0: #"Components"
				return (Pc1*1000, Pc2*1000)
			if self.c_p.get() == 1: #"Polar"
				return self.p_to_c(Pc1*1000, Pc2)
		elif self.pt_ds.get() == 1:
			try:
				Pc3 = float(self.Pc3_entry.get())
				Pc4 = float(self.Pc4_entry.get())
			except ValueError:
				return ("NaN","NaN")
			if self.c_p.get() == 0:
				return ( (Pc1*1000, Pc2*1000), (Pc3*1000, Pc4*1000) )
			if self.c_p.get() == 1:
				return ( self.p_to_c(Pc1*1000, Pc2), self.p_to_c(Pc3*1000, Pc4))
	def mag(self):
		if self.pt_ds.get() == 0:
			x,y = self.get_P()
			if self.c_p.get() == 0:
				return math.sqrt(x**2+y**2)
			if self.c_p.get() == 1:
				return x
		if self.pt_ds.get() == 1:
			(x0,y0), (x1,y1) = self.get_P()
			if self.c_p.get() == 0:
				return ( math.sqrt(x0**2+y0**2), math.sqrt(x1**2+y1**2) )
			if self.c_p.get() == 1:
				return ( x0, x1 )
	#Polar to Coords. th in deg
	@staticmethod
	def p_to_c(r, th):
		x = r*math.cos(math.radians(th))
		y = r*math.sin(math.radians(th))
		return (x,y)
	#Return true if all fields have numbers. Also returns false if mag==0
	def has_float_vals(self):
		#print(361)
		try:
			if self.mag() == 0:
				return False
			float(self.Pc1_entry.get())
			float(self.Pc2_entry.get())
			if self.pt_ds.get() == 1:
				float(self.Pc3_entry.get())
				float(self.Pc4_entry.get())
		except ValueError:
			return False
		#print(371)
		return True


#Expansion of tk.Text, used for displaying text, not receiving any input
#You can't copy text from a Label, which is a problem
class Txt_wig(tk.Text):
	ft_fam = "Helvetica"#"MS Sans Serif"#"System"
	ft_size = 11
	padw = 4
	#rep_lbl = tk.Label(root, text=report, justify=tk.LEFT)
	#tk.Text allows for copy / paste
	def __init__(self, root, txt):
		numlines = txt.count("\n") + 1
		#TO DO: One too few when there's a blank line ?????
		#print(374, numlines)
		super().__init__(root, relief="flat", height=numlines)
		
		#print(555, txt)
		self.txt = txt
		self.insert(tk.END, txt)
		#int(self.index("end-1c").split('.')[0])
		self.TWfont = font.Font(root=root, family=self.ft_fam, size=self.ft_size)
		print(self.TWfont.actual())
		#self.TWfont = font.Font(family="System", size=12)
		#print(font.families(root=root))
		print(self.TWfont)
		self.config(font=self.TWfont)
		#self.config(font=("System", 12))
		self.config(state="disabled")
		self.bind("<Configure>", self.resize)
	def packslf(self):
		self.pack(fill="both", padx=1, pady=1)
	def resize(self, event):
		#self.config(font=self.TWfont)
		#print(381, event.width, event.height)
		textw = event.width - self.padw
		#print(382, len(self.get("1.0", "end-1c").split("\n")))
		#print(390, self.TWfont.measure("abcdefghi"))
		#print(390.5, self.TWfont.measure("abc def ghi"))
		#print(391, self.TWfont.measure("abcde\tfghi"))
		#print(391, self.TWfont.measure("\tabcdefghi"))
		#print(391, self.TWfont.measure("\n\tabcdefghi"))
		#print(382, self.TWfont.measure(self.txt.split('\n')[0]))
		height = 0
		for line in self.txt.split('\n'):
			#Using max(1, l) means that empty lines still add to the height
			height += max(1, math.ceil(self.TWfont.measure(line) / textw))
		self.config(height=height)
	
def txtwigtest():
	root = tk.Tk()
	tw = Txt_wig(root, "ABCDEFGHIJKLMNOPQRSTUVWXYZ\nDEF\nGHIJKLM")
	tw.packslf()
	root.mainloop()
#txtwigtest()