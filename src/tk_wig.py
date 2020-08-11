import math
import tkinter as tk
from tkinter import font
import os

from member import Member

#Wrapper for tk.Tk() incorporating the title and icon already
class Tk_rt(tk.Tk):
	ico_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../img/phys.ico")
	def __init__(self, title):
		super().__init__()
		self.title(title)
		try:
			self.iconbitmap(self.ico_path)
		except:
			print("Error loading icon ("+ico_path+")")

#Phys.Mtrl Menu bar
#Honestly, this isn't a very useful class, since it's not reused.
#	It's just moving code from main.py to here.
class PM_Menu(tk.Menu):
	def __init__(self, root, mlab):
		super().__init__(root)
		#File Menu
		filemenu = tk.Menu(self, tearoff=0)
		filemenu.add_command(label="Open Lab", command=mlab.open_lab)
		filemenu.add_command(label="Save Lab", command=mlab.save_lab)
		self.add_cascade(label="File", menu=filemenu)
		#View Menu
		viewmenu = tk.Menu(self, tearoff=0)
		viewmenu.add_command(label="Set Scale", command=mlab.edit_scale)
		viewmenu.add_command(label="Show All Toolbars", command=mlab.show_allt)
		viewmenu.add_command(label="Hide All Toolbars", command=mlab.hide_allt)
		show_mem = tk.IntVar()
		show_mem.set(True)
		mlab.options.show_mem = show_mem
		viewmenu.add_checkbutton(label="Add Member Toolbar", 
			variable=show_mem, command=mlab.toggle_smem)
		show_sup = tk.IntVar()
		show_sup.set(True)
		mlab.options.show_sup = show_sup
		viewmenu.add_checkbutton(label="Add Support Toolbar", 
			variable=show_sup, command=mlab.toggle_ssup)
		show_ld = tk.IntVar()
		show_ld.set(True)
		mlab.options.show_ld = show_ld
		viewmenu.add_checkbutton(label="Add Load Toolbar", 
			variable=show_ld, command=mlab.toggle_sld)
		self.add_cascade(label="View", menu=viewmenu)
		#Options Menu
		optmenu = tk.Menu(self, tearoff=0)
		mem_wtls = tk.IntVar()
		mem_wtls.set(True)
		mlab.options.mem_wtls = mem_wtls
		optmenu.add_checkbutton(label="Weightless Members", 
			variable=mem_wtls, command=mlab.toggle_wtls)
		sup_endsnap = tk.IntVar()
		sup_endsnap.set(True)
		mlab.options.sup_endsnap = sup_endsnap
		optmenu.add_checkbutton(label="Snap Supports to Member Ends", 
			variable=sup_endsnap)
		self.add_cascade(label="Options", menu=optmenu)
		sup_thsnap = tk.IntVar()
		sup_thsnap.set(True)
		mlab.options.sup_thsnap = sup_thsnap
		mlab.toggle_thsnap() #Set it to True, to start out with
		optmenu.add_checkbutton(label="Auto Support Angle", 
			variable=sup_thsnap, command=mlab.toggle_thsnap)
		#Evaluation Menu
		evalmenu = tk.Menu(self, tearoff=0)
		for i, evnm in Member.eval_names.items():
			evalmenu.add_command(label=evnm, command=lambda ri=i: mlab.eval_report(ri))
		self.add_cascade(label="Evaluate", menu=evalmenu)
		delmenu = tk.Menu(self, tearoff=0)
		delmenu.add_command(label="Clear All", command=mlab.clear_all)
		delmenu.add_command(label="Member", command=mlab.del_mem_mode)
		delmenu.add_command(label="Ld/Sup/Jt", command=mlab.del_lsj_mode)
		self.add_cascade(label="Delete", menu=delmenu)

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
		#print(self.TWfont.actual())
		#self.TWfont = font.Font(family="System", size=12)
		#print(font.families(root=root))
		#print(self.TWfont)
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