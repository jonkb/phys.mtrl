import tkinter as tk

from lab import Lab
from toolbars import Add_mem, Add_sup, Add_load
from tk_wig import Tk_rt, PM_Menu


#Root window
root = Tk_rt("Materials Physics")

#Main frame: container for everything else
main_frm = tk.Frame(root)
main_frm.pack(fill=tk.BOTH, expand=1)

#Toolbars
add_mem_bar = Add_mem(main_frm)
add_sup_bar = Add_sup(main_frm)
add_load_bar = Add_load(main_frm)

#"Materials" Lab with canvas
mlab = Lab(main_frm, add_mem_bar, add_sup_bar, add_load_bar)

#Menu bar
menubar = PM_Menu(root, mlab)
root.config(menu=menubar)

def cleanup():
	mlab.cleanup()
	root.destroy()
root.protocol("WM_DELETE_WINDOW", cleanup)

#Start
root.mainloop()
