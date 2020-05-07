import tkinter as tk

from lab import Lab
from toolbars import Add_mem, Add_sup, Add_load


#Root window
root = tk.Tk()
root.title("Materials Physics")
root.iconbitmap("../img/phys.ico")


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
menubar = tk.Menu(root)
filemenu = tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Nothing Here Yet")
menubar.add_cascade(label="File", menu=filemenu)
optmenu = tk.Menu(menubar, tearoff=0)
mem_wtls = tk.IntVar()
mem_wtls.set(True)
optmenu.add_checkbutton(label="Weightless Members", variable=mem_wtls)
optmenu.add_command(label="Set Scale", command=mlab.edit_scale)
menubar.add_cascade(label="Options", menu=optmenu)
mlab.mem_wtls = mem_wtls
evalmenu = tk.Menu(menubar, tearoff=0)
evalmenu.add_command(label="Axial Stress", command=mlab.eval_axial)
evalmenu.add_command(label="Euler Buckling", command=mlab.eval_buckling)
menubar.add_cascade(label="Evaluate", menu=evalmenu)
delmenu = tk.Menu(menubar, tearoff=0)
delmenu.add_command(label="Clear All", command=mlab.clear_all)
delmenu.add_command(label="Member", command=mlab.del_mode)
#evelmenu.add_command(label="Support", command=mlab.)
#evelmenu.add_command(label="Load", command=mlab.)
menubar.add_cascade(label="Delete", menu=delmenu)
root.config(menu=menubar)

def cleanup():
	mlab.cleanup()
	root.destroy()
root.protocol("WM_DELETE_WINDOW", cleanup)

#Start
root.mainloop()
