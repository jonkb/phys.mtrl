import tkinter as tk

from lab import Lab
from member import Member
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
filemenu.add_command(label="Open Lab", command=mlab.open)
filemenu.add_command(label="Save Lab", command=mlab.save)
menubar.add_cascade(label="File", menu=filemenu)
viewmenu = tk.Menu(menubar, tearoff=0)
viewmenu.add_command(label="Show All Toolbars", command=mlab.show_allt)
viewmenu.add_command(label="Hide All Toolbars", command=mlab.hide_allt)
show_mem = tk.IntVar()
show_mem.set(True)
mlab.show_mem = show_mem
viewmenu.add_checkbutton(label="Add Member Toolbar", variable=show_mem, command=mlab.toggle_smem)
show_sup = tk.IntVar()
show_sup.set(True)
mlab.show_sup = show_sup
viewmenu.add_checkbutton(label="Add Support Toolbar", variable=show_sup, command=mlab.toggle_ssup)
show_ld = tk.IntVar()
show_ld.set(True)
mlab.show_ld = show_ld
viewmenu.add_checkbutton(label="Add Load Toolbar", variable=show_ld, command=mlab.toggle_sld)
menubar.add_cascade(label="View", menu=viewmenu)
optmenu = tk.Menu(menubar, tearoff=0)
mem_wtls = tk.IntVar()
mem_wtls.set(True)
mlab.mem_wtls = mem_wtls
optmenu.add_checkbutton(label="Weightless Members", variable=mem_wtls, command=mlab.toggle_wtls)
optmenu.add_command(label="Set Scale", command=mlab.edit_scale)
menubar.add_cascade(label="Options", menu=optmenu)
evalmenu = tk.Menu(menubar, tearoff=0)
for i, evnm in Member.eval_names.items():
	evalmenu.add_command(label=evnm, command=lambda ri=i: mlab.eval_report(ri))
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
