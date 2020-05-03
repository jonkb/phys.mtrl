import tkinter as tk

from lab import Lab
from toolbars import Add_mem, Add_sup, Add_load


#Root window
root = tk.Tk()
root.title("Materials Physics")
root.iconbitmap("../img/phys.ico")

#Menu bar
menubar = tk.Menu(root)
filemenu = tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Nothing Here Yet")
menubar.add_cascade(label="File", menu=filemenu)
optmenu = tk.Menu(menubar, tearoff=0)
optmenu.add_command(label="Set Scale")
mem_wtls = tk.IntVar()
optmenu.add_checkbutton(label="Weightless Members", variable=mem_wtls)
menubar.add_cascade(label="Options", menu=optmenu)
root.config(menu=menubar)

#Main frame: container for everything else
main_frm = tk.Frame(root)
main_frm.pack(fill=tk.BOTH, expand=1)

#Toolbars
add_mem_bar = Add_mem(main_frm)
add_sup_bar = Add_sup(main_frm)
add_load_bar = Add_load(main_frm)

#"Materials" Lab with canvas
mlab = Lab(main_frm, add_mem_bar, add_sup_bar, add_load_bar)
mlab.mem_wtls = mem_wtls

#Start
tk.mainloop()
