import tkinter as tk

from lab import Lab
from toolbars import Add_mem, Add_sup


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

#Lab with canvas
mlab = Lab(main_frm, add_mem_bar, add_sup_bar)

#Start
tk.mainloop()
