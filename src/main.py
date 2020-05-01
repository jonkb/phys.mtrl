import tkinter as tk

from lab import Lab
from toolbars import Add_mem

#Root window
root = tk.Tk()
root.title("Materials Physics")
root.iconbitmap("../img/phys.ico")

main_frm = tk.Frame(root)
main_frm.pack(fill=tk.BOTH, expand=1)

add_mem_bar = Add_mem(main_frm)
mlab = Lab(main_frm, add_mem_bar)

tk.mainloop()
