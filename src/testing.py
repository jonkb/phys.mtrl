import tkinter as tk
import load as ld
import lab
import toolbars

r = tk.Tk()
m = ld.Moment(123, 0, 0, 1, 0)
m1 = ld.Moment(123, 0, 0, -1, 0)
m2 = ld.Moment(123, 0, 0, 2000, 0)
m3 = ld.Moment(123, 5000, -5000, 2000, 0)

#lab object here
amt = toolbars.Add_mem(r)
ast = toolbars.Add_sup(r)
alt = toolbars.Add_load(r)
lb = lab.Lab(r, amt, ast, alt)
lb.px_per_m = 360
lb.px_per_kN = 20.0

m.draw(lb, 50, 50)
m1.draw(lb, 150, 50)
m2.draw(lb, 100, 150)
m3.draw(lb, 200, 150)

r.mainloop()
