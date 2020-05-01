import tkinter as tk
#My files
from member import *
from region import *

#Constants
px_per_m = 200#scale
dash_len = int(px_per_m/6)
num_e_wid = 8

#global vars...
c_wd = 800
c_ht = 500
#Remember that when c_wd = 800 pixels, this corrosponds to x=[0,799]
add_mode = False
floating_mem = None
add_vertical = True #If False, it's horizontal

#w = Member(Materials.steel, W_F_I(30.9, 15.1, 1.32, .775), 1)
#w = Member(Materials.steel, W_F_I(25, 13, 1.22, .705), 1)
#print(w)
#print(w.xarea)
#print(w.xIx)
#print(w.xIy)

def redraw_bg(event):
	global c_wd, c_ht
	c_wd, c_ht = event.width, event.height
	canv.coords(bgbox, 0, 0, c_wd-1,c_ht-1)
	#? Add something to move everything else, so it's all anchored to the bottom left
	canv.delete(*bggrid)
	#Make 1.0m^2 grid
	for linex in range(px_per_m, c_wd, px_per_m):
		bggrid.append(canv.create_line(linex, 0, linex, c_ht-1, fill="gray", dash=(dash_len,)))
	for liney in range(px_per_m, c_ht, px_per_m):
		bggrid.append(canv.create_line(0, c_ht-1-liney, c_wd-1, c_ht-1-liney, fill="gray", dash=(dash_len,)))
def update_xparam(*args):
	region = add_xsec.get()
	xparam_entries.clear()
	for widget in xparam_frm.winfo_children():
		widget.destroy()
	if(region == "circle"):
		rad_lbl = tk.Label(xparam_frm, text="Radius (mm):")
		rad_lbl.grid(row=0, column=0)
		rad_entry = tk.Entry(xparam_frm, width=num_e_wid)
		rad_entry.grid(row=1, column=0)
		xparam_entries.append(rad_entry)
	if(region == "rectangle"):
		b_lbl = tk.Label(xparam_frm, text="Base (mm):")
		b_lbl.grid(row=0, column=0)
		b_entry = tk.Entry(xparam_frm, width=num_e_wid)
		b_entry.grid(row=1, column=0)
		xparam_entries.append(b_entry)
		h_lbl = tk.Label(xparam_frm, text="Height (mm):")
		h_lbl.grid(row=0, column=1)
		h_entry = tk.Entry(xparam_frm, width=num_e_wid)
		h_entry.grid(row=1, column=1)
		xparam_entries.append(h_entry)
	if(region == "I-beam"):
		d_lbl = tk.Label(xparam_frm, text="Depth (mm):")
		d_lbl.grid(row=0, column=1)
		d_entry = tk.Entry(xparam_frm, width=num_e_wid)
		d_entry.grid(row=1, column=1)
		xparam_entries.append(d_entry)
		w_lbl = tk.Label(xparam_frm, text="Width (mm):")
		w_lbl.grid(row=0, column=0)
		w_entry = tk.Entry(xparam_frm, width=num_e_wid)
		w_entry.grid(row=1, column=0)
		xparam_entries.append(w_entry)
		tf_lbl = tk.Label(xparam_frm, text="Flange t (mm):")
		tf_lbl.grid(row=0, column=2)
		tf_entry = tk.Entry(xparam_frm, width=num_e_wid)
		tf_entry.grid(row=1, column=2)
		xparam_entries.append(tf_entry)
		tw_lbl = tk.Label(xparam_frm, text="Web t (mm):")
		tw_lbl.grid(row=0, column=3)
		tw_entry = tk.Entry(xparam_frm, width=num_e_wid)
		tw_entry.grid(row=1, column=3)
		xparam_entries.append(tw_entry)
def mouse_moved(event):
	if not add_mode:
		return
	global floating_mem, floating_x, floating_y
	x = event.x
	y = event.y
	if floating_mem == None:
		floating_x = x
		floating_y = y
		L = float(add_L_entry.get())
		floating_mem = canv.create_rectangle(*rect_coords(x, y, L*px_per_m))#(x, y-half_h(), x+L*px_per_m, y+half_h())
	else:
		canv.move(floating_mem, x-floating_x, y-floating_y)
		floating_x = x
		floating_y = y
def mouse_click(event):
	if add_mode:
		add_member(event)
	#print(canv.find_all())
def add_member(event):
	global floating_mem
	x = event.x
	y = event.y
	try:
		L = float(add_L_entry.get())
	except:#Flash the L field red here?
		return
	matl = getattr(Materials,add_matl.get())
	rgb = matl["color"]
	if add_xsec.get() == "circle":
		try:
			r = float(xparam_entries[0].get())/1000 #mm->m
		except:
			return
		xsec = Circle(r)
	if add_xsec.get() == "rectangle":
		try:
			b = float(xparam_entries[0].get())/1000 #mm->m
			h = float(xparam_entries[1].get())/1000 
		except:
			return
		xsec = Rectangle(b, h)
	if add_xsec.get() == "I-beam":
		try:
			d = float(xparam_entries[0].get())/1000 #mm->m
			w = float(xparam_entries[1].get())/1000 
			tf = float(xparam_entries[2].get())/1000 
			tw = float(xparam_entries[3].get())/1000 
		except:
			return
		xsec = W_F_I(d, w, tf, tw)
	m = Member(matl, xsec, 1)#Store this somewhere
	print("Added new " + str(m))
	canv.create_rectangle(*rect_coords(x, y, L*px_per_m), fill=rgb)
	toggle_add()
	canv.delete(floating_mem)
	floating_mem = None
def toggle_add():
	try:
		L = float(add_L_entry.get())
	except:#Flash the L field red here?
		return
	global add_mode
	add_mode = not add_mode
def half_h():
	if add_xsec.get() == "circle":
		try:
			r = float(xparam_entries[0].get())/1000 #mm
		except:
			return "NaN"
		return r*px_per_m
	if add_xsec.get() == "rectangle":
		try:
			h = float(xparam_entries[1].get())/1000
		except:
			return "NaN"
		return h/2*px_per_m
	if add_xsec.get() == "I-beam":
		try:
			d = float(xparam_entries[0].get())/1000
		except:
			return "NaN"
		return d/2*px_per_m
def rect_coords(x,y,L_px):
	hh = half_h()
	if hh=="NaN":
		return "NaN"
	if(add_vh.get()):
		x1 = x - hh
		y1 = y
		x2 = x + hh
		y2 = y - L_px
	else:
		x1 = x
		y1 = y - hh
		x2 = x + L_px
		y2 = y + hh
	return (x1, y1, x2, y2)

#Root window
root = tk.Tk()
root.title("Materials Physics")
root.iconbitmap("../img/phys.ico")

main_frm = tk.Frame(root)
main_frm.pack(fill=tk.BOTH, expand=1)
canv = tk.Canvas(main_frm, width=c_wd, height=c_ht, borderwidth=0, highlightthickness=0)

#Toolbar to add a new member
add_mem_frm = tk.Frame(main_frm)
add_mem_frm.config(highlightcolor="grey", highlightbackground="grey", highlightthickness=1)
add_mem_frm.pack(side=tk.TOP, fill=tk.X)
add_mem_frm.grid_rowconfigure(1,weight=1)

#Choose material label
matl_lbl = tk.Label(add_mem_frm, text="Material:")
matl_lbl.grid(row=0, column=0)
#Choose material pulldown
add_matl = tk.StringVar(add_mem_frm)
add_matl.set(Materials.materials[0])
add_matl_option = tk.OptionMenu(add_mem_frm, add_matl, *Materials.materials)
add_matl_option.config(width=12)
add_matl_option.grid(row=1, column=0)

#Choose material label
xsec_lbl = tk.Label(add_mem_frm, text="Cross Section Type:")
xsec_lbl.grid(row=0, column=1)
#Choose xsection pulldown
add_xsec = tk.StringVar(add_mem_frm)
add_xsec.set(Region.regions[0])
add_xsec.trace("w", update_xparam)
add_xsec_option = tk.OptionMenu(add_mem_frm, add_xsec, *Region.regions)
add_xsec_option.config(width=12)
add_xsec_option.grid(row=1, column=1)

#Frame that adjusts itself to the chosen xsection to have the needed parameters
xparam_frm = tk.Frame(add_mem_frm)
xparam_frm.config(borderwidth=2, relief=tk.SUNKEN)
xparam_frm.grid(row=0, column=2, rowspan=2, sticky=tk.N+tk.S)
xparam_frm.grid_rowconfigure(1,weight=1)
xparam_entries = []
update_xparam()

#Length label
L_lbl = tk.Label(add_mem_frm, text="Length (m):")
L_lbl.grid(row=0, column=3)
#Length entry
add_L_entry = tk.Entry(add_mem_frm)
add_L_entry.config(width=num_e_wid)
add_L_entry.grid(row=1, column=3)

#Radio buttons for vert/horiz
add_vh = tk.IntVar(add_mem_frm)
add_vh.set(1)
add_v_btn = tk.Radiobutton(add_mem_frm, variable=add_vh, value=1)
add_v_btn.config(indicatoron=0, text="Vertical")
add_v_btn.grid(row=0, column=4, sticky=tk.W+tk.E)
add_h_btn = tk.Radiobutton(add_mem_frm, variable=add_vh, value=0)
add_h_btn.config(indicatoron=0, text="Horizontal")
add_h_btn.grid(row=1, column=4, sticky=tk.W+tk.E)

#Button to add the new member
add_btn = tk.Button(add_mem_frm, text="Add")
add_btn.config(command=toggle_add)
add_btn.grid(row=0, column=5, padx=2, pady=2, ipadx=8, rowspan=2, sticky=tk.N+tk.S)

canv.pack(side = tk.BOTTOM, fill=tk.BOTH, expand=1)
canv.config(bg='white')
canv.bind("<Configure>", redraw_bg)
canv.bind("<Motion>", mouse_moved)
canv.bind("<Button-1>", mouse_click)
#canv.create_line(0, 0, c_wd, c_ht)
#canv.create_line(0, c_ht, c_wd, 0)
#canv.create_rectangle(3,3,c_wd-4,c_ht-4)
#canv.create_rectangle(6,6,c_wd-7,c_ht-7)
#canv.create_line(0, 100, 200, 0, fill="red", dash=(4, 4))
#canv.create_rectangle(200, 25, c_wd, 75, fill="blue")

bgbox = canv.create_rectangle(0,0,c_wd-1,c_ht-1)
bggrid = []
#Make 1.0m^2 grid
for linex in range(px_per_m, c_wd, px_per_m):
	bggrid.append(canv.create_line(linex, 0, linex, c_ht-1, fill="gray", dash=(dash_len,)))
for liney in range(px_per_m, c_ht, px_per_m):
	bggrid.append(canv.create_line(0, c_ht-1-liney, c_wd-1, c_ht-1-liney, fill="gray", dash=(dash_len,)))


tk.mainloop()
