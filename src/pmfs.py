""" pmfs.py
phys.mtrl file system

Define functions for saving and loading the Lab
"""

import tkinter as tk
from datetime import datetime
import xml.etree.ElementTree as ET
import os
#My files
from member import *
from region import *
from support import *
from load import *


# Constants
ftypes = [('All Files', '*.*'), ('XML Files', '*.xml')]
save_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../data")

def save(lab):
	""" Save the Lab to an xml file
	"""
	data = """<?xml version="1.0" encoding="UTF-8"?>
<!-- This is a lab file created with phys.mtrl -->
<!-- https://github.com/jonkb/phys.mtrl -->\n"""
	data += lab.to_xml() + "\n"
	#print(data)
	now = datetime.now()
	ifile = "lab_"+now.strftime("%Y%m%d-%H%M%S")+".xml"
	# Prompt the user for the save location
	file = tk.filedialog.asksaveasfile(mode='w', 
		initialdir=save_dir, initialfile=ifile,
		filetypes=ftypes, defaultextension=".xml")
	if file is None:
		return
	file.write(data)
	file.close()

def open(lab):
	""" Load a Lab from an xml file
	"""
	# Prompt the user for the file to open
	file = tk.filedialog.askopenfile(mode="r", 
		initialdir=save_dir, defaultextension=".xml")
	if file is None:
		return
	lstr = file.read()
	file.close()
	data = ET.fromstring(lstr)
	if data.tag != "lab":
		print("Error: the selected file is not a lab file.")
		return
	if data.attrib["version"] != lab.version:
		print("Warning: the selected file was made with a different version of phys.mtrl.")
	lab.reset_lab()
	options = data.find("options")
	#Weightless not included here. Include or reset.
	lab.px_per_m = int(options.find("px_per_m").text)
	lab.px_per_kN = float(options.find("px_per_kN").text)
	lab.subdivision = int(options.find("subdivision").text)
	#print(44)
	lab.redraw(int(options.find("c_wd").text), int(options.find("c_ht").text), True)
	members = data.find("members")
	jts = []
	# Place each member
	for mem in members.findall("mem"):
		mem_def = mem.find("def")
		mem_place = mem.find("place")
		x0 = float(mem_place.find("x0").text)
		y0 = float(mem_place.find("y0").text)
		th = float(mem_place.find("th").text)
		L = float(mem_def.find("length").text)
		matl = material_dict(mem_def.find("material").text)
		xsec = mem_def.find("xsec")
		xsec_reg = xsec.attrib["region"]
		xs_param = {param.tag: float(param.text) for param in list(xsec)}
		regions = Region.reg_dict()
		xsection = regions[xsec_reg](**xs_param)
		m = Member(matl, xsection, L)
		lab.place_member(m, xc=x0, yc=y0, th=th)
		# Add every support on this member
		for sup in mem.findall("sup"):
			stype = int(sup.attrib["type"])
			axd = float(sup.find("axd").text)
			th = float(sup.find("th").text)
			lab.place_support(m, stype, axd=axd, th=th)
		# Keep track of every joint on this member to add later, after all the
		#	members have been placed
		for jt in mem.findall("jt"):
			jtype = int(jt.attrib["type"])
			axd = float(jt.find("axd").text)
			th = float(jt.find("th").text)
			jts.append((m, jtype, axd, th))
		# Add every load on this member
		for ld in mem.findall("ld"):
			is_distr = int(ld.attrib["type"])
			if is_distr:
				q0x = float(ld.find("xc0").text)
				q0y = float(ld.find("yc0").text)
				axd0 = float(ld.find("axd0").text)
				q1x = float(ld.find("xc1").text)
				q1y = float(ld.find("yc1").text)
				axd1 = float(ld.find("axd1").text)
				lab.place_distr_load(m, q0x, q0y, axd0, q1x, q1y, axd1)
			else:
				Px = float(ld.find("xc").text)
				Py = float(ld.find("yc").text)
				axd = float(ld.find("axd").text)
				lab.place_load(m, Px, Py, axd)
	# Load the joints
	#	They need to be added after all the members because they join members
	for jt in jts:
		lab.load_joint(*jt)
	print("Loaded Successfully")
	
