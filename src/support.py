

class Support:
	sup_types = {
		0: "Fixed",
		1: "Pin",
		#2: "Roller"
		2: "Slot (x)",
		3: "Slot (y)"
	}
	img_w = 20
	def __init__(self, type, tag):
		self.type = type
		self.tag = tag #tag on the canvas
	
	#Returns a binary tuple (x,y,th) for which directions the joint constrains
	def constraints(self):
		#Needs dir for roller...
		#Note, something like a thrust joint is actually (0,2,0), not (0,1,1)
		if self.type == 0:
			return (1,1,1)
		if self.type == 1:
			return (1,1,0)
		if self.type == 2:
			return (0,1,0)
		if self.type == 3:
			return (1,0,0)
	#Returns the degrees of freedom allowed by the support
	#def d_of_f(self):
	#	return 3 - sum(self.constraints())
	def draw(self, canv, x, y, dir):
		self.dir = dir#Move this somewhere else
		if self.type == 0:
			self.draw_fixed(canv, x, y, dir)
		if self.type == 1:
			self.draw_pin(canv, x, y, dir)
		if self.type == 2:
			self.draw_slot(canv, x, y, False)
		if self.type == 3:
			self.draw_slot(canv, x, y, True)
	#draws the lines making up a fixed support
	def draw_fixed(self, canv, x, y, dir):
		self.img_w = 20
		if dir==0:#support is above member
			canv.create_line(x-2*self.img_w, y, x+2*self.img_w, y, width=2, tags=self.tag)
			canv.create_line(x-2*self.img_w, y-self.img_w, x-self.img_w, y, tags=self.tag)
			canv.create_line(x-self.img_w, y-self.img_w, x, y, tags=self.tag)
			canv.create_line(x, y-self.img_w, x+self.img_w, y, tags=self.tag)
			canv.create_line(x+self.img_w, y-self.img_w, x+2*self.img_w, y, tags=self.tag)
		elif dir==1:#support is left of member
			canv.create_line(x, y+2*self.img_w, x, y-2*self.img_w, width=2, tags=self.tag)
			canv.create_line(x-self.img_w, y+2*self.img_w, x, y+self.img_w, tags=self.tag)
			canv.create_line(x-self.img_w, y+self.img_w, x, y, tags=self.tag)
			canv.create_line(x-self.img_w, y, x, y-self.img_w, tags=self.tag)
			canv.create_line(x-self.img_w, y-self.img_w, x, y-2*self.img_w, tags=self.tag)
		elif dir==2:#support is below member
			canv.create_line(x-2*self.img_w, y, x+2*self.img_w, y, width=2, tags=self.tag)
			canv.create_line(x-2*self.img_w, y+self.img_w, x-self.img_w, y, tags=self.tag)
			canv.create_line(x-self.img_w, y+self.img_w, x, y, tags=self.tag)
			canv.create_line(x, y+self.img_w, x+self.img_w, y, tags=self.tag)
			canv.create_line(x+self.img_w, y+self.img_w, x+2*self.img_w, y, tags=self.tag)
		elif dir==3:#support is right of member
			canv.create_line(x, y+2*self.img_w, x, y-2*self.img_w, width=2, tags=self.tag)
			canv.create_line(x+self.img_w, y+2*self.img_w, x, y+self.img_w, tags=self.tag)
			canv.create_line(x+self.img_w, y+self.img_w, x, y, tags=self.tag)
			canv.create_line(x+self.img_w, y, x, y-self.img_w, tags=self.tag)
			canv.create_line(x+self.img_w, y-self.img_w, x, y-2*self.img_w, tags=self.tag)
	def draw_pin(self, canv, x, y, dir):
		pts = [x,y, x-self.img_w/2,y-self.img_w, x+self.img_w/2,y-self.img_w]
		if dir == 1:
			pts[2] = x-self.img_w
			pts[3] = y+self.img_w/2
			pts[4] = x-self.img_w
			pts[5] = y-self.img_w/2
		elif dir == 2:
			pts[3] = y+self.img_w
			pts[5] = y+self.img_w
		elif dir == 3:
			pts[2] = x+self.img_w
			pts[3] = y+self.img_w/2
			pts[4] = x+self.img_w
			pts[5] = y-self.img_w/2
		canv.create_polygon(pts, width=2, fill="white", outline="black", tags=self.tag)
	def draw_slot(self, canv, x, y, isvert):
		R = self.img_w/2
		gap=3
		if isvert:
			canv.create_arc(x-R, y-3*R, x+R, y-R, extent=180, start=0, 
				width=2, style="arc", tags=self.tag)
			canv.create_arc(x-R, y+R, x+R, y+3*R, extent=180, start=180, 
				width=2, style="arc", tags=self.tag)
			canv.create_line(x-R, y-2*R, x-R, y+2*R, width=2, tags=self.tag)
			canv.create_line(x+R, y-2*R, x+R, y+2*R, width=2, tags=self.tag)
		else:
			#Fill - Maybe add a stipple pattern, though that's not supported for ovals...
			#canv.create_arc(x-3*R, y-R, x-R, y+R, extent=180, start=90, 
			#	outline="", fill="white", tags=self.tag)
			#canv.create_arc(x+R, y-R, x+3*R, y+R, extent=180, start=-90, 
			#	outline="", fill="white", tags=self.tag)
			#canv.create_rectangle(x-2*R, y-R, x+2*R, y+R, 
			#	outline="", fill="white", tags=self.tag)
			#Outline
			canv.create_arc(x-3*R, y-R, x-R, y+R, extent=180, start=90, 
				width=2, style="arc", tags=self.tag)
			canv.create_arc(x+R, y-R, x+3*R, y+R, extent=180, start=-90, 
				width=2, style="arc", tags=self.tag)
			canv.create_line(x-2*R, y-R, x+2*R, y-R, width=2, tags=self.tag)
			canv.create_line(x-2*R, y+R, x+2*R, y+R, width=2, tags=self.tag)
		#Pin
		canv.create_oval(x-R+gap, y-R+gap, x+R-gap, y+R-gap, fill="black", tags=self.tag)
