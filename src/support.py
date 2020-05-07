

class Support:
	sup_types = {
		0: "Fixed",
		1: "Pin",
		#2: "Roller"
		2: "Slot (x)",
		3: "Slot (y)"
	}
	img_w = 20
	def __init__(self, tag):
		self.tag = tag #tag on the canvas
	
	#Returns a binary tuple (x,y,th) for which directions the joint constrains
	def constraints(self):
		return (0,0,0)
	
	#Draws self on the given canvas
	def draw(self, canv, x, y):
		canv.create_oval(x-10,y-10,x+10,y+10)

class Fixed(Support):
	def __init__(self, tag):
		super().__init__(tag)
	
	def constraints(self):
		return (1,1,1)
	
	def draw(self, canv, x, y, dir):
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
			
class Pin(Support):
	def __init__(self, tag):
		super().__init__(tag)
	def constraints(self):
		return (1,1,0)
	def draw(self, canv, x, y, dir):
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
		
class Slot(Support):
	R = Support.img_w/2
	gap = 3
	def __init__(self, tag, isv):
		super().__init__(tag)
		self.isv = isv
	def constraints(self):
		if self.isv:
			return (1,0,0)
		else:
			return (0,1,0)
	def draw(self, canv, x, y):
		if self.isv:
			canv.create_arc(x-self.R, y-3*self.R, x+self.R, y-self.R, extent=180, start=0, 
				width=2, style="arc", tags=self.tag)
			canv.create_arc(x-self.R, y+self.R, x+self.R, y+3*self.R, extent=180, start=180, 
				width=2, style="arc", tags=self.tag)
			canv.create_line(x-self.R, y-2*self.R, x-self.R, y+2*self.R, width=2, tags=self.tag)
			canv.create_line(x+self.R, y-2*self.R, x+self.R, y+2*self.R, width=2, tags=self.tag)
		else:
			#Fill - Maybe add a stipple pattern, though that's not supported for ovals...
			#canv.create_arc(x-3*self.R, y-self.R, x-self.R, y+self.R, extent=180, start=90, 
			#	outline="", fill="white", tags=self.tag)
			#canv.create_arc(x+self.R, y-self.R, x+3*self.R, y+self.R, extent=180, start=-90, 
			#	outline="", fill="white", tags=self.tag)
			#canv.create_rectangle(x-2*self.R, y-self.R, x+2*self.R, y+self.R, 
			#	outline="", fill="white", tags=self.tag)
			#Outline
			canv.create_arc(x-3*self.R, y-self.R, x-self.R, y+self.R, extent=180, start=90, 
				width=2, style="arc", tags=self.tag)
			canv.create_arc(x+self.R, y-self.R, x+3*self.R, y+self.R, extent=180, start=-90, 
				width=2, style="arc", tags=self.tag)
			canv.create_line(x-2*self.R, y-self.R, x+2*self.R, y-self.R, width=2, tags=self.tag)
			canv.create_line(x-2*self.R, y+self.R, x+2*self.R, y+self.R, width=2, tags=self.tag)
		#Pin
		canv.create_oval(x-self.R+self.gap, y-self.R+self.gap, x+self.R-self.gap, 
			y+self.R-self.gap, fill="black", tags=self.tag)