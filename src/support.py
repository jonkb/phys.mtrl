

class Support:
	sup_types = {
		0: "Fixed",
		1: "Pin",
		2: "Roller"
	}
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
	#Returns the degrees of freedom allowed by the support
	#def d_of_f(self):
	#	return 3 - sum(self.constraints())
	def draw(self, canv, x, y, dir):
		self.dir = dir#Move this somewhere else
		if self.type == 0:
			self.draw_fixed(canv, x, y, dir)
	#draws the lines making up a fixed support
	def draw_fixed(self, canv, x, y, dir):
		img_w = 20
		if dir==0:#support is above member
			canv.create_line(x-2*img_w, y, x+2*img_w, y, width=2, tags=self.tag)
			canv.create_line(x-2*img_w, y-img_w, x-img_w, y, tags=self.tag)
			canv.create_line(x-img_w, y-img_w, x, y, tags=self.tag)
			canv.create_line(x, y-img_w, x+img_w, y, tags=self.tag)
			canv.create_line(x+img_w, y-img_w, x+2*img_w, y, tags=self.tag)
		if dir==1:#support is left of member
			canv.create_line(x, y+2*img_w, x, y-2*img_w, width=2, tags=self.tag)
			canv.create_line(x-img_w, y+2*img_w, x, y+img_w, tags=self.tag)
			canv.create_line(x-img_w, y+img_w, x, y, tags=self.tag)
			canv.create_line(x-img_w, y, x, y-img_w, tags=self.tag)
			canv.create_line(x-img_w, y-img_w, x, y-2*img_w, tags=self.tag)
		if dir==2:#support is below member
			canv.create_line(x-2*img_w, y, x+2*img_w, y, width=2, tags=self.tag)
			canv.create_line(x-2*img_w, y+img_w, x-img_w, y, tags=self.tag)
			canv.create_line(x-img_w, y+img_w, x, y, tags=self.tag)
			canv.create_line(x, y+img_w, x+img_w, y, tags=self.tag)
			canv.create_line(x+img_w, y+img_w, x+2*img_w, y, tags=self.tag)
		if dir==3:#support is right of member
			canv.create_line(x, y+2*img_w, x, y-2*img_w, width=2, tags=self.tag)
			canv.create_line(x+img_w, y+2*img_w, x, y+img_w, tags=self.tag)
			canv.create_line(x+img_w, y+img_w, x, y, tags=self.tag)
			canv.create_line(x+img_w, y, x, y-img_w, tags=self.tag)
			canv.create_line(x+img_w, y-img_w, x, y-2*img_w, tags=self.tag)
			
