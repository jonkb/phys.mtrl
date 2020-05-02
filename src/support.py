

class Support:
	sup_types = {
		0: "Fixed",
		1: "Pin",
		2: "Roller"
	}
	def __init__(self, type, tag):
		self.type = type
		self.tag = tag #tag on the canvas
	
	def draw(self, canv, x, y, dir):
		if self.type == 0:
			self.draw_fixed(canv, x, y, dir)
	#draws the lines making up a fixed support
	def draw_fixed(self, canv, x, y, dir):
		if dir==0:#support is above member
			canv.create_line(x-20, y, x+20, y, width=2, tags=self.tag)
			canv.create_line(x-20, y-10, x-10, y, tags=self.tag)
			canv.create_line(x-10, y-10, x, y, tags=self.tag)
			canv.create_line(x, y-10, x+10, y, tags=self.tag)
			canv.create_line(x+10, y-10, x+20, y, tags=self.tag)
		if dir==1:#support is left of member
			canv.create_line(x, y+20, x, y-20, width=2, tags=self.tag)
			canv.create_line(x-10, y+20, x, y+10, tags=self.tag)
			canv.create_line(x-10, y+10, x, y, tags=self.tag)
			canv.create_line(x-10, y, x, y-10, tags=self.tag)
			canv.create_line(x-10, y-10, x, y-20, tags=self.tag)
		if dir==2:#support is below member
			canv.create_line(x-20, y, x+20, y, width=2, tags=self.tag)
			canv.create_line(x-20, y+10, x-10, y, tags=self.tag)
			canv.create_line(x-10, y+10, x, y, tags=self.tag)
			canv.create_line(x, y+10, x+10, y, tags=self.tag)
			canv.create_line(x+10, y+10, x+20, y, tags=self.tag)
		if dir==3:#support is right of member
			canv.create_line(x, y+20, x, y-20, width=2, tags=self.tag)
			canv.create_line(x+10, y+20, x, y+10, tags=self.tag)
			canv.create_line(x+10, y+10, x, y, tags=self.tag)
			canv.create_line(x+10, y, x, y-10, tags=self.tag)
			canv.create_line(x+10, y-10, x, y-20, tags=self.tag)