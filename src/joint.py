import support as sup


class Joint:
	#Each Joint instance should inherit also from a Support.
	#This sounds dumb, but it kind of makes sense to me.
	#Again, maybe not the best way to do this... IDK
	def __init__(self, m0, m1, axd0=None, axd1=None):
		self.m0 = m0
		self.m1 = m1
		if axd0 is None or axd1 is None:
			intsx = self.m0.intersection(m1)
			if intsx is None: print("ERROR: Non-intersecting members.")
			(_, axd0, axd1) = intsx
		self.axd0 = axd0
		self.axd1 = axd1
	
	def __str__(self):
		assert self.tag is not None
		return "{} joint ({}) @axd0={}".format(type(self).__name__, 
			self.tag, sup.m_u.m_str(self.axd0))
	
	def other_mem(self, m):
		return self.m0 if m is self.m1 else self.m1
		
	def axd(self, m):
		if m is self.m0:
			return self.axd0
		if m is self.m1:
			return self.axd1
		return None
		
	def jtype(self):
		return self.stype

class Fixed(Joint, sup.Fixed):
	def __init__(self, tag, m0, m1, axd0=None, axd1=None):
		Joint.__init__(self, m0, m1, axd0, axd1)
		sup.Fixed.__init__(self, tag)
	def draw(self, canv, x, y):
		r = self.img_w * 0.5
		canv.create_rectangle(x-r, y-r, x+r, y+r, tags=self.tag)
		canv.create_line(x-r, y-r, x+r, y+r, tags=self.tag)
		canv.create_line(x-r, y+r, x+r, y-r, tags=self.tag)

class Pin(Joint, sup.Pin):
	def __init__(self, tag, m0, m1, axd0=None, axd1=None):
		Joint.__init__(self, m0, m1, axd0, axd1)
		sup.Pin.__init__(self, tag)
	def draw(self, canv, x, y):
		ro = self.img_w/2
		ri = ro - sup.Slot.gap
		canv.create_oval(x-ri, y-ri, x+ri, y+ri, fill="black", tags=self.tag)
		canv.create_oval(x-ro, y-ro, x+ro, y+ro, width=2, tags=self.tag)

class Slot(Joint, sup.Slot):
	def __init__(self, tag, m0, m1, axd0=None, axd1=None, th=0):
		Joint.__init__(self, m0, m1, axd0, axd1)
		sup.Slot.__init__(self, tag, th=th)

class Thrust(Joint, sup.Thrust):
	def __init__(self, tag, m0, m1, axd0=None, axd1=None, th=0):
		Joint.__init__(self, m0, m1, axd0, axd1)
		sup.Thrust.__init__(self, tag, th=th)
	
