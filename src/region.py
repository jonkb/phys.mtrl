import math

#A closed 2D Region.
#Used for cross sections of members.
class Region:
	regions = ("circle", "rectangle", "I-beam", "annulus")

	#Ix,Iy - second moment of area about (horizontal,vertical) axes
	def __init__(self, area, Ix, Iy):
		self.area = area
		self.Ix = Ix
		self.Iy = Iy
	def __repr__(self):
		return "Region({},{},{})".format(self.area, self.Ix, self.Iy)
	@property
	def Ip(self):
		return self.Ix+self.Iy
	@property
	def Imin(self):
		#This isn't necessarily correct
		#If x and y aren't principle axes, this needs to be calculated
		#https://calcresource.com/moment-of-inertia-rotation.html
		return min(self.Ix, self.Iy)

class Circle(Region):
	def __init__(self, radius):
		self.radius = radius
		#Do I need the super() thing?
	def __str__(self):
		return "circle with radius=" + str(self.radius)
	@property #getter
	def area(self):
		return math.pi*self.radius**2
	@property
	def Ix(self):
		return math.pi/4*self.radius**4
	@property
	def Iy(self):
		return self.Ix
	#First moment of area above point of interest. y1 is measured from centroid.
	#For a circle, this really only matters at y1=0.
	def Q(self, y1):
		if y1 == 0:
			return 2/3*self.radius**3
		else:
			return "VQ/Ib is invalid anywhere but the centroid"
	#returns Q/Ib for the cross section
	def Q_div_Ib(self, y1):
		if y1 == 0:
			return self.Q(0) / (self.Ix*2*self.radius)
		else:
			return "VQ/Ib is invalid anywhere but the centroid"

class Rectangle(Region):
	def __init__(self, base, height):
		self.base = base
		self.height = height
	def __str__(self):
		return "rectangle with base=" +str(self.base)+ ", and height=" +str(self.height)
	@property
	def area(self):
		return self.base * self.height
	@property
	def Ix(self):
		return self.base * self.height**3 / 12
	@property
	def Iy(self):
		return self.height * self.base**3 / 12
	#First moment of area above point of interest. y1 is measured from centroid.
	def Q(self, y1):
		if abs(y1) > self.height/2:
			return "Out of domain"
		h1 = self.height/2-y1
		return self.base*h1 * (y1+h1/2)
	#returns Q/Ib for the cross section
	def Q_div_Ib(self, y1):
		Q = self.Q(y1)
		if isinstance(Q, str):
			return Q
		return Q / (self.Ix * self.base)

#Wide Flanged I-beam
class W_F_I(Region):
	#Major depth(height), flange width, flange thickness, web thickness
	def __init__(self, depth, width, tflg, tweb):
		self.depth = depth
		self.width = width
		self.tflg = tflg
		self.tweb = tweb
	def __str__(self):
		return "wide flanged I-beam with height=" +str(self.depth)
	@property
	def area(self):
		return 2*self.tflg*self.width + self.tweb*(self.depth-2*self.tflg)
	@property
	def Ix(self):
		Ix1 = self.tweb/3*(self.depth/2-self.tflg)**3
		Ix2 = self.width/3*((self.depth/2)**3 - (self.depth/2-self.tflg)**3)
		return 2*(Ix1+Ix2)
	@property
	def Iy(self):
		return (self.depth-2*self.tflg)*self.tweb**3/12 + self.tflg*self.width**3/6
	#First moment of area above point of interest. y1 is measured from centroid.
	def Q(self, y1):
		if abs(y1) > self.depth/2:
			return "Out of domain"
		if abs(y1) > self.depth/2 - self.tflg:
			return "VQ/Ib is invalid in the flanges"
		dw1 = self.depth/2-self.tflg-y1
		Aflg = self.width*self.tflg
		return Aflg*(self.depth-self.tflg)/2 + dw1*self.tweb * (y1+dw1/2)
	#returns Q/Ib for the cross section
	def Q_div_Ib(self, y1):
		Q = self.Q(y1)
		if isinstance(Q, str):
			return Q
		return Q / (self.Ix * self.tweb)

class Annulus(Region):
	def __init__(self, ro, ri):
		self.ro = ro
		self.ri = ri
	def __str__(self):
		return "annulus with outer radius=" +str(self.ro)+", and inner radius=" +str(self.ri)
	@property
	def area(self):
		return math.pi * (self.ro**2 - self.ri**2)
	@property
	def Ix(self):
		return math.pi/4* (self.ro**4 - self.ri**4)
	@property
	def Iy(self):
		return self.Ix
	#First moment of area above point of interest. y1 is measured from centroid.
	#For a circle, this really only matters at y1=0.
	def Q(self, y1):
		if y1 == 0:
			return 2/3*(self.ro**3 - self.ri**3)
		else:
			return "VQ/Ib is invalid anywhere but the centroid"
	#returns Q/Ib for the cross section
	def Q_div_Ib(self, y1):
		if y1 == 0:
			return self.Q(0) / (self.Ix * 2*(self.ro-self.ri))
