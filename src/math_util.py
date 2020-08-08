import math
import numpy as np
import sympy as sym
#from sympy.sets import Interval

eps = 1e-14
default_sigfig = 6
grav = 9.81 #m/s^2

#round x to n significant figures
def sigfig(x, n=default_sigfig):
	if isinstance(x, np.ndarray):
		#There seems to be a bug here that causes things like 1.234999999999999
		x_log10 = np.log10(abs(x))	
		x_log10 = np.nan_to_num(x_log10, nan=0, posinf=0, neginf=0)
		round_dec = np.floor(-x_log10+n)
		return np.round(x*10.**round_dec) * 10.**-round_dec
	else:
		if eps_round(x) == 0:
			return 0
		round_dec = int(-math.log10(abs(x))+n)
		return round(x, round_dec)

#Takes an iterable and rounds those elements that are very close to int
#Use for a list, not for an np.array
def eps_round_iter(A):
	for i,e in enumerate(A):
		if (abs(e)+eps/2) % 1 < eps:
			A[i] = round(e)
	return A

def eps_round(n):
	if (abs(n)+eps/2) % 1 < eps:
		n = round(n)
	return n

#Check if two numbers are equal, allowing for small errors
def eps_eq(a, b):
	return abs(a - b) < eps

#Use for a list, not for an np.array
def sigfig_iter(A, n=default_sigfig):
	for i,e in enumerate(A):
		A[i] = sigfig(e,n)
	return A

def N_to_kN_str(N, n=default_sigfig):
	k = float(N/1e3)
	return str(sigfig(k, n)) + "kN"

def Nm_to_kNm_str(Nm, n=default_sigfig):
	km = float(Nm/1e3)
	return str(sigfig(km, n))+"kN-m"

def Pa_to_MPa_str(P, n=default_sigfig):
	M = float(P/1e6)
	return str(sigfig(M, n))+"MPa"

def m_str(m, n=default_sigfig):
	return str(sigfig(float(m), n))+"m"

def coords_str(x,y, n=default_sigfig):
	return "("+str(sigfig(float(x),n))+","+str(sigfig(float(y),n))+")"

#Returns the limits of the subdomains of an expression that may contain piecewise
def pw_sdls(f, x):
	#NOTE: should have used sym.piecewise_fold to simplify this
	lims = []
	#look for all piecewise func in the expr tree 
	#(see https://docs.sympy.org/latest/tutorial/manipulation.html)
	def rec_search(expr):
		nonlocal lims
		if type(expr) == sym.Piecewise:
			sdoms = expr._intervals(x)
			for sdom in sdoms:
				if math.isfinite(sdom[0]):
					lims = np.append(lims, float(sdom[0]))
				if math.isfinite(sdom[1]):
					lims = np.append(lims, float(sdom[1]))
		for arg in expr.args:
			rec_search(arg)
	rec_search(f)
	return lims

#Returns the roots of a piecewise expression
def pw_roots(f, x):
	roots = []
	#This combines any addition or multiplication into the piecewise
	f = sym.piecewise_fold(f)
	def rec_search(expr):
		nonlocal roots
		if type(expr) == sym.Piecewise:
			sdoms = expr.as_expr_set_pairs()
			for sdom in sdoms:
				subf = sdom[0]
				limits = list(sdom[1].boundary)
				#print(65, limits)
				sf_eq0 = sym.Eq(subf, 0)
				if sf_eq0 == True: #function is zero
					numlims = len(limits)
					if numlims == 0: #(-inf, inf)
						#Arbitrarily return zero if the range is infinite
						root = 0
					elif numlims == 1:
						root = limits[0]
					else:
						root = float(limits[1]+limits[0])/2
					roots = np.append(roots, root)
				else:
					roots = np.append(roots, sym.solve(subf, x))
		for arg in expr.args:
			rec_search(arg)
	if type(f) == sym.Piecewise:
		rec_search(f)
	else:
		roots = sym.solve(f, x)
	return roots

#Return a simplified version of the function, restricted to the given domain
#dom = (xmin, xmax)
def f_restrict(f, x, dom, lims=True):
	condition = (x>=dom[0]) & (x<=dom[1]) if lims else (x>dom[0]) & (x<dom[1])
	domf = sym.Piecewise((1, condition))
	return (f*domf).simplify()

#Returns the discontinuities in the function along the interval
#Supports everything that sym.singularities supports as well as Piecewise
def discontinuities(f, x, dom=None):
	discont = np.array([])
	#Try the pre-made singularities function to see if it's enough
	try: sings = sym.singularities(f, x)#, domain=Interval.open(*dom))
	except NotImplementedError: pass #singularities throws exceptions for Piecewise
	else: discont = np.array(list(sings)).astype(np.float64)
	#Now inspect the limits of any piecewise subdomains
	pw_sd_limits = pw_sdls(f,x)
	flam = sym.lambdify(x, f, "numpy")
	for sdlim in pw_sd_limits:
		dx = eps
		f_xl = flam(sdlim - dx)
		f_xll = flam(sdlim - dx*2)
		dfdxl = (f_xl - f_xll)/dx
		#Linear extrapolation 2*dx to the right to skip over the discontinuity
		expected_f_xr = f_xl + dfdxl*(dx*2)
		f_xr = flam(sdlim + dx)
		eps_y = abs(f_xr - expected_f_xr)
		#Anything with 2nd derivitave greater than 1e3 will trigger as discontinuous
		#Anything with a discontinuous derivitave, like sqrt(0), may also fire
		if eps_y > 2e3*eps:
			discont = np.append(discont, sdlim)
	if dom is None:
		return discont
	d_in_dom = []
	for pt in discont:
		if pt > dom[0] and pt < dom[1]:
			d_in_dom = np.append(d_in_dom, pt)
	return d_in_dom
	
#Returns (Max, Max_x), (Min, Min_x)
def max1d(f, x, dom):
	#interval = Interval(*dom)
	#f_max = sym.maximum(f, x, interval)
	#x_max = sym.solve(sym.Eq(f, f_max))[0]
	#f_min = sym.minimum(f, x, interval)
	#x_min = sym.solve(sym.Eq(f, f_min))[0]
	f_lam = sym.lambdify(x, f, "numpy")
	dfdx = f.diff(x)
	try: critpts = sym.solve(dfdx, x)
	except NotImplementedError: critpts = pw_roots(f,x)
	critpts = np.array(critpts).astype(np.float64)
	discont = discontinuities(f, x, dom)
	if len(discont) > 0:
		critpts = np.append(critpts, [discont, discont-eps, discont+eps])
	#critpts = np.array( [pt for pt in critpts if (pt > dom[0] and pt < dom[1])] )
	cpts_in_dom = []
	for pt in critpts:
		if pt > dom[0] and pt < dom[1]:
			cpts_in_dom = np.append(cpts_in_dom, pt)
	critpts = cpts_in_dom
	critpts = np.append(critpts, dom)
	critvals = f_lam(critpts)
	f_max = np.max(critvals)
	x_max = critpts[np.where(critvals == f_max)[0][0]]
	f_min = np.min(critvals)
	x_min = critpts[np.where(critvals == f_min)[0][0]]
	return (f_max, x_max), (f_min, x_min)

#Calculates the max and min of a 2d function on the given domain
#First, it tries to do it analytically, then if that fails, numerically.
#Returns (Max, Max_x, Max_y), (Min, Min_x, Min_y)
def max2d(f, x, y, xdom, ydom, xres=100, yres=100):
	try:
		return max2d_grad(f, x, y, xdom, ydom)
	except Exception as e:
		print(160, e)
		f_lam = sym.lambdify((x,y), f, "numpy")
		xls = np.linspace(*xdom, xres)
		yls = np.linspace(*ydom, yres)
		X, Y = np.meshgrid(xls, yls)
		F = f_lam(X, Y)
		fmax = np.max(F)
		fmin = np.min(F)
		max_coords = np.where(F == fmax)
		min_coords = np.where(F == fmin)
		max_x = max_coords[1][0] * (xdom[1]-xdom[0])/(xres-1) + xdom[0]
		max_y = max_coords[0][0] * (ydom[1]-ydom[0])/(yres-1) + ydom[0]
		min_x = min_coords[1][0] * (xdom[1]-xdom[0])/(xres-1) + xdom[0]
		min_y = min_coords[0][0] * (ydom[1]-ydom[0])/(yres-1) + ydom[0]
		return (fmax, max_x, max_y), (fmin, min_x, min_y)

#Symbolically calculate the max and min of the given 2d function on the given domain
#Returns (Max, Max_x, Max_y), (Min, Min_x, Min_y)
def max2d_grad(f, x, y, xdom, ydom):
	x0, x1 = xdom
	y0, y1 = ydom
	#print(61, xdom, ydom)
	f_lam = sym.lambdify((x,y), f, "numpy")
	#print("f: ", f)
	#print("dfdx: ", f.diff(x))
	#print("dfdy: ", f.diff(y))
	grad0 = ( sym.Eq(f.diff(x),0), sym.Eq(f.diff(y),0) )
	#print(125, grad0 )
	#print(126, sym.solve(grad0, (x,y)) )
	#Likely to fail with Piecewise. :(
	critpts = np.array(sym.solve(grad0, (x,y))).astype(np.float64)
	#print(69, critpts)
	cpts_in_dom = np.empty((0,2), float)
	for pt in critpts:
		if pt[0] > xdom[0] and pt[0] < xdom[1] and pt[1] > ydom[0] and pt[1] < ydom[1]:
			cpts_in_dom = np.append(cpts_in_dom, [pt], axis=0)
	critpts = cpts_in_dom
	#critpts = np.array( [pt for pt in critpts if 
	#	(pt[0] > xdom[0] and pt[0] < xdom[1] and pt[1] > ydom[0] and pt[1] < ydom[1])] )
	#print(77, critpts)
	critvals = []
	if len(critpts) > 0:
		critvals = f_lam(critpts[:,0], critpts[:,1])
	else:
		critvals = []
	N_bound = f.subs(y, y1)
	(N_max, N_max_x), (N_min, N_min_x) = max1d(N_bound, x, (x0, x1))
	critvals = np.append(critvals, [N_max, N_min])
	#print(92, N_bound)
	#print(93, critpts)
	critpts = np.append(critpts, [[N_max_x, y1], [N_min_x, y1]], axis=0)
	#print(94, critpts)
	E_bound = f.subs(x, x1)
	(E_max, E_max_y), (E_min, E_min_y) = max1d(E_bound, y, (y0, y1))
	critvals = np.append(critvals, [E_max, E_min])
	critpts = np.append(critpts, [[x1, E_max_y], [x1, E_min_y]], axis=0)
	#print(99, critpts)
	S_bound = f.subs(y, y0)
	(S_max, S_max_x), (S_min, S_min_x) = max1d(S_bound, x, (x0, x1))
	critvals = np.append(critvals, [S_max, S_min])
	critpts = np.append(critpts, [[S_max_x, y0], [S_min_x, y0]], axis=0)
	#print(104, critpts)
	W_bound = f.subs(x, x0)
	(W_max, W_max_y), (W_min, W_min_y) = max1d(W_bound, y, (y0, y1))
	critvals = np.append(critvals, [W_max, W_min])
	critpts = np.append(critpts, [[x0, W_max_y], [x0, W_min_y]], axis=0)
	#print(226, critpts)
	f_max = np.max(critvals)
	#print(89, np.where(critvals == f_max))
	xy_max = critpts[np.where(critvals == f_max)[0][0]]
	f_min = np.min(critvals)
	xy_min = critpts[np.where(critvals == f_min)[0][0]]
	return (f_max, xy_max[0], xy_max[1]), (f_min, xy_min[0], xy_min[1])

#Accepts any number of points and rotates them the given angle
#th is in degrees
#To be used in graphics, so it's inverted.
def rot_pts(th, *pts, origin=(0,0)):
	#The negative is placed here since down is positive in graphics world.
	thr = - math.radians(th)
	rot = np.array([[math.cos(thr), -math.sin(thr)], [math.sin(thr), math.cos(thr)]])
	rotated = []
	for pt in pts:
		npt = np.vstack(pt) - np.vstack(origin)
		rotated.append(np.dot(rot, npt) + np.vstack(origin))
	return np.array(rotated)



#TESTING:

def test_maxmin():
	x = sym.symbols('x')
	f = sym.sin(x)
	print("max1d('sin(x)', x, [0, pi+.001])")
	print(max1d(f, x, [0, math.pi+.001]))
	
	f1 = 1/x
	print("max1d('1/x', x, [-1, 1])")
	print(max1d(f1, x, [-1, 1]))
	
	f2 = sym.Piecewise((x, x<0), (2-.1*x, True))
	print("max1d('{(x, x<0), (2-.1*x, True)', x, [-1, 1])")
	print(max1d(f2, x, [-1, 1]))
	
	y = sym.symbols('y')
	fxy = sym.sin(x)+1.1*sym.sin(y)
	print("max2d(sin(x)+1.1*sin(y), x, y, [0, math.pi+.001], [0, math.pi+.001])")
	print(max2d(fxy, x, y, [0, math.pi+.001], [0, math.pi+.001]))
#test_maxmin()

def test_solve_pw():
	d, h = sym.symbols('d,h')
	eq0 = sym.Eq(-75.0*sym.Piecewise((-2000.0*d, d < 1.0), (-2000.0, True)) - 150000.0, 0)
	#print(eq0)
	#print(sym.solve(eq0), d)
	
	eq1 = sym.Eq(-75.0*h*sym.Piecewise((-2000.0, d < 1.0), (0, True)), 0)
	#print(eq1)
	#print(sym.solve((eq0, eq1), (d,h)))
	
	x = sym.symbols('x')
	f = sym.Piecewise((x, x<1), (1, x<2), (3-x, True))
	try:
		sym.solve(sym.Eq(sym.diff(f),0), x)
	except NotImplementedError:
		print(211)

#IDEA: first strip Piecewise of anything outside of the domain of interest
#Piecewise((x, x < L and x > 0), (0, True)) --> f = x (on (0,L))
#I think all the boundaries should be considered anyway
#	--> Use f_restrict()

#Another problem is when an interval is constant "solve cannot represent interval solutions"

#test_solve_pw()




