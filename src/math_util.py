import math
import numpy as np
import sympy as sym
#from sympy.sets import Interval

eps = 1e-14
default_sigfig = 6
grav = 9.81 #m/s^2

#round x to n significant figures
def sigfig(x, n=default_sigfig):
	if abs(x) < eps:
		return 0
	return round(x, int(-math.log10(abs(x))+n))
#Takes an iterable and rounds those elements that are very close to int
def eps_round(A):
	for i,e in enumerate(A):
		if (abs(e)+eps/2) % 1 < eps:
			A[i] = round(e)
	return A
def sigfig_iter(A, n=default_sigfig):
	for i,e in enumerate(A):
		A[i] = sigfig(e,n)
	return A
def N_to_kN_str(N, n=default_sigfig):
	k = N/1e3
	return str(sigfig(k, n))+"kN"
def Nm_to_kNm_str(Nm, n=default_sigfig):
	km = Nm/1e3
	return str(sigfig(km, n))+"kN-m"
def Pa_to_MPa_str(P, n=default_sigfig):
	M = P/1e6
	return str(sigfig(M, n))+"MPa"
def coords_str(x,y, n=default_sigfig):
	return "("+str(sigfig(x,n))+","+str(sigfig(y,n))+")"

#Returns the limits of the subdomains of an expression that may contain piecewise
def pw_sdls(f, x):
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
#Returns the discontinuities in the function along the interval
#Supports everything that singularities supports as well as Piecewise
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
	critpts = np.array(sym.solve(sym.Eq(dfdx,0), x)).astype(np.float64)
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

#Symbolically calculate the max and min of the given 2d function on the given domain
#Returns (Max, Max_x, Max_y), (Min, Min_x, Min_y)
def max2d(f, x, y, xdom, ydom):
	x0, x1 = xdom
	y0, y1 = ydom
	#print(61, xdom, ydom)
	f_lam = sym.lambdify((x,y), f, "numpy")
	#print("f: ", f)
	#print("dfdx: ", f.diff(x))
	#print("dfdy: ", f.diff(y))
	grad0 = ( sym.Eq(f.diff(x),0), sym.Eq(f.diff(y),0) )
	#print(63, grad0 )
	#print(64, sym.solve(grad0, (x,y)) )
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
	print(109, critpts)
	f_max = np.max(critvals)
	#print(89, np.where(critvals == f_max))
	xy_max = critpts[np.where(critvals == f_max)[0][0]]
	f_min = np.min(critvals)
	xy_min = critpts[np.where(critvals == f_min)[0][0]]
	return (f_max, xy_max[0], xy_max[1]), (f_min, xy_min[0], xy_min[1])




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
