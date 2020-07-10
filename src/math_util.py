import math
import numpy as np
import sympy as sym
#from sympy.sets import Interval

eps = 1e-14
default_sigfig = 4
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
	critpts = np.array( [pt for pt in critpts if (pt > dom[0] and pt < dom[1])] )
	critpts = np.append(np.append(critpts, dom[0]), dom[1])
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
	f_lam = sym.lambdify((x,y), f, "numpy")
	#print("f: ", f)
	#print("dfdx: ", f.diff(x))
	#print("dfdy: ", f.diff(y))
	grad0 = ( sym.Eq(f.diff(x),0), sym.Eq(f.diff(y),0) )
	#print(63, grad0 )
	#print(64, sym.solve(grad0, (x,y)) )
	critpts = np.array(sym.solve(grad0, (x,y))).astype(np.float64)
	#print(69, critpts)
	critpts = np.array( [pt for pt in critpts if 
		(pt[0] > xdom[0] and pt[0] < xdom[1] and pt[1] > ydom[0] and pt[1] < ydom[1])] )
	#print(72, critpts)
	critvals = f_lam(critpts[:,0], critpts[:,1])
	N_bound = f.subs(y, y1)
	(N_max, N_max_x), (N_min, N_min_x) = max1d(N_bound, x, (x0, x1))
	critvals = np.append(np.append(critvals, N_max), N_min)
	critpts = np.append(np.append(critpts, [[N_max_x, y1]], axis=0), [[N_min_x, y1]], axis=0)
	E_bound = f.subs(x, x1)
	(E_max, E_max_y), (E_min, E_min_y) = max1d(E_bound, y, (y0, y1))
	critvals = np.append(np.append(critvals, E_max), E_min)
	critpts = np.append(np.append(critpts, [[x1, E_max_y]], axis=0), [[x1, E_min_y]], axis=0)
	S_bound = f.subs(y, y0)
	(S_max, S_max_x), (S_min, S_min_x) = max1d(S_bound, x, (x0, x1))
	critvals = np.append(np.append(critvals, S_max), S_min)
	critpts = np.append(np.append(critpts, [[S_max_x, y0]], axis=0), [[S_min_x, y0]], axis=0)
	W_bound = f.subs(x, x0)
	(W_max, W_max_y), (W_min, W_min_y) = max1d(W_bound, y, (y0, y1))
	critvals = np.append(np.append(critvals, W_max), W_min)
	critpts = np.append(np.append(critpts, [[x0, W_max_y]], axis=0), [[x0, W_min_y]], axis=0)
	#print(56, critpts)
	f_max = np.max(critvals)
	#print(89, np.where(critvals == f_max))
	xy_max = critpts[np.where(critvals == f_max)[0][0]]
	f_min = np.min(critvals)
	xy_min = critpts[np.where(critvals == f_min)[0][0]]
	return (f_max, xy_max[0], xy_max[1]), (f_min, xy_min[0], xy_min[1])










def test():
	x = sym.symbols('x')
	f = sym.sin(x)
	print("max1d('sin(x)', x, [0, pi+.001])")
	print(max1d(f, x, [0, math.pi+.001]))
	y = sym.symbols('y')
	fxy = sym.sin(x)+1.1*sym.sin(y)
	print(max2d(fxy, x, y, [0, math.pi+.001], [0, math.pi+.001]))
#test()
