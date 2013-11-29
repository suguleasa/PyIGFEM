import math
import numpy

# Commputing the Legendre-Gauss weights w and nodes x on an interval [a,b]
# with truncation order N
# Code based on the Matlab code written by Greg von Winckel
# written by Elena Caraba - Jan 5, 2013


def lgwt(n,a,b):
	N = n
	N = N - 1
	N1 = N + 1
	N2 = N + 2;
	xu = numpy.linspace(-1,1,N1)	
	xu = xu.transpose()

	t = range(0,N+1)
	cos_arg = numpy.array(t)
	cos_arg = cos_arg.transpose()
	cos_arg = (cos_arg * 2 + 1) * math.pi / float(2 * N + 2)
	sin_arg = math.pi * xu * float(N)/N2

	# initial guess:
	y = numpy.cos(cos_arg) + 0.27 / float(N1) * numpy.sin(sin_arg)

	# Legendre-Gauss Vandermode Matrix
	L = numpy.zeros((N1,N2))
	
	# Dertivative of LGVM
	Lp = numpy.zeros((N1,2))

	# Compute the zeros of the N+1 Legendre Polynomial
	# using the recursion relation and the Newton-Raphson method
	y0=2;

	searching  = (abs(y-y0)).max() > 10**(-16)


	# Iterate until new points are uniformly within epsilon of old points
	while searching:
		L[:,0] = 1;
		Lp[:,0] = 0;

		L[:,1] = y;
		Lp[:,1] = 1;
	
		for k in range(1,N1):
			L[:,k+1] = ( (2*(k+1)-1) * y * L[:,k] -((k+1)-1) * L[:,k-1] )/float(k+1)

		y.astype('float64')
		tLp= (N2)*( L[:,N1-1] - y * L[:,N2-1] ) / (1-y**2)
		y0 = y;
		tLp.astype('float64')
		y = y0 - L[:,N2-1]/(tLp)

		searching  = (abs(y-y0)).max()>2.2204*10**(-16)

	# computing the nodes:
	x = (a * (1-y) + b*(1+y)) / float(2);
	# computing the weights:
	w = (b-a) / ( (1-y**2) * tLp**2) * (N2/float(N1))**2;
	return [x,w]
