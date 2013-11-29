import numpy
import scipy

## This function computes the double integral of f(x,y) dx dy, where 
## x is in (a,b) and y is in (c(x),d(x)), using Legendre-Gaussian nodes ui and weights wi
## 
## Implemented by Elena Caraba on Jan 5, 2013

def quad2d(f,a,b,c,d,ui,wi):
	I = 0
	alpha1 = float(b-a)/2
	beta1 = float(a+b)/2

	nx = len(ui)
	ny = nx
	vj = ui[:]
	wj = wi[:]

	x = []
	y = []

	alpha2 = float( d - c ) /2
	beta2 = float( d + c ) /2
	for i in range(0,nx):

		x = x + [ alpha1 * ui[i] + beta1 ]

		J = 0

		for j in range(0,ny):
			y = y + [ alpha2 * vj[j] + beta2 ] 
			yf = f(x[i],y[j])
			J = J + wj[j] * yf

		if c <= d:
			I = I + wi[i] * alpha2 * J
		else:
			I = 0


	I = float(b-a)/2 * I
	return I



