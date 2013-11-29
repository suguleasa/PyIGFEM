import numpy


class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y

#A = Point(0.5, 1.0/3.0)
#B = Point(2.0/3.0, 0.5)
#C = Point(1.0/3.0, 0.5)
#D = Point(0.5, 2.0/3.0)
#rhombus = [ (A.x, A.y), (B.x,B.y), (D.x,D.y), (C.x,C.y) ]

#definition of hexagon corners
#A = Point(1.0/4.0, 1.0/6.0)
#B = Point(3.0/4.0, 1.0/6.0)
#C = Point(5.0/6.0, 1.0/2.0)
#D = Point(3.0/4.0, 5.0/6.0)
#E = Point(1.0/4.0, 5.0/6.0)
#F = Point(1.0/6.0, 1.0/2.0)
#hexagon = [ (A.x,A.y), (B.x,B.y), (C.x,C.y), (D.x,D.y), (E.x,E.y), (F.x,F.y) ]

# rhombus 3
#A = Point(0.5, 0.25)
#B = Point(0.75, 0.5)
#C = Point(0.5, 0.75)
#D = Point(0.25, 0.5)

#rhombus 2
#A = Point(0.5, 0.25)
#B = Point(2.0/3.0, 0.5)
#C = Point(0.5, 0.75)
#D = Point(1.0/3.0, 0.5)

#rhombus 4
#A = Point(0.5, 0.25)
#B = Point(5.0/6.0, 0.5)
#C = Point(0.5, 0.75)
#D = Point(1.0/6.0, 0.5)


#rhombus = [ (A.x, A.y), (B.x,B.y), (C.x,C.y), (D.x,D.y) ]


# SE corner triangle
#domainInclusion = [(0.9,0), (1.0,0.0),(1.0,0.1)]
# NE corner triangle
#domainInclusion = [(0.9,1.0), (1.0,1.0),(1.0,0.9)]
# NW corner triangle
#domainInclusion = [(0.0,0.9),(0.1,1.0),(0.0,1.0)]
# SW corner triangle
#domainInclusion = [(0.0,0.0),(0.1,0.0),(0.0,0.1)]
# SE/NW triangle
#domainInclusion = [(2.0/3.0,0), (1.0,0.0),(1.0,1.0/3.0)]
#domainInclusion = rhombus
#domainInclusion = hexagon

# vertical interface
#domainInclusion = [(0.57,0.0),(0.6,1.0),(1.0,1.0),(1.0,0.0)]
# horizontal interface
#domainInclusion = [(0.0,0.57),(1.0,0.6),(1.0,1.0),(0.0,1.0)]

# semi-circle def:
def f_circle_s(x,y):
  return (x-1.0)**2 + (y-0.0)**2

# circle 1 def:
def f_circle1(x,y):
  return (x-0.25)**2 + (y-0.25)**2

#circle 2 def:
def f_circle2(x,y):
  return (x-0.5)**2 + (y-0.75)**2



# definition of circle
def f_circle(x,y):
  return (x-0.5)**2 + (y-0.5)**2
	#return (x-1.0)**2 + (y-0.0)**2


def uexact(x,loc,k1,k2):
	Amat = numpy.array([[ 2*k1, 0,0,0,0,0],[0,0,0,2*k2,0,0],[0,0,1,0,0,0],[0,0,0,1,1,1],[loc*loc, loc, 1, -loc*loc, -loc, -1],[2*k1*loc, k1, 0, -2*k2*loc, -k2, 0]])
	bvec = numpy.array([1,1,0,0,0,0])
	usol = numpy.linalg.solve(Amat,bvec)

	A = usol[0]
	B = usol[1]
	C = usol[2]
	a = usol[3]
	b = usol[4]
	c = usol[5]

	if x <= loc:
		return A*x*x + B*x + C
	else:
		return a*x*x + b*x + c

#definition of the interface inside the material
def x_fcn(y):
	return 0.5#0.03 * y + 0.57
#	return 0.03 * y + 0.57
#	return 0.65

def interface_fcn(y):
	return 0.03* y+ 0.57
	#return y + 0.85
  #interface is straigh up vertically x = f(y), at location 0.65
#	return 0.65

#definition of the interface inside the material
def y_fcn(x):
	return x

def ufct(x,loc,k1,k2):
	Amat=[[2*k1,0,0,0,0,0],
			[ 0,0, 0, 2*k2, 0, 0],
			[ 0, 0, 1, 0, 0, 0],
			[0, 0, 0, 1, 1, 1],
			[loc**2, loc, 1, -loc**2, -loc, -1],
			[k1*2*loc, k1, 0, -k2*2*loc, -k2, 0]]
	bvec = [-1,-1,0,0,0,0]
	uvec = numpy.linalg.solve(Amat,bvec)
	c1 = uvec[0]
	c2 = uvec[1]
	c3 = uvec[2]
	c4 = uvec[3]
	c5 = uvec[4]
	c6 = uvec[5]
	if x <= loc:
		return c1*x**2 + c2 * x + c3
	else:
		return c4*x**2 + c5 * x + c6
    
  

# Right-hand-side of the PDE:
def rhs(x,y):
	#return 0
	return -1


