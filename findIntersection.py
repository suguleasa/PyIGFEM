import sys
from meshgeneration import *
import collections
import math
# import auxfcts
from itertools import groupby

#define a class to represent the x and y coordinates of points
class Point:
	def __init__(self,x,y):
		self.x = x
		self.y = y

# intersection between diamond/rhombus inside the domain and the mesh lines
def hexagon_intersection(m,A,B,C,D,E,F,pp):
	
	hx = 1.0 / m
	hy = 1.0 / m 
	ccoords = []
	for i in range(1,m):
		x = i*hx
		y = i*hx

		[x_AB,y_AB] = line_seg(A,B,Point(x,0),Point(x,1))
		if x_AB != 0 and y_AB != 0:
			ccoords = ccoords + [[x_AB,y_AB]]
		[x_AB,y_AB] = line_seg(A,B,Point(0,y),Point(1,y))
		if x_AB != 0 and y_AB != 0:
			ccoords = ccoords + [[x_AB,y_AB]]


		[x_BC,y_BC] = line_seg(B,C,Point(x,0),Point(x,1))
		if x_BC != 0 and y_BC != 0:
			ccoords = ccoords + [[x_BC,y_BC]]
		[x_BC,y_BC] = line_seg(B,C,Point(0,y),Point(1,y))
		if x_BC != 0 and y_BC != 0:
			ccoords = ccoords + [[x_BC,y_BC]]


		[x_DC,y_DC] = line_seg(D,C,Point(x,0),Point(x,1))
		if x_DC != 0 and y_DC != 0:
			ccoords = ccoords + [[x_DC,y_DC]]
		[x_DC,y_DC] = line_seg(D,C,Point(0,y),Point(1,y))
		if x_DC != 0 and y_DC != 0:
			ccoords = ccoords + [[x_DC,y_DC]]


		[x_DE,y_DE] = line_seg(E,D,Point(x,0),Point(x,1))
		if x_DE != 0 and y_DE != 0:
			ccoords = ccoords + [[x_DE,y_DE]]
		[x_DE,y_DE] = line_seg(E,D,Point(0,y),Point(1,y))
		if x_DE != 0 and y_DE != 0:
			ccoords = ccoords + [[x_DE,y_DE]]


		[x_FE,y_FE] = line_seg(F,E,Point(x,0),Point(x,1))
		if x_FE != 0 and y_FE != 0:
			ccoords = ccoords + [[x_FE,y_FE]]
		[x_FE,y_FE] = line_seg(F,E,Point(0,y),Point(1,y))
		if x_FE != 0 and y_FE != 0:
			ccoords = ccoords + [[x_FE,y_FE]]



		[x_AF,y_AF] = line_seg(A,F,Point(x,0),Point(x,1))
		if x_AF != 0 and y_AF != 0:
			ccoords = ccoords + [[x_AF,y_AF]]
		[x_AF,y_AF] = line_seg(A,F,Point(0,y),Point(1,y))
		if x_AF != 0 and y_AF != 0:
			ccoords = ccoords + [[x_AF,y_AF]]


	# sort the coordinates first by y, and then by x
	ccoords = sorted(ccoords, key=lambda x: (x[1],x[0]))
	# remove duplicates from list
	ccoords = [ key for key,_ in groupby(ccoords)]
	pp = numpy.vstack([pp,ccoords])

	return pp 
		
# intersection between diamond/rhombus inside the domain and the mesh lines
def diamond_intersection(m,A,B,C,D,pp):
	
	b = numpy.copy(pp)

	hx = 1.0 / m
	hy = 1.0 / m 
	ccoords = []
	for i in range(0,m):
		x = i*hx
		y = i*hx

		[x_AB,y_AB] = line_seg(A,B,Point(x,0),Point(x,1))
		if x_AB != 0 and y_AB != 0:
			ccoords = ccoords + [[x_AB,y_AB]]
		[x_AB,y_AB] = line_seg(A,B,Point(0,y),Point(1,y))
		if x_AB != 0 and y_AB != 0:
			ccoords = ccoords + [[x_AB,y_AB]]

		[x_BC,y_BC] = line_seg(B,C,Point(x,0),Point(x,1))
		if x_BC != 0 and y_BC != 0:
			ccoords = ccoords + [[x_BC,y_BC]]
		[x_BC,y_BC] = line_seg(B,C,Point(0,y),Point(1,y))
		if x_BC != 0 and y_BC != 0:
			ccoords = ccoords + [[x_BC,y_BC]]


		[x_DC,y_DC] = line_seg(D,C,Point(x,0),Point(x,1))
		if x_DC != 0 and y_DC != 0:
			ccoords = ccoords + [[x_DC,y_DC]]
		[x_DC,y_DC] = line_seg(D,C,Point(0,y),Point(1,y))
		if x_DC != 0 and y_DC != 0:
			ccoords = ccoords + [[x_DC,y_DC]]


		[x_AD,y_AD] = line_seg(A,D,Point(x,0),Point(x,1))
		if x_AD != 0 and y_AD != 0:
			ccoords = ccoords + [[x_AD,y_AD]]
		[x_AD,y_AD] = line_seg(A,D,Point(0,y),Point(1,y))
		if x_AD != 0 and y_AD != 0:
			ccoords = ccoords + [[x_AD,y_AD]]

		if A.y == 0:
			ccoords = ccoords + [[A.x,A.y]]
		if B.y == 0:
			ccoords = ccoords + [[B.x,B.y]]


	# sort the coordinates first by y, and then by x
	ccoords = sorted(ccoords, key=lambda x: (x[1],x[0]))
	# remove duplicates from list
	ccoords = [ key for key,_ in groupby(ccoords)]
	pp = numpy.vstack([pp,ccoords])


	
	#excludedIndex = []
	#nrOrig = (m+1)*(m+1)
	#nrEnr = len(pp) - nrOrig
	#for i in range(0,nrEnr):
	#	for j in range(0,nrOrig):
	#		if pp[j,0] == pp[nrOrig+i,0] and pp[j,1] == pp[nrOrig+i,1]:
	#			excludedIndex = excludedIndex + [nrOrig + i]
#
#	pp = numpy.delete(pp,excludedIndex,0)
	return pp 
		
def find_intersection_semicircle(m,a,b,r,pp):
	hx = 1.0 / m
	hy = 1.0 / m 
	ccoords = []
	for i in range(0,m+1):
		x = i*hx
		temp =  (r**2 - (x-a)**2) 
		if temp >= 0:
			y1 = b - temp ** 0.5
			y2 = b + temp ** 0.5
			if y1>= 0 and y1 <= 1:
				ccoords = ccoords + [[x,y1]]
			if y2>= 0 and y2 <= 1:
				ccoords = ccoords + [[x,y2]]

		y = i * hy
		temp1 = (r**2 - (y-b)**2) 
		if temp1 >= 0:
			x1 = a - temp1 ** 0.5
			x2 = a + temp1 ** 0.5
			if x1 >= 0 and x1<=1 :
				ccoords = ccoords + [[x1,y]]
			if x2 >= 0 and x2<=1 :
				ccoords = ccoords + [[x2,y]]

	# sort the coordinates first by y, and then by x
	ccoords = sorted(ccoords, key=lambda x: (x[1],x[0]))
	#pp = numpy.vstack([pp,ccoords])

	return ccoords#pp 
		
def find_intersection(m,a,b,r,pp):
	hx = 1.0 / m
	hy = 1.0 / m 
	ccoords = []
	for i in range(1,m):
		x = i*hx
		temp =  (r**2 - (x-a)**2) 
		if temp >= 0:
			y1 = b - temp ** 0.5
			y2 = b + temp ** 0.5
			ccoords = ccoords + [[x,y1],[x,y2]]


		y = i * hy
		temp1 = (r**2 - (y-b)**2) 
		if temp1 >= 0:
			x1 = a - temp1 ** 0.5
			x2 = a + temp1 ** 0.5
			ccoords = ccoords + [[x1,y],[x2,y]]

	# sort the coordinates first by y, and then by x
	ccoords = sorted(ccoords, key=lambda x: (x[1],x[0]))
	#pp = numpy.vstack([pp,ccoords])

	return ccoords#pp 

def line_seg(A,B,C,D):

	X4_X3 = (D.x-C.x)
	Y1_Y3 = (A.y-C.y)
	Y4_Y3 = (D.y-C.y)
	X1_X3 = (A.x-C.x)
	X2_X1 = (B.x-A.x)
	Y2_Y1 = (B.y-A.y)

	numerator_a = X4_X3 * Y1_Y3 - Y4_Y3 * X1_X3
	numerator_b = X2_X1 * Y1_Y3 - Y2_Y1 * X1_X3
	denominator = Y4_Y3 * X2_X1 - X4_X3 * Y2_Y1

	if denominator != 0.0:
		u_a = numerator_a / denominator
		u_b = numerator_b / denominator

		# Find the adjacency matrix A of intersecting lines.
		INT_X = A.x + X2_X1 * u_a
		INT_Y = A.y + Y2_Y1 * u_a
		INT_B = (u_a >= 0) & (u_a <= 1) & (u_b >= 0) & (u_b <= 1)
	
		return [round(INT_X * INT_B,15), round(INT_Y * INT_B,15)]
	else:
		return [round(A.x,15),round(A.y,15)]

def main():

	k1 = 1
	k2 = 10

	m = sys.argv[1]	
	m = int(m)
	pp = fake_reg_grid(m,m)
	tt = create_corners_list(m+1,m+1,pp,k1,k2,lambda x: 0.5)

	# Circle-meshgrid intersection
	# Circle defined at center (a,b) and radius r
	a = 0.5
	b = 0.5
	r = 1.0/3.0#.25
	# pp = find_intersection(m,a,b,r,pp)				
	
	# Rhombus - meshgrid intersection
	A = Point(0.5, 1.0/3.0)
	B = Point(2.0/3.0, 0.5)
	C = Point(1.0/3.0, 0.5)
	D = Point(0.5, 2.0/3.0) 
	pp = diamond_intersection(m,A,B,C,D,pp)
			
	tt = number_nodes(m,pp,tt)

def number_nodes(m,pp,tt):
 	nr_of_nodes = (m+1)**2 - 1
	T = m*m

	for e in range(0,T):
		nodes = tt[e]
		coords = pp[nodes,:]

		x0 = coords[0,0]
		x1 = coords[1,0]
		y0 = coords[0,1]
		y1 = coords[2,1]	

		for k in range(nr_of_nodes+1,len(pp)):
			in_poly = point_in_square(pp[k][0],pp[k][1],x0,y0,x1,y1)
			if in_poly == True:
				tt[e] =  tt[e] + [k]#numpy.hstack((tt[e],k))

#	for i in range(0,len(tt) ):
#		if len(tt[i]) == 5 :
#			tt[i] = tt[i][0:4]

	return tt

def point_in_square(x,y,x0,y0,x1,y1):

	if (x0 <= x and x <= x1) and (y0<=y and y <= y1):# and 
		#	(x != x0 and y != y0) and :
		return True
	else:
		return False

if __name__ == "__main__":
	main()
