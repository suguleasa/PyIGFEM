import sys
import numpy

def grid_coordinates(m,n,loc_x_fcn):

	p = numpy.zeros(((n+1)*(m+1),2))
	tmp = 0
	# in the y direction
	for i in range(0,m+1):

		# since there is no interface along y, the spacing hy is 1/m
		# where m is the number of elements in the y direction
		y = i * 1.0 / m

		# compute the x coordinate of the interface
		xloc = lambda e: loc_x_fcn(e) 

		# in the x direction
		for j in range(0,n+1):
			# on the left side of the vertically slanted interface
			if j<n/2.0:
				x = 0 + j * ( xloc(y) - 0 ) / (n/2.0)
			# on the right side of the vertically slanted interface
			else:
				x = xloc(y)  + (j-n/2) * ( 1 - xloc(y) ) / (n/2.0)
			p[tmp,:] = [x,y]
			tmp = tmp + 1

	return p

def fake_reg_grid(m,n):
	p = numpy.zeros(((n+1)*(m+1),2))
	tmp = 0
	for i in range(0,n+1):
		y = i * 1.0/n
		for j in range(0,m+1):
			x = j * 1.0/m
			p[tmp,:] = [x,y]
			tmp = tmp + 1
	return p			 

def g_of_y(y):
	return y + 0.9 # SE
	
def f_of_x(x):
	return x - 0.9 # SE
	
def fake_data_triangle_right_corner(pvec,hx,hy):
#lower right corner

	for i in range(0,len(pvec)):
		x_co = pvec[i,0]
		y_co = pvec[i,1]
		if (pvec[i,0] <= g_of_y(y_co) and g_of_y(y_co) <= pvec[i+1,0]):
			x_pi = g_of_y(y_co)#0.85
			y_pi = y_co
			newRow = [x_pi,y_pi]
			pvec = numpy.vstack([pvec,newRow])
		if (pvec[i,1] <= f_of_x(x_co) and f_of_x(x_co) <= pvec[i,1]+hy and f_of_x(x_co)<=1):
			x_pi = x_co
			y_pi = f_of_x(x_co)#0.15
			newRow = [x_pi,y_pi]
			pvec = numpy.vstack([pvec,newRow])
	return pvec

def check_interface_alignment_H(p,loc_y_fcn,m):
	xr = p[:,0]

	y_fcn_array = numpy.array(loc_y_fcn(xr))
	locIndex = numpy.nonzero( numpy.in1d(xr,y_fcn_array))

	# the interface is not along element edges
	if numpy.array(locIndex).size == 0:
		for i in range(0, len(p)):
			if (p[i,1] <= loc_y_fcn(0.0) and loc_y_fcn(1.0) <= p[i%m+m*(i/m+1),1]): #or ( (p[i,1] <= loc_y_fcn(p[i,0]) and loc_y_fcn(p[i,0]) <= p[i%m + m *(i/m+1), 1]) and ( p[i+1,1] <= loc_y_fcn( p[i+1,0]) and loc_y_fcn( p[i+1,0]) <= p[i%m + m *(i/m+1)+1,1] ) and loc_y_fcn(p[i+1,0]<=1.0) ):
				x_pi = p[i,0]
				y_pi = loc_y_fcn(x_pi)
				newRow = [x_pi,y_pi]
				p = numpy.vstack([p,newRow])
	return p
	

def check_interface_alignment(p,k1,k2,loc_x_fcn):
	
	xr = p[:,0]

	# search to see if there are any points in the grid that have the x coordinate
	# the same with the x coordinate of the interface -  assuming interface is vertical
	x_fcn_array = numpy.array(loc_x_fcn(xr))
	locIndex = numpy.nonzero( numpy.in1d(xr,x_fcn_array))

	# the interface is not along element edges
	if numpy.array(locIndex).size == 0 and k1!=k2:
		for i in range(0, len(p)):
			if p[i,0] < loc_x_fcn(0.0) and loc_x_fcn(1.0) < p[i+1,0]:
				y_pi = p[i,1]
				x_pi = loc_x_fcn(y_pi)
				newRow = [x_pi,y_pi]
				p = numpy.vstack([p,newRow])
	return p

def fake_corners_list(m,n,p):
	enr_node_iter = 0
	t = []

	for j in range(1,n):
		for i in range(1,m):
			c1 = i + (j-1) * m
			c2 = c1 + 1
			c3 = i + j * m
			c4 = c3 + 1

			c1 -= 1
			c2 -= 1
			c3 -= 1
			c4 -= 1
	
			part4 = p[c3,0] <= g_of_y(p[c3,1]) and g_of_y(p[c4,1]) <= p[c4,0]
			part1 = p[c1,0] <= g_of_y(p[c1,1]) and g_of_y(p[c2,1]) <= p[c2,0]
			part2 = p[c1,1] <= f_of_x(p[c1,0]) and f_of_x(p[c3,0]) <= p[c3,1]
			part3 = p[c2,1] <= f_of_x(p[c2,0]) and f_of_x(p[c4,0]) <= p[c4,1]
			comp = (part1 or part4) and ( part2 or part3)

			if comp or ( (part1 and part4) or (part2 and part3) ):
				en1 = m*n + enr_node_iter
				en2 = en1 + 1
				t = t + [[c1,c2,c4,c3, en1,en2]]
				enr_node_iter += 1
			else:
				t = t + [[c1,c2,c4,c3]]

	return t				

def pt_inside_element(c1,c2,c3,c4,p,pt):
	inside = 0
	c1x = p[c1,0]
	c2x = p[c2,0]
	c1y = p[c1,1]
	c4y = p[c4,1]
	ptx = pt[0]
	pty = pt[1]

	if (f_of_x(c1x) <= ptx and ptx <= f_of_x(c2x)) and ( g_of_y(c1y) <= pty and pty <= g_of_y(c4y)):
		inside = 1
	return inside
 
def create_corners_list_H(m,n,p,loc_y_fcn):

	enr_node_iter = 0;
	t = []
	for j in range(1,n):
		for i in range(1,m):
			c1 = i + (j-1) * m 
			c2 = c1 + 1  
			c3 = i + j * m 
			c4 = c3 + 1 

			# Start with indexing at 0
			c1 -= 1
			c2 -= 1
			c3 -= 1
			c4 -= 1

			if ( (p[c1,1] <= loc_y_fcn(0.0) and loc_y_fcn(1.0) <= p[c3,1]) or (p[c2,1] <= loc_y_fcn(0.0) and loc_y_fcn(1.0) <= p[c4,1]) ):
				en1 = m*n + enr_node_iter
				en2 = en1 + 1
				t = t + [[c1,c2,c4,c3, en1,en2]]
				enr_node_iter += 1
			else:
				t = t + [[c1,c2,c4,c3]]
	return t

# create a list with all the corners of an element per row
def create_corners_list(m,n,p,k1,k2,loc_x_fcn):

	enr_node_iter = 0;
	t = []
	for j in range(1,n):
		for i in range(1,m):
			c1 = i + (j-1) * m 
			c2 = c1 + 1  
			c3 = i + j * m 
			c4 = c3 + 1 

			# Start with indexing at 0
			c1 -= 1
			c2 -= 1
			c3 -= 1
			c4 -= 1

			# on the elements with enrichment nodes en1, en2, add the nodes to the element's corners list: [c1, c2, c4, c3, en1, en2] 
			if ( (p[c1,0] < loc_x_fcn(0.0) and loc_x_fcn(1.0) < p[c2,0]) or (p[c3,0] < loc_x_fcn(0.0) and loc_x_fcn(1.0) < p[c4,0]) ) and k1!=k2:
				en1 = m*n + enr_node_iter
				en2 = en1 + 1
				t = t + [[c1,c2,c4,c3, en1,en2]]
				enr_node_iter += 1
			else:
				t = t + [[c1,c2,c4,c3]]

	return t

