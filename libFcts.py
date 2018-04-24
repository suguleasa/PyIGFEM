import SimpleITK as sitk
from globalVars import *
from math import sqrt, floor
from numpy import *
import numpy

class Coordinate(object):
    def __init__(self,x=-1,y=-1):
        self.x = x
        self.y = y
            
def search_in(my_list,pi,pj,inImage):
    Lk_list1 = [[x[0],x[1]] in [[pi,pj],] for x in my_list] #[True, False, False, True, etc]
    Lk_list2 = [[x[0],x[1]] in [[pj,pi],] for x in my_list] #[True, False, False, True, etc]
    if True in Lk_list1: #if we found it in the list:
        Lk_ind = Lk_list1.index(True)
        Lk = my_list[Lk_ind][2]
    else: 
        if True in Lk_list2:
            Lk_ind = Lk_list2.index(True)
            Lk = my_list[Lk_ind][2]
        else:
            Lk = linear_search(inImage,pi,pj)
            my_list.append([pi,pj,Lk])
    return Lk

## This function draws a line between two points in space: pinit and pend
def draw_line_normals(image,pinit,pend):
    pix_col = 2
    return
    if (pinit.x == pend.x):
        # print "vertical line"
        if (pinit.y <= pend.y):
            startLine = pinit.y;
            endLine = pend.y;
        else:
            startLine = pend.y;
            endLine = pinit.y;
        # draw the line
        for i in range(startLine,endLine+1):
            image.SetPixel(pinit.x,i,0,pix_col);

    elif pinit.y == pend.y:
    #    print "horizontal line"
        if (pinit.x <= pend.x):
            startLine = pinit.x;
            endLine = pend.x;
        else:
            startLine = pend.x;
            endLine = pinit.x;
        # draw the line
        for i in range(startLine,endLine+1):
            image.SetPixel(i,pinit.y,0,pix_col);
    else:
        #print "oblique line"
        dx = pend.x - pinit.x;
        dy = pend.y - pinit.y;

        # line expressed as y = f(x)
        if abs(dx) > abs(dy):
        #    print 'y=f(x)'
            slope = float(dy)/dx;
            if (pinit.x < pend.x):
                bbegin = pinit.x;
                eend = pend.x;
            else:
                bbegin = pend.x;
                eend = pinit.x;
            for i in range(bbegin,eend+1):
#                type(i) is int
                yloc = slope * (i - pend.x) + pend.y;
                yloc = int(yloc)
#                type(yloc) is int
                image.SetPixel(i,yloc,0,pix_col);

        # line expressed as x = f(y)
        else:
        #    print 'x = f(y)'
            slope = float(dx)/dy;
            if (pinit.y < pend.y):
                bbegin = pinit.y;
                eend = pend.y;
            else:
                bbegin = pend.y;
                eend = pinit.y;
            for i in range(bbegin,eend+1):
                xloc = slope * ( i - pinit.y) + pinit.x;
                xloc = int(xloc);
                image.SetPixel(xloc,i,0,pix_col);
                
## This function draws a line between two points in space: pinit and pend
def draw_line(image,pinit,pend):
    pix_col = 1
    if (pinit.x == pend.x):
		# print "vertical line"
		if (pinit.y <= pend.y):
			startLine = pinit.y;
			endLine = pend.y;
		else:
			startLine = pend.y;
			endLine = pinit.y;
		# draw the line
		for i in range(startLine,endLine+1):
			image.SetPixel(pinit.x,i,0,pix_col);

    elif pinit.y == pend.y:
	#	print "horizontal line"
		if (pinit.x <= pend.x):
			startLine = pinit.x;
			endLine = pend.x;
		else:
			startLine = pend.x;
			endLine = pinit.x;
		# draw the line
		for i in range(startLine,endLine+1):
			image.SetPixel(i,pinit.y,0,pix_col);
    else:
		#print "oblique line"
		dx = pend.x - pinit.x;
		dy = pend.y - pinit.y;

		# line expressed as y = f(x)
		if abs(dx) > abs(dy):
		#	print 'y=f(x)'
			slope = float(dy)/dx;
			if (pinit.x < pend.x):
				bbegin = pinit.x;
				eend = pend.x;
			else:
				bbegin = pend.x;
				eend = pinit.x;
			for i in range(bbegin,eend+1):
#				type(i) is int
				yloc = slope * (i - pend.x) + pend.y;
				yloc = int(yloc)
#				type(yloc) is int
				image.SetPixel(i,yloc,0,pix_col);

		# line expressed as x = f(y)
		else:
		#	print 'x = f(y)'
			slope = float(dx)/dy;
			if (pinit.y < pend.y):
				bbegin = pinit.y;
				eend = pend.y;
			else:
				bbegin = pend.y;
				eend = pinit.y;
			for i in range(bbegin,eend+1):
				xloc = slope * ( i - pinit.y) + pinit.x;
				xloc = int(xloc);
				image.SetPixel(xloc,i,0,pix_col);


## This function returns 1 if the pixel value is within a given
## bin with boundaries: [lowLim,highLim], and 0 otherwise
def in_bin_i(pxVal, lowLim, highLim):
    return (lowLim <= pxVal and pxVal < highLim)

## This function checks if two given pixels belong in the same bin
def is_in_same_bin(pxVal1, pxVal2):
    isHomogeneous = False

    for i in range(1,len(binBnd)):
        lowLim = binBnd[i-1]
        highLim = binBnd[i]
        isHomogeneous = isHomogeneous or ( in_bin_i(pxVal1, lowLim, highLim) and in_bin_i(pxVal2, lowLim, highLim))
    return isHomogeneous



## Check if the ends of a line are from the same bin
def ends_in_same_bin(image, p1, p2):
	pxVal1 = image.GetPixel(int(p1.x), int(p1.y));
	pxVal2 = image.GetPixel(int(p2.x), int(p2.y));
	val1 = pxVal1 == pxVal2
	val2 = is_in_same_bin(pxVal1, pxVal2)
	return (val1 or val2)


## Evaluate the Newton polynomial p at x
def newt_fct_eval(a,xData,x):
    n = len(xData) - 1  # Degree of polynomial
    p = a[n]
    for k in range(1,n+1):
        p = a[n-k] + (x - xData[n-k]) * p
    return p 
 
# Compute the coefficients in the Newton polynomial in the a vector
def newt_coef(xData,yData):
    m = len(xData)  # Number of data points
    a = yData[:]

    for k in range(1,m):
        for i in range(m-1,k-1,-1):
            a[i] = (a[i] - a[i-1]) / (xData[i] - xData[i-k]);

    return a

## This function finds the midpoint between two points in space: p1 and p2
def find_mid_point(p1, p2):
    xMid = int( (p1.x + p2.x)/2.0 )
    yMid = int( (p1.y + p2.y)/2.0 )
    return Coordinate(xMid,yMid)


## This function checks if the four corners belong to the same
## pixel value (same bin that has a range of pixel values)
## check if the 4 corners are all in the same bin of colors or not
def four_corners_test(pxVal1, pxVal2, pxVal3, pxVal4):

    isHomogeneous = False
    for i in range(len(binBnd)-1):
        lowLim = binBnd[i-1]
        highLim = binBnd[i]
        isHomogeneous = isHomogeneous or (in_bin_i(pxVal1,lowLim,highLim) and in_bin_i(pxVal2,lowLim,highLim) and in_bin_i(pxVal3,lowLim,highLim) and in_bin_i(pxVal4,lowLim,highLim))

    return isHomogeneous

## This function draws a polynomial between points p1 and p2
def draw_curve(image,p1,p2,d,x):
	for i in range(p1.x,p2.x+1):
		yloc = (newt_fct_eval(d,x,i)); 
		yloc = int(yloc)
		image.SetPixel(i,yloc,0,1);

## This function computes the probability that there is an inclusion
## of a different material inside thie element.
## Given that all 4 corners are found to be homogeneous, or to belong in to
## the same bin, then is there any inclusion inside this element that is
## larger than a certain area? If so, return true, otherwise return false
def has_inclusions(image,p1,p2,p3,p4):

    pxVal1 = image.GetPixel(p1.x,p1.y);

    xHigh = p2.x 
    xLow = p1.x	
    yHigh = p4.y
    yLow = p1.y
    areaElem = abs( (p4.y - p1.y) * (p2.x - p1.x))
    nr_samples = int(log(PROB)/log (abs(areaElem - AREA_INCLUSION)/areaElem))
	
    for i in range (1,nr_samples):
        rx = random.randint(xLow,xHigh)
        ry = random.randint(yLow,yHigh)
        samplePixel = image.GetPixel(rx,ry,0)
        if (0 == is_in_same_bin(pxVal1,samplePixel)):
            return True

    return False

## Find distance between two points
def find_distance(p1, p2):
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    return ( sqrt(dx**2 + dy**2) )


## Doing a log search along a line to discover possible intersections
## with interfaces
## Note: it returns the last pixel in material 1, and the next one is 
## in material 2.
def ilog_search(image, begin, end):
    mid = find_mid_point(begin, end)
    dist = find_distance(begin, mid)
    searching = True

    while searching:
        if ends_in_same_bin(image, begin, mid):
            begin.x = mid.x
            begin.y = mid.y
            mid.x = (begin.x + end.x)/2
            mid.y = (begin.y + end.y)/2
        else:
            end.x = mid.x
            end.y = mid.y
            mid = find_mid_point(begin, end)

        dist = find_distance(begin, end)
        searching = (dist > 1.0)

    return mid 

## Finding the coordinates of point at the intersection between two lines
## that are given by end points: L1, L3 and p1, p3 
## x= f(y)
def	line_line_intersection_x(image,p1,p3,L1,L3):
	mL = float( (L3.x - L1.x)) / (L3.y - L1.y) ;
	mP = float( (p3.x - p1.x)) / (p3.y - p1.y);
	
	yC = ( mP * p3.y - mL * L3.y - p3.x + L3.x ) / (mP - mL);
	xC = mL * (yC - L3.y) +  L3.x;
	ptIntersection = Coordinate( floor(xC), floor(yC) );

	return ptIntersection

## Finding the coordinates of point at the intersection between two lines
## that are given by the end points: L2, L4 and p2, p4
## y = f(x)
def line_line_intersection_y(image,p2,p4,L2,L4):
	mL = float( (L4.y - L2.y)) / (L4.x - L2.x);
	mP = float( (p4.y - p2.y)) / (p4.x - p2.x);

	xC = ( mP * p2.x - mL * L2.x - p2.y + L2.y ) / (mP - mL);
	yC = mP * (xC - p2.x) +  p2.y;

	ptIntersection = Coordinate( floor(xC), floor(yC) );

	return ptIntersection

def linear_search(image,bbegin,eend):
		
		list_nodes = []
		begin = Coordinate(bbegin.x,bbegin.y)
		end = Coordinate(eend.x,eend.y)
		
		old = Coordinate(bbegin.x,bbegin.y)
		dist = find_distance(begin,end)
		if bbegin.x == eend.x and dist>2: # vertical line
			
			next = Coordinate(begin.x,begin.y+1)
			while next.y <= end.y:
				if not(ends_in_same_bin(image,next,old)):
					
					list_nodes = list_nodes + [Coordinate(next.x,next.y)]
				old = Coordinate(next.x,next.y)
				next = Coordinate(next.x,next.y+1)
		
		elif bbegin.y == eend.y and dist>2: # horizontal line
		
			next = Coordinate(begin.x+1,begin.y)			
			
			while next.x <= end.x :

				if not(ends_in_same_bin(image, next,old)):
					list_nodes = list_nodes + [Coordinate(next.x,next.y)]
                    
				old = Coordinate(next.x,next.y)
				next = Coordinate(next.x+1,next.y)		
	
		return list_nodes

def log_search(image,bbegin,eend):
		
		begin = Coordinate(bbegin.x,bbegin.y);
		end = Coordinate(eend.x,eend.y);
		mid = Coordinate();

		mid.x = (begin.x + end.x)/2.0;
		mid.y = (begin.y + end.y)/2.0;
		dist = find_distance(begin,end)
		while dist>2 and not(ends_in_same_bin(image,begin,end)):
			if ends_in_same_bin(image,begin, mid):
				begin.x = mid.x
				begin.y = mid.y;
				mid.x = (begin.x + end.x)/2.0;
				mid.y = (begin.y + end.y)/2.0;
			else:
				end.x = mid.x;
				end.y = mid.y;
				mid.x = (begin.x + end.x)/2.0;
				mid.y = (begin.y + end.y)/2.0;
		
			dist = find_distance(begin,end)

		mid.x = int(mid.x)
		mid.y = int(mid.y)
		return mid

## case 1: 3:1 -- P1 the outsider
def case_NW_polynomial_test(image,p1,p2,p3,p4,L1,L4, poly_opt = 0):
    x_is_F_of_y = False;
    		
    if POL_APPROX != 3:    
        
        if ( find_distance(L1,L4) <= 2.0):
    		pt = Coordinate(-1,-1);
    		vecCoord = [];
    		vecCoord.append(pt);
    		return vecCoord;
    	    	
    	# LINEAR POLYNOMIAL
    	# Step 1. Search for the interface along the 45 degree line
        A = log_search(image,p1,p3);
    
    	# Step 2. Find intersection of line between L1,L4 and the 45 degree line
        if (abs(L1.x-p1.x) < abs(L4.y - p1.y)):
    		B = line_line_intersection_x(image,p1,p3,L1,L4);
    		x_is_F_of_y = True;
        else:
    		B = line_line_intersection_y(image,p1,p3,L1,L4);
    		x_is_F_of_y = False;
    			
    	# Step 3. If a linear is close enough, return coordinates.
        if (find_distance(A,B) <= TOL_LINEARS):
            vecCoord = [L1, L4]
            return vecCoord
    
        if TOL_QUADRATICS >= 0:
        	# QUADRATIC POLYNOMIAL
        	# Step 1. Build the quadratic polynomial
            if A.x==L4.x or A.x==L1.x:
        		pt = Coordinate(-1,-1);
        		vecCoord = []; 
        		vecCoord.append(pt);
        		return vecCoord;
        
            xx = [L4.x,A.x,L1.x];
            yy = [L4.y,A.y,L1.y];
            d = newt_coef(xx,yy); 
            if (x_is_F_of_y == False):
        		# Step 2. Search for the interface along a vertical line
        		p12 = Coordinate( (p1.x + L1.x)/2.0, p1.y);
        		p34 = Coordinate( p12.x, p4.y);
        		C = log_search(image,p12,p34);
        
        		# Step 3. Find intersection between quadratic and vertical line
        		D = Coordinate(p12.x, floor(newt_fct_eval(d,xx,p12.x)) );
            else:
        		# Step 2. Search for the interface along a vertical line
        		p14 = Coordinate( p1.x, (p1.y + L4.y) / 2.0 );
        		p23 = Coordinate( p2.x, p14.y );
        		C = log_search(image,p14,p23);
        
        		# Step 3. Find intersection between quadratic and vertical line
        		D = Coordinate( p14.y, floor(newt_fct_eval(d,xx,p14.y)) );
        
        		temp = D.x;
        		D.x = int(D.y);
        		D.y = int(temp);
           
            CD_dist = find_distance(C,D)
           
        	# Step 4. Check distance
            if (CD_dist <= TOL_QUADRATICS):
                vecCoord = [L4, L1, A];
                return vecCoord;
    
        if TOL_CUBICS >= 0:
        	# CUBIC POLYNOMIAL
            if (x_is_F_of_y == False):
        		# Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        		p1third12 = Coordinate( double ((L1.x - p1.x) / 3.0) + p1.x, p1.y);
        		p1third34 = Coordinate( p1third12.x, p4.y);
        		E = log_search(image,p1third12,p1third34);
        
        		p2thirds12 = Coordinate( double (2.0 * (L1.x - p1.x) / 3.0) + p1.x, p1.y);
        		p2thirds34 = Coordinate(p2thirds12.x, p4.y);
        		F = log_search(image,p2thirds12,p2thirds34);
        
        		# Step 2. Build the cubic polynomial
        		xi = [L4.x, E.x, F.x, L1.x];
        		yi = [L4.y, E.y, F.y, L1.y];
        		di = newt_coef(xi,yi);
        
        		# Step 3. Find intersection of interface with  this line
        		G = Coordinate( p1third12.x, floor(newt_fct_eval(di,xi,p1third12.x)) );
        		H = Coordinate( p2thirds12.x, floor(newt_fct_eval(di,xi,p2thirds12.x)) );
        
            else:
        		# Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        		p1third14 = Coordinate( p1.x, float ((L4.y - p1.y) / 3.0) + p1.y);
        		p1third23 = Coordinate( p2.x, p1third14.y);
        		E = log_search(image,p1third14,p1third23);
        
        		p2thirds14 = Coordinate( p1.x, float (2.0 * (L4.y - p1.y) / 3.0) + p1.y);
        		p2thirds23 = Coordinate( p2.x, p2thirds14.y);
        		F = log_search(image,p2thirds14,p2thirds23);
        
        		# Step 2. Build the cubic polynomial
        		xi = [L4.y, E.y, F.y, L1.y];
        		yi = [L4.x, E.x, F.x, L1.x];
        		di = newt_coef(xi,yi);
        
        		# Step 3. Find intersection of interface with  this line
        		G = Coordinate(p1third14.y, floor( newt_fct_eval(di,xi,p1third14.y)) );
        		H = Coordinate(p2thirds14.y, floor( newt_fct_eval(di,xi,p2thirds14.y)) );
        		
        		temp2 = G.x
        		G.x = int(G.y)
        		G.y = int(temp2)
        
        		temp3 = H.x
        		H.x = int(H.y)
        		H.y = int(temp3)
        
            EG_dist = find_distance(E,G)
            FH_dist = find_distance(F,H)
            
            # Step 4. Check distance
            if (EG_dist <= TOL_CUBICS) and (FH_dist <= TOL_CUBICS):       
        		vecCoord = [L4, L1, E, F];
        		return vecCoord;
    
        if poly_opt == 1:
            vec_list =  [L4,L1,A]
            return vec_list
        
        if poly_opt == 2:
            vec_list = [L4,L1,E,F] 
            return vec_list
        
    pt = Coordinate(-1,-1);
    vecCoord = [];
    vecCoord.append(pt);
    return vecCoord;


 
# case 2: 3:1 -- P2 the outsider
def case_NE_polynomial_test(image,p1,p2,p3,p4,L1,L2, poly_opt = 0):
    x_is_F_of_y = False;
		
    if ( find_distance(L1,L2) <= 2):
		pt = Coordinate(-1,-1);
		vecCoord = [];
		vecCoord.append(pt);
		return vecCoord;
		

	# LINEAR POLYNOMIAL
	# Step 1. Search for the interface along the 45 degree line
    A = log_search(image,p4,p2)
    if L2.x == p3.x and L2.y == p3.y:
		A = log_search(image,p1,p3)
 
	# Step 2. Find intersection of line between L1,L2 and the 45 degree line
    if (abs(p2.x-L1.x) < abs(L2.y - p2.y)):
		B = line_line_intersection_x(image,p4,p2,L1,L2);
		x_is_F_of_y = True;
    else:
		B = line_line_intersection_y(image,p4,p2,L1,L2);
		x_is_F_of_y = False;
	
 	# Step 3. If a linear is close enough, return coordinates.
    if(find_distance(A,B) <= TOL_LINEARS ):
# 		draw_line(imageOut,L1,L2);
		vecCoord = [L1, L2]; 
		return vecCoord;

    if TOL_QUADRATICS >= 0:
    	# QUADRATIC POLYNOMIAL
    	# Step 1. Build the quadratic polynomial
        xx = [L1.x, A.x, L2.x];
        yy = [L1.y, A.y, L2.y];
        if A.x==L2.x or A.x==L1.x:
    		pt = Coordinate(-1,-1);
    		vecCoord = [];
    		vecCoord.append(pt);
    		return vecCoord;
    
        d = newt_coef(xx,yy);
    
        if (x_is_F_of_y == False):
    		# Step 2. Search for the interface along a vertical line
    		p12 = Coordinate( (L1.x + p2.x)/2.0, p2.y );
    		p34 = Coordinate( p12.x, p3.y );
    		C = log_search(image, p12, p34);
    
    		# Step 3. Find intersection between quadratic and vertical line
    		D = Coordinate( p12.x, floor( newt_fct_eval(d,xx,p12.x)) );
    
        else:
    		# Step 2. Search for the interface along a horizontal line
    		p14 = Coordinate( p1.x, (p2.y + L2.y) / 2.0 );
    		p23 = Coordinate( p2.x, p14.y );
    		C = log_search(image, p14, p23);
    
    		# Step 3. Find intersection between quadratic and horizontal line
    		D = Coordinate( p14.y, floor( newt_fct_eval(d,xx,p14.y)) );
    
    		tempor = D.x;
    		D.x = D.y;
    		D.y = tempor;
       
        CD_dist = find_distance(C,D)
       
    	# Step 4. Check distance
        if (CD_dist <= TOL_QUADRATICS ):
            vecCoord = [L1, L2, A];
            return vecCoord; 

    if TOL_CUBICS >= 0:
    	# CUBIC POLYNOMIAL
        if( x_is_F_of_y == False):
    		# Step 1. Search for the interface along vertical lines at 1/3 and 2/3
    		p1third12 = Coordinate( float ((p2.x - L1.x) / 3.0) + L1.x, p1.y);
    		p1third34 = Coordinate( p1third12.x, p4.y );
    		E = log_search(image,p1third12,p1third34);
    
    		p2thirds12  = Coordinate( float (2.0 * (p2.x - L1.x) / 3.0) + L1.x, p1.y);
    		p2thirds34 = Coordinate( p2thirds12.x, p4.y);
    		F = log_search(image,p2thirds12,p2thirds34);
    
    		# Step 2. Build the cubic polynomial
    		xi = [L1.x, E.x, F.x, L2.x]; 
    		yi = [L1.y, E.y, F.y, L2.y];
    		di = newt_coef(xi,yi);
    	
    		# Step 3. Find intersection of interface with  this line
    		G = Coordinate( p1third12.x, floor( newt_fct_eval(di,xi,p1third12.x)) );	
    		H = Coordinate( p2thirds12.x, floor( newt_fct_eval(di,xi,p2thirds12.x)) );	
        else:
    		# Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
    		p1third14 = Coordinate( p1.x, float((L2.y - p2.y) / 3.0) + p2.y );
    		p1third23 = Coordinate( p2.x, p1third14.y);
    		E = log_search(image,p1third14,p1third23);
    
    		p2thirds14 = Coordinate( p1.x, float (2.0 * (L2.y - p2.y) / 3.0) + p2.y );
    		p2thirds23 = Coordinate( p2.x, p2thirds14.y );
    		F = log_search(image,p2thirds14,p2thirds23);
    
    		# Step 2. Build the cubic polynomial
    		xi = [L1.y, E.y, F.y, L2.y];
    		yi = [L1.x, E.x, F.x, L2.x];
    		di = newt_coef(xi,yi)
    
    		# Step 3. Find intersection of interface with  this line
    		G = Coordinate( p1third14.y, floor(newt_fct_eval(di,xi,p1third14.y)) );
    		H = Coordinate( p2thirds14.y, floor(newt_fct_eval(di,xi,p2thirds14.y)) );
    
    		tempor2 = G.x
    		G.x = int(G.y)
    		G.y = int(tempor2)
    
    		tempor3 = H.x
    		H.x = int(H.y)
    		H.y = int(tempor3)
    
        EG_dist = find_distance(E,G)
        FH_dist = find_distance(F,H)
        
        # Step 4. Check distance
        if (EG_dist <= TOL_CUBICS) and (FH_dist <= TOL_CUBICS):
       		vecCoord = [L1, L2, E, F];
    		return vecCoord;
    
    if poly_opt == 1:
        vec_list = [L1,L2,A]
        return vec_list
    
    if poly_opt == 2:
        vec_list = [L1,L2,E,F] 
        return vec_list
        
    pt = Coordinate(-1,-1);
    vecCoord = [];
    vecCoord.append(pt);
    return vecCoord;


## case 3: 3:1 - P3 is the outsider
def case_SE_polynomial_test(image,p1,p2,p3,p4,L2,L3, poly_opt = 0):
    is_x_F_of_y = False;

    if ( find_distance(L2,L3) <= 2):
		pt = Coordinate(-1,-1);
		vecCoord = [];
		vecCoord.append(pt);
		return vecCoord;
		

	# LINEAR POLYNOMIAL
	# Step 1. Search for the interface along the 45 degree line
    A = log_search(image,p1,p3);
	
	# Step 2. Find intersection of line between L2,L3 and the 45 degree line
    if( abs(L3.x - p3.x) < abs(p3.y - L2.y) ):
		B = line_line_intersection_x(image,p1,p3,L2,L3);
		is_x_F_of_y = True;

    else:
        B = line_line_intersection_y(image,p1,p3,L2,L3);
        is_x_F_of_y = False;

	# Step 3. If a linear is close enough, return coordinates.
    if(find_distance(A,B) <= TOL_LINEARS ):
		vecCoord = [L2, L3];
		return vecCoord;

    if TOL_QUADRATICS >= 0:
    	# QUADRATIC POLYNOMIAL
    	# Step 1. Build the quadratic polynomial
    	if A.x==L2.x or A.x==L3.x:
    		pt = Coordinate(-1,-1);
    		vecCoord = [];
    		vecCoord.append(pt);
    		return vecCoord;
    
        xx = [L3.x, L2.x, A.x];
        yy = [L3.y, L2.y, A.y];
        d = newt_coef(xx,yy);
    
        if (is_x_F_of_y == False): 
    		# Step 2. Search for the interface along a vertical line
    		p12 = Coordinate( (L3.x + p3.x)/2.0, p2.y);
    		p34 = Coordinate( p12.x, p3.y);
    		C = log_search(image,p12,p34);
    
    		# Step 3. Find intersection between quadratic and vertical line
    		D = Coordinate(p12.x, floor( newt_fct_eval(d,xx,p12.x)) );
        else:
    		# Step 2. Search for the interface along a vertical line
    		p14 = Coordinate( p4.x, (L2.y + p3.y) / 2.0 );
    		p23 = Coordinate( p3.x, p14.y);	
    		C = log_search(image,p14,p23);
    	
    		# Step 3. Find intersection between quadratic and horizontal line
    		D = Coordinate(p14.y, floor(newt_fct_eval(d,xx,p14.y)) );
    	
    		tmp = D.x;
    		D.x = D.y;
    		D.y = tmp;
    
        CD_dist = find_distance(C,D)
        
    	# Step 4. Check distance
        if (CD_dist <= TOL_QUADRATICS):
            vecCoord = [L3, L2, A];
            return vecCoord;

    if TOL_CUBICS >= 0:
    	# CUBIC POLYNOMIAL
        if (is_x_F_of_y == False):
    		# Step 1. Search for the interface along vertical lines at 1/3 and 2/3
    		p1third12 = Coordinate( float ((p3.x - L3.x) / 3.0) + L3.x, p1.y);
    		p1third34 = Coordinate( p1third12.x, p4.y);
    		E = log_search(image,p1third12,p1third34);
    
    		p2thirds12  = Coordinate( float (2.0 * (p3.x - L3.x) / 3.0) + L3.x, p1.y );
    		p2thirds34 = Coordinate( p2thirds12.x, p4.y );
    		F = log_search(image,p2thirds12,p2thirds34);
    
    		# Step 2. Build the cubic polynomial
    		xi = [L3.x, E.x, F.x, L2.x];
    		yi = [L3.y, E.y, F.y, L2.y];
    		di = newt_coef(xi,yi);
    	
    		# Step 3. Find intersection of interface with  this line
    		G = Coordinate( p1third12.x, floor( newt_fct_eval(di,xi,p1third12.x)) );
    		H = Coordinate( p2thirds12.x, floor( newt_fct_eval(di,xi,p2thirds12.x)) );
        else:
    		# Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
    		p1third14 = Coordinate( p1.x, float ((p3.y - L2.y) / 3.0) + L2.y );
    		p1third23 = Coordinate( p2.x, p1third14.y);
    		E = log_search(image,p1third14,p1third23);
    
    		p2thirds14 = Coordinate( p1.x, float (2.0 * (p3.y - L2.y) / 3.0) + L2.y );
    		p2thirds23 = Coordinate( p2.x, p2thirds14.y );
    		F = log_search(image,p2thirds14,p2thirds23);
    
    		# Step 2. Build the cubic polynomial
    		xi = [ L3.y, E.y, F.y, L2.y];
    		yi = [ L3.x, E.x, F.x, L2.x];
    		di = newt_coef(xi,yi);
    
    		G = Coordinate(p1third14.y, floor(newt_fct_eval(di,xi,p1third14.y)) );
    		H = Coordinate(p2thirds14.y, floor(newt_fct_eval(di,xi,p2thirds14.y)) );
    
    		tmp2 = G.x
    		G.x = int(G.y)
    		G.y = int(tmp)
    
    		tmp3 = H.x
    		H.x = int(H.y)
    		H.y = int(tmp3)
    
        EG_dist = find_distance(E,G)
        FH_dist = find_distance(F,H)
        
        # Step 4. Check distance
        if (EG_dist <= TOL_CUBICS) and (FH_dist <= TOL_CUBICS):
     		vecCoord = [L3, L2, E, F]
    		return vecCoord;
    

    if poly_opt == 1:
        vec_list = [L3,L2,A]
        return vec_list
    
    if poly_opt == 2:
        vec_list = [L3,L2,E,F] 
        return vec_list
            
    pt = Coordinate(-1,-1);
    vecCoord = [];
    vecCoord.append(pt);
    return vecCoord;


# case 4: 3:1 -- P4 the outsider
def case_SW_polynomial_test(image, p1, p2, p3, p4,L3,L4, poly_opt = 0):
    is_x_F_of_y = False;

    if ( find_distance(L4,L3) <= 2):
		pt = Coordinate(-1,-1);
		vecCoord = [];
		vecCoord.append(pt);
		return vecCoord;


	# LINEAR POLYNOMIAL
	# Step 1. Search for the interface along the 45 degree line
    A = log_search(image,p4,p2);
	
	# Step 2. Find intersection of line between L4,L3 and the 45 degree line
    if ( abs(L3.x - p4.x) < abs(p4.y - L4.y) ):
		B = line_line_intersection_x(image,p4,p2,L4,L3);
		is_x_F_of_y = True;
    else: 
		B = line_line_intersection_y(image,p4,p2,L4,L3);
		is_x_F_of_y = False;

	# Step 3. If a linear is close enough, return coordinates.
    if(find_distance(A,B) <= TOL_LINEARS ):
		vecCoord = [L4, L3];
		return vecCoord;

    if TOL_QUADRATICS >= 0:
    	# QUADRATIC POLYNOMIAL
    	# Step 1. Build the quadratic polynomial
        if A.x==L4.x or A.x==L3.x:
    		pt = Coordinate(-1,-1);
    		vecCoord = []; 
    		vecCoord.append(pt);
    		return vecCoord;
    
        xx = [L4.x, A.x, L3.x];
        yy = [L4.y, A.y, L3.y];
        d = newt_coef(xx,yy);
    
        if (is_x_F_of_y == False):
    		# Step 2. Search for the interface along a vertical line
    		p12 = Coordinate( (p4.x + L3.x)/2.0, p1.y);
    		p34 = Coordinate( p12.x, p4.y );
    		C = log_search(image,p12,p34);
    
    		# Step 3. Find intersection between quadratic and vertical line
    		D = Coordinate( p12.x, floor( newt_fct_eval(d,xx,p12.x)) );
        else:
    		# Step 2. Search for the interface along a vertical line
    		p14 = Coordinate( p1.x, (L4.y + p4.y) / 2.0) ;
    		p23 = Coordinate(p3.x, p14.y);
    		C = log_search(image,p14,p23);
    		 
    		# Step 3. Find intersection between quadratic and vertical line
    		D = Coordinate( p14.y, floor( newt_fct_eval(d,xx,p14.y)) );
    		
    		tmp = D.x;
    		D.x = D.y;	
    		D.y = tmp;
    
        CD_dist = find_distance(C,D)
        if (CD_dist <= TOL_QUADRATICS):
            vecCoord = [L4, L3, A];
            return vecCoord;
    
    if TOL_CUBICS >= 0:
    	# CUBIC POLYNOMIAL
        if (is_x_F_of_y == False):
            # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
            p1third12 = Coordinate( float ((p4.x - L3.x) / 3.0) + L3.x, p1.y);
            p1third34 = Coordinate( p1third12.x, p4.y);
            E = log_search(image,p1third12,p1third34);
            
            p2thirds12 = Coordinate(float (2.0 * (p4.x - L3.x) / 3.0) + L3.x, p1.y);
            p2thirds34 = Coordinate( p2thirds12.x, p4.y);
            F = log_search(image, p2thirds12, p2thirds34);
		
            # Step 2. Build the cubic polynomial
            xi = [L4.x, E.x, F.x, L3.x];		
            yi = [L4.y, E.y, F.y, L3.y];		
            di = newt_coef(xi,yi);
            
            # Step 3. Find intersection of interface with  this line
            G = Coordinate(p1third12.x, floor( newt_fct_eval(di,xi,p1third12.x)) );
            H = Coordinate(p2thirds12.x, floor( newt_fct_eval(di,xi,p2thirds12.x)) );
        else:
    		# Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
    		p1third14 = Coordinate( p1.x, float ((p4.y - L4.y) / 3.0) + L4.y);
    		p1third23 = Coordinate( p2.x, p1third14.y);
    		E = log_search(image,p1third14,p1third23);
    
    		p2thirds14 = Coordinate(p1.x, float (2.0 * (p4.y - L4.y) / 3.0) + L4.y );
    		p2thirds23 = Coordinate(p2.x, p2thirds14.y);
    		F = log_search(image, p2thirds14,p2thirds23);
    
    		# Step 2. Build the cubic polynomial
    		xi = [L4.y, E.y, F.y, L3.y];
    		yi = [L4.x, E.x, F.x, L3.x];
    		di = newt_coef(xi,yi);
    
    		# Step 3. Find intersection of interface with  this line
    		G = Coordinate(p1third14.y,floor( newt_fct_eval(di,xi,p1third14.y)) );
    		H = Coordinate(p2thirds14.y,floor( newt_fct_eval(di,xi,p2thirds14.y)) );
    
    		tmp2 = G.x
    		G.x = int(G.y)
    		G.y = int(tmp2)
    
    		tmp3 = H.x
    		H.x = int(H.y)
    		H.y = int(tmp3)
    
        EG_dist = find_distance(E,G)
        FH_dist = find_distance(F,H)
        
        # Step 4. Check distance
        if (EG_dist <= TOL_CUBICS) and (FH_dist <= TOL_CUBICS):
    		vecCoord = [L4, L3, E, F];
    		return vecCoord;

    if poly_opt == 1:
        vec_list = [L4,L3,A]
        return vec_list
    
    if poly_opt == 2:
        vec_list = [L4,L3,E,F] 
        return vec_list
            
    pt = Coordinate(-1,-1);
    vecCoord = [];
    vecCoord.append(pt);
    return vecCoord;
		

# case 5: 2:2 vertical crossing
def case_vertical_polynomial_test(image, p1, p2, p3, p4,L1,L3, poly_opt = 0):

    if ( find_distance(L1,L3) <= 2):
        pt = Coordinate(-1,-1);
        vecCoord = [];
        vecCoord.append(pt);
        return vecCoord;

    p14mid = find_mid_point(p1,p4)
    p23mid = find_mid_point(p2,p3)
	
	# LINEAR POLYNOMIAL
	# Step 1. Search for the interface along the mid line between p1 and p2
    A = log_search(image, p14mid, p23mid)
	
	# Step 2. Find intersection of line between L1,L3 and the 45 degree line
    B = line_line_intersection_x(image,p1,p3,L1,L3);

	# Step 3. If a linear is close enough, return coordinates.
    if( find_distance(A,B) <= TOL_LINEARS):
		vecCoord = [L1, L3];
		return vecCoord;

    if TOL_QUADRATICS >= 0:
    	# QUADRATIC POLYNOMIAL
    	# Step 1. Build the quadratic polynomial
        xx = [L1.y, A.y, L3.y];
        yy = [L1.x, A.x, L3.x];
        if A.x==L1.x or A.x==L3.x:
    		pt = Coordinate(-1,-1);
    		vecCoord = []; 
    		vecCoord.append(pt);
    		return vecCoord;
    
        d = newt_coef(xx,yy);
    	
    	#Step 2. Search for the interface along an horizontal line
        p14 = Coordinate(p1.x, floor( (p1.y + p4.y)/2.0) );
        p23 = Coordinate(p2.x, floor( (p2.y + p3.y)/2.0) );
        C = log_search(image, p14, p23);
    
    	# Step 3. Find intersection between quadratic and horizontal line
        D = Coordinate(p14.y, floor(newt_fct_eval(d,xx,p14.y)) );
    	
        tempor = D.x;
        D.x = D.y;
        D.y = tempor;
        
        CD_dist = find_distance(C,D)
        
        if( CD_dist <= TOL_QUADRATICS):
            vecCoord = [L1, L3, A];
            return vecCoord;
    
    if TOL_CUBICS >= 0:
    	# CUBIC POLYNOMIAL
    	# Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third14 = Coordinate( p1.x, float ((p4.y - p1.y) / 3.0) + p1.y);
        p1third23 = Coordinate( p2.x, float ((p3.y - p2.y) / 3.0) + p2.y);
        E = log_search(image,p1third14,p1third23);
    
        p2thirds14 = Coordinate(p1.x, float (2.0 * (p4.y - p1.y) / 3.0) + p1.y );
        p2thirds23 = Coordinate(p2.x, float (2.0 * (p3.y - p2.y) / 3.0) + p2.y );
        F = log_search(image,p2thirds14,p2thirds23);
    
    	# Step 2. Build the cubic polynomial
        xi = [L1.y, E.y, F.y, L3.y];
        yi = [L1.x, E.x, F.x, L3.x];
        di = newt_coef(xi,yi);
    
    	# Step 3. Find intersection of interface with  this line
        G = Coordinate(p1third14.y, floor( newt_fct_eval(di,xi,p1third14.y)) );
        H = Coordinate(p2thirds14.y, floor( newt_fct_eval(di,xi,p2thirds14.y)) );
    
        tmp = G.x
        G.x = int(G.y)
        G.y = int(tmp)
    
        tmp2 = H.x
        H.x = int(H.y)
        H.y = int(tmp2)
    
        EG_dist = find_distance(E,G)
        FH_dist = find_distance(H,F)
    	# Step 4. Check distance
        if(EG_dist <= TOL_CUBICS) and (FH_dist <= TOL_CUBICS):
    		vecCoord = [L1, L3, E, F]
    		return vecCoord;	
    
    if poly_opt == 1:
        vec_list = [L1,L3,A]
        return vec_list
    
    if poly_opt == 2:
        vec_list = [L1,L3,E,F]
        return vec_list
    
    pt = Coordinate(-1,-1);
    vecCoord = [];
    vecCoord.append(pt);
    return vecCoord;

# case 6: 2:2 horizontal crossing
def case_horizontal_polynomial_test(image,p1,p2,p3,p4,L2,L4, poly_opt = 0):
	
    if ( find_distance(L2,L4) <= 2):
		pt = Coordinate(-1,-1);
		vecCoord = [];
		vecCoord.append(pt);
		return vecCoord;

    p12mid = find_mid_point(p1,p2)
    p34mid = find_mid_point(p4,p3)
	
	#LINEAR POLYNOMIAL
	# Step 1. Search for the interface along the mid line between p1 and p2
    A = log_search(image, p12mid, p34mid)
    if L4.x == A.x and L4.y == A.y:
		A = log_search(image,p1,p3)

	# Step 2. Find intersection of line between L2, L4 and the 45 degrees line.
    B = line_line_intersection_y(image, p2, p4, L2, L4);
    
    AB_dist = find_distance(A,B)
	
    # Step 3. If a linear is close enough, return coordinates.
    if (AB_dist <= TOL_LINEARS):
		vecCoord = [L2, L4];
		return vecCoord;

    if TOL_QUADRATICS >= 0:
    	# QUADRATIC POLYNOMIAL
    	# Step 1. Build the quadratic polynomial
        xx = [L4.x, A.x, L2.x];
        yy = [L4.y, A.y, L2.y];
    
        d = newt_coef(xx,yy);
    
    	# Step 2. Search for the interface along a vertical line
        p12 = Coordinate( (p1.x + p2.x)/2.0, p1.y);
        p34 = Coordinate( (p3.x + p4.x)/2.0, p4.y);
        C = log_search(image,p12,p34);
    
    	# Step 3. Find intersection between quadratic and vertical line
        D = Coordinate(p12.x, floor( newt_fct_eval(d,xx,p12.x)) );
    	
        CD_dist = find_distance(C,D)
         
    	# Step 4. Check distance
        if(CD_dist <= TOL_QUADRATICS):
    		vecCoord = [L2, L4, A];
    		return vecCoord;
    
    if TOL_CUBICS >= 0:    
    	# CUBIC POLYNOMIAL
    	# Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate( float ((p2.x - p1.x) / 3.0) + p1.x, p1.y);
        p1third34 = Coordinate( float ((p3.x - p4.x) / 3.0) + p4.x, p4.y);
        E = log_search(image,p1third12,p1third34);
    
        p2thirds12 = Coordinate( float (2.0 * (p2.x - p1.x) / 3.0) + p1.x, p1.y);
        p2thirds34 = Coordinate( float (2.0 * (p3.x - p4.x) / 3.0) + p4.x, p4.y);
        F = log_search(image,p2thirds12,p2thirds34);
    
    	# Step 2. Build the cubic polynomial
        xi = [L4.x, E.x, F.x, L2.x];
        yi = [L4.y, E.y, F.y, L2.y];
        di = newt_coef(xi,yi);
    
    	# Step 3. Find intersection of interface with  this line
        G = Coordinate( p1third12.x, floor( newt_fct_eval(di,xi,p1third12.x)) );
        H = Coordinate( p2thirds12.x, floor( newt_fct_eval(di,xi,p2thirds12.x)) );
    
        EG_dist = find_distance(E,G)
        FH_dist = find_distance(F,H)
    	# Step 4. Check distance
        if (EG_dist <= TOL_CUBICS) and (FH_dist <= TOL_CUBICS):
    		vecCoord = [L2, L4, E, F]
    		return vecCoord;

    if poly_opt == 1:
        vec_list = [L2,L4,A]
        return vec_list
    
    if poly_opt == 2:
        vec_list = [L2,L4,E,F]
        return vec_list

        
    pt = Coordinate(-1,-1)
    vecCoord = []
    vecCoord.append(pt)
    return vecCoord

def Nurbs_control_points(Q):
    R = numpy.matrix([[1.0, 0.0, 0.0, 0.0],[ 1.0/9.0, 2.0/3.0, 2.0/9.0, 0.0 ], [  0.0, 2.0/9.0, 2.0/3.0, 1.0/9.0],[ 0.0, 0.0, 0.0, 1.0] ] )
    Q = numpy.matrix( [ [Q[0][0], Q[0][1]], [Q[1][0], Q[1][1]], [Q[2][0], Q[2][1]], [Q[3][0], Q[3][1]] ])
    
    P = numpy.linalg.solve(R,Q)
    return P

def Nurbs_basis_fcts(t,P):

    if t < 1.0/2.0:
        N = lambda t: P[0] * (1 - 2 * t) * (1 - 2 * t) + P[1] * 2 * t * (2 - 3 * t) + P[2] * 2 * t * t
    else:
        N = lambda t: P[1] * 2 * (1 - t) * (1 - t) + P[2] * (8*t - 6*t*t - 2) + P[3] * (2*t - 1) * (2*t -1) 
    return N

def Nurbs_NW_case(image,p1,p2,p3,p4,L1,L4):
    
    x_is_F_of_y = False
    
    if (abs(L1.x-p1.x) < abs(L4.y - p1.y)):
        x_is_F_of_y = True
    else:
        x_is_F_of_y = False
    
    if (x_is_F_of_y == False):
        
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate( double ((L1.x - p1.x) / 3.0) + p1.x, p1.y);
        p1third34 = Coordinate( p1third12.x, p4.y);
        E = log_search(image,p1third12,p1third34);
	
        p2thirds12 = Coordinate( double (2.0 * (L1.x - p1.x) / 3.0) + p1.x, p1.y);
        p2thirds34 = Coordinate(p2thirds12.x, p4.y);
        F = log_search(image,p2thirds12,p2thirds34);
        
        # Step 2. Build the sample points vector   
        Q = [ [L4.x, L4.y], [E.x, E.y], [F.x, F.y], [L1.x, L1.y]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(L1.x-p1.x)
        
        # Step 4. Search for the interface along a vertical line
        p12 = Coordinate( (p1.x + L1.x)/2.0, p1.y)
        p34 = Coordinate( p12.x, p4.y)
        C = log_search(image,p12,p34)
         
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
               
    else:
        
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate( p1.x, float ((L4.y - p1.y) / 3.0) + p1.y);
        p1third23 = Coordinate( p2.x, p1third14.y);
        E = log_search(image,p1third14,p1third23);
        
        p2thirds14 = Coordinate( p1.x, float (2.0 * (L4.y - p1.y) / 3.0) + p1.y);
        p2thirds23 = Coordinate( p2.x, p2thirds14.y);
        F = log_search(image,p2thirds14,p2thirds23);

        # Step 2. Build the sample points vector   
        Q = [ [L4.y, L4.x], [F.y, F.x], [E.y, E.x], [L1.y, L1.x]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)

        t_step = 1.0 / abs(L4.y - p1.y)
        
        # Step 4. Search for the interface along a vertical line
        p14 = Coordinate( p1.x, (p1.y + L4.y) / 2.0 )
        p23 = Coordinate( p2.x, p14.y )
        C = log_search(image,p14,p23)
        
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
        
    t = arange(0, 1+t_step, t_step)
    
    if find_distance(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]


def Nurbs_NE_case(image,p1,p2,p3,p4,L1,L2):
    
    x_is_F_of_y = False
        
    if (abs(p2.x-L1.x) < abs(L2.y - p2.y)):
        x_is_F_of_y = True;
    else:
        x_is_F_of_y = False;

    if( x_is_F_of_y == False):
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate( float ((p2.x - L1.x) / 3.0) + L1.x, p1.y);
        p1third34 = Coordinate( p1third12.x, p4.y );
        E = log_search(image,p1third12,p1third34);

        p2thirds12  = Coordinate( float (2.0 * (p2.x - L1.x) / 3.0) + L1.x, p1.y);
        p2thirds34 = Coordinate( p2thirds12.x, p4.y);
        F = log_search(image,p2thirds12,p2thirds34);
        
        # Step 2. Build the sample points vector   
        Q = [ [L1.x, L1.y], [E.x, E.y], [F.x, F.y], [L2.x, L2.y]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(p2.x - L1.x)
        
        p12 = Coordinate( (L1.x + p2.x)/2.0, p2.y )
        p34 = Coordinate( p12.x, p3.y )
        C = log_search(image, p12, p34)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
        
    else:
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate( p1.x, float((L2.y - p2.y) / 3.0) + p2.y );
        p1third23 = Coordinate( p2.x, p1third14.y);
        E = log_search(image,p1third14,p1third23);

        p2thirds14 = Coordinate( p1.x, float (2.0 * (L2.y - p2.y) / 3.0) + p2.y );
        p2thirds23 = Coordinate( p2.x, p2thirds14.y );
        F = log_search(image,p2thirds14,p2thirds23);
        
        # Step 2. Build the sample points vector   
        Q = [ [L1.y, L1.x], [E.y, E.x], [F.y, F.x], [L2.y, L2.x]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(L2.y - p2.y)
        
        p14 = Coordinate( p1.x, (p2.y + L2.y) / 2.0 )
        p23 = Coordinate( p2.x, p14.y )
        C = log_search(image, p14, p23)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
        
    t = arange(0, 1+t_step, t_step)
    
    if find_distance(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
    
def Nurbs_SE_case(image,p1,p2,p3,p4,L2,L3):
    
    x_is_F_of_y = False
    if( abs(L3.x - p3.x) < abs(p3.y - L2.y) ):
        x_is_F_of_y = True;

    else:
        x_is_F_of_y = False;

    if (x_is_F_of_y == False):
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate( float ((p3.x - L3.x) / 3.0) + L3.x, p1.y);
        p1third34 = Coordinate( p1third12.x, p4.y);
        E = log_search(image,p1third12,p1third34);

        p2thirds12  = Coordinate( float (2.0 * (p3.x - L3.x) / 3.0) + L3.x, p1.y );
        p2thirds34 = Coordinate( p2thirds12.x, p4.y );
        F = log_search(image,p2thirds12,p2thirds34);

        # Step 2. Build the sample points vector   
        Q = [ [L3.x, L3.y], [E.x, E.y], [F.x, F.y], [L2.x, L2.y]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(p3.x - L3.x)
        
        p12 = Coordinate( (L3.x + p3.x)/2.0, p2.y)
        p34 = Coordinate( p12.x, p3.y)
        C = log_search(image,p12,p34)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
        
    else:
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate( p1.x, float ((p3.y - L2.y) / 3.0) + L2.y );
        p1third23 = Coordinate( p2.x, p1third14.y);
        E = log_search(image,p1third14,p1third23);

        p2thirds14 = Coordinate( p1.x, float (2.0 * (p3.y - L2.y) / 3.0) + L2.y );
        p2thirds23 = Coordinate( p2.x, p2thirds14.y );
        F = log_search(image,p2thirds14,p2thirds23);

        # Step 2. Build the sample points vector   
        Q = [ [L3.y, L3.x], [F.y, F.x], [E.y, E.x], [L2.y, L2.x]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
    
        t_step = 1.0 / abs(p3.y - L2.y)
        
        p14 = Coordinate( p4.x, (L2.y + p3.y) / 2.0 )
        p23 = Coordinate( p3.x, p14.y)    
        C = log_search(image,p14,p23)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
    
    t = arange(0, 1+t_step, t_step)
    
    if find_distance(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
                    
def Nurbs_SW_case(image,p1,p2,p3,p4,L3,L4):
    
    x_is_F_of_y = False
    
    # Step 2. Find intersection of line between L4,L3 and the 45 degree line
    if ( abs(L3.x - p4.x) < abs(p4.y - L4.y) ):
        x_is_F_of_y = True;
    else: 
        x_is_F_of_y = False;

    if (x_is_F_of_y == False):
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate( float ((p4.x - L3.x) / 3.0) + L3.x, p1.y);
        p1third34 = Coordinate( p1third12.x, p4.y);
        E = log_search(image,p1third12,p1third34);
        
        p2thirds12 = Coordinate(float (2.0 * (p4.x - L3.x) / 3.0) + L3.x, p1.y);
        p2thirds34 = Coordinate( p2thirds12.x, p4.y);
        F = log_search(image, p2thirds12, p2thirds34);

        # Step 2. Build the sample points vector   
        Q = [ [L4.x, L4.y], [F.x, F.y], [E.x, E.y], [L3.x, L3.y]]

        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(p4.x - L3.x)
        
        p12 = Coordinate( (p4.x + L3.x)/2.0, p1.y)
        p34 = Coordinate( p12.x, p4.y )
        C = log_search(image,p12,p34)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
        
    else:
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate( p1.x, float ((p4.y - L4.y) / 3.0) + L4.y);
        p1third23 = Coordinate( p2.x, p1third14.y);
        E = log_search(image,p1third14,p1third23);

        p2thirds14 = Coordinate(p1.x, float (2.0 * (p4.y - L4.y) / 3.0) + L4.y );
        p2thirds23 = Coordinate(p2.x, p2thirds14.y);
        F = log_search(image, p2thirds14,p2thirds23);

        # Step 2. Build the sample points vector   
        Q = [ [L3.y, L3.x], [F.y, F.x], [E.y, E.x], [L4.y, L4.x]]
    
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
    
        t_step = 1.0 / abs(p4.y - L4.y)
        
        p14 = Coordinate( p1.x, (L4.y + p4.y) / 2.0)
        p23 = Coordinate(p3.x, p14.y)
        C = log_search(image,p14,p23)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )

    
    t = arange(0, 1+t_step, t_step)
    
    if find_distance(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
    
def Nurbs_vertical_case(image,p1,p2,p3,p4,L1,L3):
    
    x_is_F_of_y = True
    
    # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
    p1third14 = Coordinate( p1.x, float ((p4.y - p1.y) / 3.0) + p1.y);
    p1third23 = Coordinate( p2.x, float ((p3.y - p2.y) / 3.0) + p2.y);
    E = log_search(image,p1third14,p1third23);

    p2thirds14 = Coordinate(p1.x, float (2.0 * (p4.y - p1.y) / 3.0) + p1.y );
    p2thirds23 = Coordinate(p2.x, float (2.0 * (p3.y - p2.y) / 3.0) + p2.y );
    F = log_search(image,p2thirds14,p2thirds23);

    # Step 2. Build the sample points vector   
    Q = [ [L3.y, L3.x], [F.y, F.x], [E.y, E.x], [L1.y, L1.x]]
    
    # Step 3. Determine the control points
    P = Nurbs_control_points(Q)

    t_step = 1.0 / abs(p1.y - p4.y)
    
    p14 = Coordinate(p1.x, floor( (p1.y + p4.y)/2.0) )
    p23 = Coordinate(p2.x, floor( (p2.y + p3.y)/2.0) )
    C = log_search(image, p14, p23)
        
    Nx = Nurbs_basis_fcts(0.5,P[:,1])
    Ny = Nurbs_basis_fcts(0.5,P[:,0])
    NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
        
    t = arange(0, 1+t_step, t_step)
    
    if find_distance(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
                        
def Nurbs_horizontal_case(image,p1,p2,p3,p4,L2,L4):
    
    x_is_F_of_y = False    

    # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
    p1third12 = Coordinate( float ((p2.x - p1.x) / 3.0) + p1.x, p1.y);
    p1third34 = Coordinate( float ((p3.x - p4.x) / 3.0) + p4.x, p4.y);
    E = log_search(image,p1third12,p1third34);

    p2thirds12 = Coordinate( float (2.0 * (p2.x - p1.x) / 3.0) + p1.x, p1.y);
    p2thirds34 = Coordinate( float (2.0 * (p3.x - p4.x) / 3.0) + p4.x, p4.y);
    F = log_search(image,p2thirds12,p2thirds34);

    # Step 2. Build the cubic polynomial
    xi = [L4.x, E.x, F.x, L2.x];
    yi = [L4.y, E.y, F.y, L2.y];
        
            
    # Step 2. Build the sample points vector   
    Q = [ [L4.x, L4.y], [E.x, E.y], [F.x, F.y], [L2.x, L2.y]]
    
    # Step 3. Determine the control points
    P = Nurbs_control_points(Q)

    t_step = 1.0 / abs(p1.x - p2.x)
    
    p12 = Coordinate( (p1.x + p2.x)/2.0, p1.y)
    p34 = Coordinate( (p3.x + p4.x)/2.0, p4.y)
    C = log_search(image,p12,p34)
        
    Nx = Nurbs_basis_fcts(0.5,P[:,0])
    Ny = Nurbs_basis_fcts(0.5,P[:,1])
    NC = Coordinate( int(Nx(0.5)) , int(Ny(0.5)) )
            
    t = arange(0, 1+t_step, t_step)
    
    if find_distance(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
                    
def draw_nurbs(image,t,P,x_is_F_of_y,p1,p2,p4):
    px_col = 255
    
    imageSize = image.GetSize()
    
    Px = P[:,0]
    Py = P[:,1]
    
    for i in range(0,len(t)):
        if x_is_F_of_y == False:
            Nx = Nurbs_basis_fcts(t[i],Px)
            Ny = Nurbs_basis_fcts(t[i],Py)
        else:
            Nx = Nurbs_basis_fcts(t[i],Py)
            Ny = Nurbs_basis_fcts(t[i],Px)

        xloc = int(Nx(t[i]))
        yloc = int(Ny(t[i]))
        
        if ( (p1.x <= xloc and xloc <= p2.x and p1.y <= yloc and yloc <= p4.y) or 
         (p1.x <= yloc and yloc <= p2.x and p1.y <= xloc and xloc <= p4.y) ):
                image.SetPixel(xloc,yloc,0,px_col)
        
    
