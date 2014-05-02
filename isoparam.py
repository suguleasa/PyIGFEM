import numpy
import sympy
import scipy

def iso_basis_fct():
	# defining basis functions for an isoparametric element: square centered at (0,0)
	# sides are of length 2
	# x,y get values from (-1,-1) to (1,1)

	# the south-west corner
	N1 = lambda e,n: 1.0/4.0 * (1 - e) * (1 - n) 
	# the south-east corner
	N2 = lambda e,n: 1.0/4.0 * (1 + e) * (1 - n)
	# the north-east corner
	N3 = lambda e,n: 1.0/4.0 * (1 + e) * (1 + n)
	# the north-west corner
	N4 = lambda e,n: 1.0/4.0 * (1 - e) * (1 + n)

	
	return [N1,N2,N3,N4]

def iso_basis_fct_HN(N,S,E,W):
	# for hanging nodes
	# defining basis functions for an isoparametric element: square centered at (0,0)
	# sides are of length 2
	# x,y get values from (-1,-1) to (1,1)

	# the South edge:
	if S == 1:
		N5 = lambda e,n: 1.0/2.0 * (1 - abs(e)) * (1 - n)
	else:
		N5 = lambda e,n: 0.0
		
	#the East edge:
	if E == 1:
		N6 = lambda e,n: 1.0/2.0 * (1 + e) * (1 - abs(n))
	else:
		N6 = lambda e,n: 0.0
		
	# the North edge:
	if N == 1:
		N7 = lambda e,n: 1.0/2.0 * (1 - abs(e)) * ( 1 + n)
	else:
		N7 = lambda e,n: 0.0
		
	# the West edge:
	if W == 1:
		N8 = lambda e,n: 1.0/2.0 * (1 - e) * (1 - abs(n))
	else:
		N8 = lambda e,n: 0.0
	
	# the south-west corner
	N1 = lambda e,n: 1.0/4.0 * (1 - e) * (1 - n) - 1.0/2.0 * (N5(e,n) + N8(e,n)) 
	# the south-east corner
	N2 = lambda e,n: 1.0/4.0 * (1 + e) * (1 - n) - 1.0/2.0 * (N5(e,n) + N6(e,n))
	# the north-east corner
	N3 = lambda e,n: 1.0/4.0 * (1 + e) * (1 + n) - 1.0/2.0 * (N6(e,n) + N7(e,n))
	# the north-west corner
	N4 = lambda e,n: 1.0/4.0 * (1 - e) * (1 + n) - 1.0/2.0 * (N7(e,n) + N8(e,n))

	
	return [N1,N2,N3,N4]

def tri_iso_basis_fct():
	# defining basis functions for an isoparametri TRIANGLE
	N1 = lambda e,n: (1 - e - n)
	N2 = lambda e,n: e 
	N3 = lambda e,n: n

	return [N1,N2,N3]





def iso_deriv_wrt_to_e():
	# defining the derivatives with respect to eta of the isoparametric basis functions

	# the south-west corner
	N1e = lambda e,n: 1.0/4.0 * (-1) * (1 - n)
	# the south-east corner
	N2e = lambda e,n: 1.0/4.0 * (1 - n)
	# the north-east corner
	N3e = lambda e,n: 1.0/4.0 * (1 + n)
	# the north-west corner
	N4e = lambda e,n: 1.0/4.0 * (-1) * (1 + n)

	
	return [N1e,N2e,N3e,N4e]

def iso_deriv_wrt_to_n():
	# defining the derivatives with respect to niu of the isoparametric basis functions

	# the south-west corner
	N1n = lambda e,n: 1.0/4.0 * (1 - e) * (-1) 
	# the south-east corner
	N2n = lambda e,n: 1.0/4.0 * (1 + e) * (-1)
	# the north-east corner
	N3n = lambda e,n: 1.0/4.0 * (1 + e)
	# the north-west corner
	N4n = lambda e,n: 1.0/4.0 * (1 - e)

	return [N1n,N2n,N3n,N4n]

def iso_deriv_wrt_to_e_HN(N,S,E,W):
	# with hanging nodes
	# defining the derivatives with respect to eta of the isoparametric basis functions
	
	# the South edge:
	if S == 1:
		N5e = lambda e,n: 1.0/2.0 * (-1) * (1 - n)
	else:
		N5e = lambda e,n: 0.0
		
	#the East edge:
	if E == 1:		
		N6e = lambda e,n: 1.0/2.0 * (1 - abs(n))
	else:
		N6e = lambda e,n: 0.0
		
	# the North edge:
	if N == 1:
		N7e = lambda e,n: 1.0/2.0 * (-1) * ( 1 + n)
	else:
		N7e = lambda e,n: 0.0
		
	# the West edge:
	if W == 1:
		N8e = lambda e,n: 1.0/2.0 * (-1) * (1 - abs(n))
	else:
		N8e = lambda e,n: 0.0

	# the south-west corner
	N1e = lambda e,n: 1.0/4.0 * (-1) * (1 - n) - 1.0/2.0 * (N5e(e,n) + N8e(e,n))
	# the south-east corner
	N2e = lambda e,n: 1.0/4.0 * (1 - n) - 1.0/2.0 * (N5e(e,n) + N6e(e,n))
	# the north-east corner
	N3e = lambda e,n: 1.0/4.0 * (1 + n) - 1.0/2.0 * (N6e(e,n) + N7e(e,n))
	# the north-west corner
	N4e = lambda e,n: 1.0/4.0 * (-1) * (1 + n) - 1.0/2.0 * (N7e(e,n) + N8e(e,n))
	
	return [N1e,N2e,N3e,N4e]

def iso_deriv_wrt_to_n_HN(N,S,E,W):
	# with hanging nodes
	# defining the derivatives with respect to niu of the isoparametric basis functions

	# the South edge:
	if S == 1:
		N5n = lambda e,n: 1.0/2.0 * (1 - abs(e)) * (-1)
	else:
		N5n = lambda e,n: 0.0
		
	#the East edge:
	if E == 1:
		N6n = lambda e,n: 1.0/2.0 * (1 + e) * (-1)
	else:
		N6n = lambda e,n: 0.0
		
	# the North edge:
	if N == 1:
		N7n = lambda e,n: 1.0/2.0 * (1 - abs(e)) * ( 1)
	else:
		N7n = lambda e,n: 0.0
		
	# the West edge:
	if W == 1:
		N8n = lambda e,n: 1.0/2.0 * (1 - e) * (-1)
	else:
		N8n = lambda e,n: 0.0	
	
	# the south-west corner
	N1n = lambda e,n: 1.0/4.0 * (1 - e) * (-1) - 1.0/2.0 * (N5n(e,n) + N8n(e,n))
	# the south-east corner
	N2n = lambda e,n: 1.0/4.0 * (1 + e) * (-1) - 1.0/2.0 * (N5n(e,n) + N6n(e,n))
	# the north-east corner
	N3n = lambda e,n: 1.0/4.0 * (1 + e) - 1.0/2.0 * (N6n(e,n) + N7n(e,n))
	# the north-west corner
	N4n = lambda e,n: 1.0/4.0 * (1 - e) - 1.0/2.0 * (N7n(e,n) + N8n(e,n))
	


	return [N1n,N2n,N3n,N4n]

def get_dN_xy(x,y):

	dNe = iso_deriv_wrt_to_e()
	dNn = iso_deriv_wrt_to_n()

	J11 = lambda e,n: dNe[0](e,n) * x[0] + dNe[1](e,n) * x[1] + dNe[2](e,n) * x[2] + dNe[3](e,n) * x[3]
	J12 = lambda e,n: dNe[0](e,n) * y[0] + dNe[1](e,n) * y[1] + dNe[2](e,n) * y[2] + dNe[3](e,n) * y[3]
	J22 = lambda e,n: dNn[0](e,n) * y[0] + dNn[1](e,n) * y[1] + dNn[2](e,n) * y[2] + dNn[3](e,n) * y[3]
	J21 = lambda e,n: dNn[0](e,n) * x[0] + dNn[1](e,n) * x[1] + dNn[2](e,n) * x[2] + dNn[3](e,n) * x[3]

	detJ = lambda e,n: J11(e,n) * J22(e,n) - J12(e,n) * J21(e,n)
	
	dNx1 = lambda e,n: ( J22(e,n) * dNe[0](e,n) - J12(e,n) * dNn[0](e,n) ) * 1/detJ(e,n)
	dNx2 = lambda e,n: ( J22(e,n) * dNe[1](e,n) - J12(e,n) * dNn[1](e,n) ) * 1/detJ(e,n)
	dNx3 = lambda e,n: ( J22(e,n) * dNe[2](e,n) - J12(e,n) * dNn[2](e,n) ) * 1/detJ(e,n)
	dNx4 = lambda e,n: ( J22(e,n) * dNe[3](e,n) - J12(e,n) * dNn[3](e,n) ) * 1/detJ(e,n)
	dNx = [dNx1,dNx2,dNx3,dNx4]
	
	dNy1 = lambda e,n: ( -J21(e,n) * dNe[0](e,n) + J11(e,n) * dNn[0](e,n) ) * 1/detJ(e,n)
	dNy2 = lambda e,n: ( -J21(e,n) * dNe[1](e,n) + J11(e,n) * dNn[1](e,n) ) * 1/detJ(e,n)
	dNy3 = lambda e,n: ( -J21(e,n) * dNe[2](e,n) + J11(e,n) * dNn[2](e,n) ) * 1/detJ(e,n)
	dNy4 = lambda e,n: ( -J21(e,n) * dNe[3](e,n) + J11(e,n) * dNn[3](e,n) ) * 1/detJ(e,n)
	dNy = [dNy1,dNy2,dNy3,dNy4]
	
	return [dNx,dNy,detJ]

def get_dN_xy_HN(x,y,N,S,E,W):
	# with hanging nodes
	dNe = iso_deriv_wrt_to_e_HN(N,S,E,W)
	dNn = iso_deriv_wrt_to_n_HN(N,S,E,W)

	J11 = lambda e,n: dNe[0](e,n) * x[0] + dNe[1](e,n) * x[1] + dNe[2](e,n) * x[2] + dNe[3](e,n) * x[3]
	J12 = lambda e,n: dNe[0](e,n) * y[0] + dNe[1](e,n) * y[1] + dNe[2](e,n) * y[2] + dNe[3](e,n) * y[3]
	J22 = lambda e,n: dNn[0](e,n) * y[0] + dNn[1](e,n) * y[1] + dNn[2](e,n) * y[2] + dNn[3](e,n) * y[3]
	J21 = lambda e,n: dNn[0](e,n) * x[0] + dNn[1](e,n) * x[1] + dNn[2](e,n) * x[2] + dNn[3](e,n) * x[3]

	detJ = lambda e,n: J11(e,n) * J22(e,n) - J12(e,n) * J21(e,n)
	
	dNx1 = lambda e,n: ( J22(e,n) * dNe[0](e,n) - J12(e,n) * dNn[0](e,n) ) * 1/detJ(e,n)
	dNx2 = lambda e,n: ( J22(e,n) * dNe[1](e,n) - J12(e,n) * dNn[1](e,n) ) * 1/detJ(e,n)
	dNx3 = lambda e,n: ( J22(e,n) * dNe[2](e,n) - J12(e,n) * dNn[2](e,n) ) * 1/detJ(e,n)
	dNx4 = lambda e,n: ( J22(e,n) * dNe[3](e,n) - J12(e,n) * dNn[3](e,n) ) * 1/detJ(e,n)
	dNx = [dNx1,dNx2,dNx3,dNx4]
	
	dNy1 = lambda e,n: ( -J21(e,n) * dNe[0](e,n) + J11(e,n) * dNn[0](e,n) ) * 1/detJ(e,n)
	dNy2 = lambda e,n: ( -J21(e,n) * dNe[1](e,n) + J11(e,n) * dNn[1](e,n) ) * 1/detJ(e,n)
	dNy3 = lambda e,n: ( -J21(e,n) * dNe[2](e,n) + J11(e,n) * dNn[2](e,n) ) * 1/detJ(e,n)
	dNy4 = lambda e,n: ( -J21(e,n) * dNe[3](e,n) + J11(e,n) * dNn[3](e,n) ) * 1/detJ(e,n)
	dNy = [dNy1,dNy2,dNy3,dNy4]
	
	return [dNx,dNy,detJ]

def xy_fct(x,y):
	# expressing x and y as functions of e and n

	N_e = iso_basis_fct()
	x_fct = lambda e,n: N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2] + N_e[3](e,n) * x[3]
	y_fct = lambda e,n: N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2] + N_e[3](e,n) * y[3]

	return [x_fct,y_fct]

def xy_fct_HN(x,y,N,S,E,W):
	# with hanging nodes
	# expressing x and y as functions of e and n
	N_e = iso_basis_fct_HN(N,S,E,W)
	x_fct = lambda e,n: N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2] + N_e[3](e,n) * x[3]
	y_fct = lambda e,n: N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2] + N_e[3](e,n) * y[3]

	return [x_fct,y_fct]

def tri_xy_fct(x,y):
	N_e = tri_iso_basis_fct()
	x_fct = lambda e,n: N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2]
	y_fct = lambda e,n: N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2]
	return [x_fct,y_fct]




def jacobian_mat(x,y):
	# defining the Jacobian matrix for the transformation from a quadrilateral element 
	# to an isoparametric element: [ dx/de dy/de; dx/dn dy/dn]

	# dx / de
	J_11 = lambda e,n: 1.0/4.0 * ( -(1 - n) * x[0] + (1 - n) * x[1] + (1 + n) * x[2] - (1 + n) * x[3] )
	# dy / de
	J_12 = lambda e,n: 1.0/4.0 * ( -(1 - n) * y[0] + (1 - n) * y[1] + (1 + n) * y[2] - (1 + n) * y[3] )
	# dx / dn
	J_21 = lambda e,n: 1.0/4.0 * ( -(1 - e) * x[0] - (1 + e) * x[1] + (1 + e) * x[2] + (1 - e) * x[3] )
	# dy / dn
	J_22 = lambda e,n: 1.0/4.0 * ( -(1 - e) * y[0] - (1 + e) * y[1] + (1 + e) * y[2] + (1 - e) * y[3] )

	J = numpy.zeros((2,2))
	J = [[J_11, J_12],[J_21, J_22]]

	return J

def jacobian_mat_HN(x,y,N,S,E,W):
	# defining the Jacobian matrix for the transformation from a quadrilateral element 
	# to an isoparametric element: [ dx/de dy/de; dx/dn dy/dn]
	dNe = iso_deriv_wrt_to_e_HN(N,S,E,W)
	dNn = iso_deriv_wrt_to_n_HN(N,S,E,W)
	# dx / de
	J_11 = lambda e,n: dNe[0](e,n) * x[0] + dNe[1](e,n) * x[1] + dNe[2](e,n) * x[2] + dNe[3](e,n) * x[3]
	# dy / de
	J_12 = lambda e,n: dNe[0](e,n) * y[0] + dNe[1](e,n) * y[1] + dNe[2](e,n) * y[2] + dNe[3](e,n) * y[3]
	# dx / dn
	J_22 = lambda e,n: dNn[0](e,n) * y[0] + dNn[1](e,n) * y[1] + dNn[2](e,n) * y[2] + dNn[3](e,n) * y[3]
	# dy / dn
	J_21 = lambda e,n: dNn[0](e,n) * x[0] + dNn[1](e,n) * x[1] + dNn[2](e,n) * x[2] + dNn[3](e,n) * x[3]
	
#	J_11 = lambda e,n: 1.0/4.0 * ( -(1 - n) * x[0] + (1 - n) * x[1] + (1 + n) * x[2] - (1 + n) * x[3] )
#	J_12 = lambda e,n: 1.0/4.0 * ( -(1 - n) * y[0] + (1 - n) * y[1] + (1 + n) * y[2] - (1 + n) * y[3] )
#	J_21 = lambda e,n: 1.0/4.0 * ( -(1 - e) * x[0] - (1 + e) * x[1] + (1 + e) * x[2] + (1 - e) * x[3] )
#	J_22 = lambda e,n: 1.0/4.0 * ( -(1 - e) * y[0] - (1 + e) * y[1] + (1 + e) * y[2] + (1 - e) * y[3] )

	J = numpy.zeros((2,2))
	J = [[J_11, J_12],[J_21, J_22]]

	return J

def tri_jacobian_mat(x,y):
#	J_11 = lambda e,n: - x[0] * (1-n) + x[1] * (1-n) - x[2] * n
#	J_12 = lambda e,n: - y[0] * (1-n) + y[1] * (1-n) - y[2] * n
#	J_21 = lambda e,n: - x[0] * (1-e) - x[1] * e + x[2] * (1-e)
#	J_22 = lambda e,n: - y[0] * (1-e) - y[1] * e + y[2] * (1-e)

	J_11 = lambda e,n: -1 * x[0] + 1 * x[1] + 0 * x[2] 
	J_12 = lambda e,n: -1 * y[0] + 1 * y[1] + 0 * y[2]
	J_21 = lambda e,n: -1 * x[0] + 0 * x[1] + 1 * x[2]
	J_22 = lambda e,n: -1 * y[0] + 0 * y[1] + 1 * y[2]

	J = numpy.zeros((2,2))
	J = [[J_11,J_12], [J_21,J_22]]
	return J





def jacobian_mat1(x,y):
  # defining the Jacobian matrix for the transformation from a quadrilateral element 
	# to an isoparametric element: [ dx/de dy/de; dx/dn dy/dn]

	# dx / de
	J_11 = lambda e,n: 1.0/4.0 * ( -(1 - n) * x[0] + (1 - n) * x[1] + (1 + n) * x[2] - (1 + n) * x[3] )
	# dy / de
	J_12 = lambda e,n: 1.0/4.0 * ( -(1 - n) * y[0] + (1 - n) * y[1] + (1 + n) * y[2] - (1 + n) * y[3] )
	# dx / dn
	J_21 = lambda e,n: 1.0/4.0 * ( -(1 - e) * x[0] - (1 + e) * x[1] + (1 + e) * x[2] + (1 - e) * x[3] )
	# dy / dn
	J_22 = lambda e,n: 1.0/4.0 * ( -(1 - e) * y[0] - (1 + e) * y[1] + (1 + e) * y[2] + (1 - e) * y[3] )
	#J = numpy.zeros((2,2))
	J = lambda e,n: [[J_11(e,n), J_12(e,n)],[J_21(e,n), J_22(e,n)]]
	return J


def determinant(J):
	# compute the determinant of the Jacobian matrix J

	#jdet = lambda e,n: J[0][0](e,n) * J[1][1](e,n) - J[0][1](e,n) * J[1][0](e,n)
	#return jdet(0,0)

	return lambda e,n: J[0][0](e,n) * J[1][1](e,n) - J[0][1](e,n) * J[1][0](e,n)

def inverse(J):
	# computes the inverse of the Jacobian matrix J

	j_det = lambda e,n: determinant(J)(e,n)
	j_inv_11 = lambda e,n: J[1][1](e,n) * 1.0/j_det(e,n)
	j_inv_12 = lambda e,n: - J[0][1](e,n) * 1.0/j_det(e,n)
	j_inv_21 = lambda e,n: - J[1][0](e,n) * 1.0/j_det(e,n)
	j_inv_22 = lambda e,n: J[0][0](e,n) * 1.0/j_det(e,n)

	J_inv = numpy.zeros((2,22))
	J_inv = [[ j_inv_11, j_inv_12],[j_inv_21, j_inv_22]]

	return J_inv

def inverse1(J):
	# computes the inverse of the Jacobian matrix J

	j_det = determinant(J)
	j_inv_11 = lambda e,n: J[1][1](e,n) * 1.0/j_det
	j_inv_12 = lambda e,n: - J[0][1](e,n) * 1.0/j_det
	j_inv_21 = lambda e,n: - J[1][0](e,n) * 1.0/j_det
	j_inv_22 = lambda e,n: J[0][0](e,n) * 1.0/j_det

	J_inv = lambda e,n: [[ j_inv_11(e,n), j_inv_12(e,n)],[j_inv_21(e,n), j_inv_22(e,n)]]
	return J_inv

def transformed_deriv(x,y):
	# defining the derivatives of the isoparametric element through the transformation relationship
	# J_inv_11 * dN_wrt_e + J_inv_12 * dN_wrt_n
	# J_inv_21 * dN_wrt_e + J_inv_22 * dN_wrt_n

	# derivatives with respect to eta
	Nie = iso_deriv_wrt_to_e()
	# derivatives with respect to niu
	Nin = iso_deriv_wrt_to_n()

	# computing the Jacobian matrix
	J = jacobian_mat(x,y)
	# inverse of the Jacobian
	J_inv = inverse(J)
	

	dN1x = lambda e,n: (J_inv[0][0](e,n) * Nie[0](e,n) + J_inv[0][1](e,n) * Nin[0](e,n) )
	dN2x = lambda e,n: (J_inv[0][0](e,n) * Nie[1](e,n) + J_inv[0][1](e,n) * Nin[1](e,n) )
	dN3x = lambda e,n: (J_inv[0][0](e,n) * Nie[2](e,n) + J_inv[0][1](e,n) * Nin[2](e,n) )
	dN4x = lambda e,n: (J_inv[0][0](e,n) * Nie[3](e,n) + J_inv[0][1](e,n) * Nin[3](e,n) )

	dN1y = lambda e,n: (J_inv[1][0](e,n) * Nie[0](e,n) + J_inv[1][1](e,n) * Nin[0](e,n) )
	dN2y = lambda e,n: (J_inv[1][0](e,n) * Nie[1](e,n) + J_inv[1][1](e,n) * Nin[1](e,n) )
	dN3y = lambda e,n: (J_inv[1][0](e,n) * Nie[2](e,n) + J_inv[1][1](e,n) * Nin[2](e,n) )
	dN4y = lambda e,n: (J_inv[1][0](e,n) * Nie[3](e,n) + J_inv[1][1](e,n) * Nin[3](e,n) )

	T_mat = [[dN1x, dN2x, dN3x, dN4x],[dN1y, dN2y, dN3y, dN4y]]
	return [T_mat,determinant(J)]

def transformed_deriv_HN(x,y,N,S,E,W):
	# with hanging nodes
	# defining the derivatives of the isoparametric element through the transformation relationship
	# J_inv_11 * dN_wrt_e + J_inv_12 * dN_wrt_n
	# J_inv_21 * dN_wrt_e + J_inv_22 * dN_wrt_n

	# derivatives with respect to eta
	Nie = iso_deriv_wrt_to_e_HN(N,S,E,W)
	# derivatives with respect to niu
	Nin = iso_deriv_wrt_to_n_HN(N,S,E,W)

	# computing the Jacobian matrix
	J = jacobian_mat_HN(x,y,N,S,E,W)
	# inverse of the Jacobian
	J_inv = inverse(J)
	

	dN1x = lambda e,n: (J_inv[0][0](e,n) * Nie[0](e,n) + J_inv[0][1](e,n) * Nin[0](e,n) )
	dN2x = lambda e,n: (J_inv[0][0](e,n) * Nie[1](e,n) + J_inv[0][1](e,n) * Nin[1](e,n) )
	dN3x = lambda e,n: (J_inv[0][0](e,n) * Nie[2](e,n) + J_inv[0][1](e,n) * Nin[2](e,n) )
	dN4x = lambda e,n: (J_inv[0][0](e,n) * Nie[3](e,n) + J_inv[0][1](e,n) * Nin[3](e,n) )

	dN1y = lambda e,n: (J_inv[1][0](e,n) * Nie[0](e,n) + J_inv[1][1](e,n) * Nin[0](e,n) )
	dN2y = lambda e,n: (J_inv[1][0](e,n) * Nie[1](e,n) + J_inv[1][1](e,n) * Nin[1](e,n) )
	dN3y = lambda e,n: (J_inv[1][0](e,n) * Nie[2](e,n) + J_inv[1][1](e,n) * Nin[2](e,n) )
	dN4y = lambda e,n: (J_inv[1][0](e,n) * Nie[3](e,n) + J_inv[1][1](e,n) * Nin[3](e,n) )

	T_mat = [[dN1x, dN2x, dN3x, dN4x],[dN1y, dN2y, dN3y, dN4y]]
	return [T_mat,determinant(J)]



def tribasisFct(C):
	N1 = lambda x,y: C[0,0] + C[1,0]*x + C[2,0]*y
	N2 = lambda x,y: C[0,1] + C[1,1]*x + C[2,1]*y
	N3 = lambda x,y: C[0,2] + C[1,2]*x + C[2,2]*y
	return [N1,N2,N3]

def basisFct(C):
	N1 = lambda x,y: C[0,0] + C[1,0]*x + C[2,0]*y + C[3,0]*x*y
	N2 = lambda x,y: C[0,1] + C[1,1]*x + C[2,1]*y + C[3,1]*x*y
	N3 = lambda x,y: C[0,2] + C[1,2]*x + C[2,2]*y + C[3,2]*x*y
	N4 = lambda x,y: C[0,3] + C[1,3]*x + C[2,3]*y + C[3,3]*x*y
	return [N1,N2,N3,N4]

def triderivX(C):
	N1x = lambda x,y: C[1,0] 
	N2x = lambda x,y: C[1,1]
	N3x = lambda x,y: C[1,2]
	return [N1x, N2x, N3x]

def derivX(C):
	N1x = lambda x,y: C[1,0] + C[3,0]*y
	N2x = lambda x,y: C[1,1] + C[3,1]*y
	N3x = lambda x,y: C[1,2] + C[3,2]*y
	N4x = lambda x,y: C[1,3] + C[3,3]*y
	return [N1x, N2x, N3x, N4x]

def triderivY(C):
	N1y = lambda x,y:  C[2,0]
	N2y = lambda x,y:  C[2,1]
	N3y = lambda x,y:  C[2,2]
	return [N1y, N2y, N3y]

def derivY(C):
	N1y = lambda x,y:  C[2,0] + C[3,0]*x 
	N2y = lambda x,y:  C[2,1] + C[3,1]*x 
	N3y = lambda x,y:  C[2,2] + C[3,2]*x 
	N4y = lambda x,y:  C[2,3] + C[3,3]*x 
	return [N1y, N2y, N3y, N4y]

def mytriderivX(C):
	print ' inside 1', C
	N1x = lambda x,y: C[1,0] 
	N2x = lambda x,y: C[1,1]
	N3x = lambda x,y: C[1,2]
	print N1x(0,324)
	return [N1x, N2x, N3x]

def mytriderivY(C):
	print 'inside 2', C
	N1y = lambda x,y:  C[2,0]
	N2y = lambda x,y:  C[2,1]
	N3y = lambda x,y:  C[2,2]
	print N1y(0,324)
	return [N1y, N2y, N3y]


# TRIANGLE QUADRATIC

def tri_iso_basis_fct_quadratic():
	# defining basis functions for an isoparametri TRIANGLE
	N1 = lambda e,n: 2 * (1 - e - n) * (1.0/2.0 - e - n)
	N2 = lambda e,n: 2 * e * (e - 1.0/2.0) 
	N3 = lambda e,n: 2 * n * (n - 1.0/2.0)
	N4 = lambda e,n: 4 * (1 - e - n) * e
	N5 = lambda e,n: 4 * e * n
	N6 = lambda e,n: 4 * n * (1 - e - n) 

	return [N1,N2,N3,N4,N5,N6]

def tri_xy_fct_quadratic(x,y):
	N_e = tri_iso_basis_fct_quadratic()
	x_fct = lambda e,n: N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2] + N_e[3](e,n) * x[3] + N_e[4](e,n) * x[4] + N_e[5](e,n) * x[5]
	y_fct = lambda e,n: N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2] + N_e[3](e,n) * y[3] + N_e[4](e,n) * y[4] + N_e[5](e,n) * y[5]
	return [x_fct,y_fct]

def tri_jacobian_mat_quadratic(x,y):
#	J_11 = lambda e,n: - x[0] * (1-n) + x[1] * (1-n) - x[2] * n
#	J_12 = lambda e,n: - y[0] * (1-n) + y[1] * (1-n) - y[2] * n
#	J_21 = lambda e,n: - x[0] * (1-e) - x[1] * e + x[2] * (1-e)
#	J_22 = lambda e,n: - y[0] * (1-e) - y[1] * e + y[2] * (1-e)

	J_11 = lambda e,n: ( 2 * (-1.5 + 2 * e + 2 * n) * x[0] + 2 * (2 * e - 0.5) * x[1] + 2* 0 * x[2] +
	 					4 * ( 1 - 2 * e - n) * x[3] + 4 * n * x[4] + 4 * (-n) * x[5])
	J_12 = lambda e,n: ( 2 * (-1.5 + 2 * e + 2 * n) * y[0] + 2 * (2 * e - 0.5) * y[1] + 2* 0 * y[2] +
	 					4 * ( 1 - 2 * e - n) * y[3] + 4 * n * y[4] + 4 * (-n) * y[5])
	J_21 = lambda e,n: ( 2 * ( -1.5 + 2 * n + 2 * e ) * x[0]+ 2 * 0 * x[1] + 2 * (2 * n - 0.5) * x[2] +
						4 * (-e) * x[3] + 4 * e * x[4] + 4 * (1 - e - 2 * n) * x[5])
	J_22 = lambda e,n: ( 2 * ( -1.5 + 2 * n + 2 * e ) * y[0]+ 2 * 0 * y[1] + 2 * (2 * n - 0.5) * y[2] +
						4 * (-e) * y[3] + 4 * e * y[4] + 4 * (1 - e - 2 * n) * y[5])

	J = numpy.zeros((2,2))
	J = [[J_11,J_12], [J_21,J_22]]
	return J


# TRIANGLE CUBIC
def tri_iso_basis_fct_cubic():
	# defining basis functions for an isoparametri TRIANGLE
	G1 = lambda e,n: 1 - e- n
	G2 = lambda e,n: e
	G3 = lambda e,n: n
	
	N1 = lambda e,n: 1.0/2.0 * G1(e,n) * (3 * G1(e,n) - 1) * (3 * G1(e,n) - 2)
	N2 = lambda e,n: 1.0/2.0 * G2(e,n) * (3 * G2(e,n) - 1) * (3 * G2(e,n) - 2)
	N3 = lambda e,n: 1.0/2.0 * G3(e,n) * (3 * G3(e,n) - 1) * (3 * G3(e,n) - 2)
	
	N4 = lambda e,n: 9.0/2.0 * G1(e,n) * G2(e,n) * ( 3 * G1(e,n) - 1)
	N5 = lambda e,n: 9.0/2.0 * G1(e,n) * G2(e,n) * ( 3 * G2(e,n) - 1)
	
	N6 = lambda e,n: 9.0/2.0 * G2(e,n) * G3(e,n) * ( 3 * G2(e,n) - 1)
	N7 = lambda e,n: 9.0/2.0 * G2(e,n) * G3(e,n) * ( 3 * G3(e,n) - 1)
	
	N8 = lambda e,n: 9.0/2.0 * G3(e,n) * G1(e,n) * ( 3 * G3(e,n) - 1)
	N9 = lambda e,n: 9.0/2.0 * G3(e,n) * G1(e,n) * ( 3 * G1(e,n) - 1)

	N10 = lambda e,n: 27 * G2(e,n) * G3(e,n) * G1(e,n)
	
	return [N1,N2,N3,N4,N5,N6, N7, N8, N9, N10]

def tri_xy_fct_cubic(x,y):
	N_e = tri_iso_basis_fct_cubic()
	x_fct = lambda e,n: ( N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2] + N_e[3](e,n) * x[3] + N_e[4](e,n) * x[4] + N_e[5](e,n) * x[5] +
						 N_e[6](e,n) * x[6] + N_e[7](e,n) * x[7] + N_e[8](e,n) * x[8] + N_e[9](e,n) * x[9])
	y_fct = lambda e,n: ( N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2] + N_e[3](e,n) * y[3] + N_e[4](e,n) * y[4] + N_e[5](e,n) * y[5] + 
						N_e[6](e,n) * y[6] + N_e[7](e,n) * y[7] + N_e[8](e,n) * y[8] + N_e[9](e,n) * y[9]  )
	return [x_fct,y_fct]

def tri_iso_deriv_wrt_to_e_cubic():
	# triangle
	# defining the derivatives with respect to eta of the isoparametric basis functions
	
	N1e = lambda e,n: -5.5 - 13.5 * e*e + e  * (18 - 27 * n) + 18 * n - 13.5 * n * n
	N2e = lambda e,n: 13.5 * e * e - 9 * e + 1
	N3e = lambda e,n: 0.0
	
	N4e = lambda e,n: (9 * (2 + 3 * n * n - 10 * e + 9 * e * e + n * (-5 + 12* e))) / 2.0
	N5e = lambda e,n: (-9 * (1 - 8 * e + 9 * e * e + n * (-1 + 6 * e))) / 2.0
	
	N6e = lambda e,n: (9 * n * (-1 + 6 * e)) / 2.0
	N7e = lambda e,n: (9 * n * (-1 + 3 * n)) / 2.0
	
	N8e = lambda e,n: (-9 * n *  (-1 + 3 * n)) / 2.0
	N9e = lambda e,n: (9 * n * (-5 + 6 * n + 6 * e)) / 2.0
	
	N10e = lambda e,n: -27 * n *  (-1 + n + 2 * e)
	
	return [N1e,N2e,N3e,N4e, N5e,N6e,N7e,N8e, N9e, N10e]

def tri_iso_deriv_wrt_to_n_cubic():
	# defining the derivatives with respect to niu of the isoparametric basis functions

	N1n = lambda e,n: -5.5 - 13.5 * n * n + n * (18 - 27 * e) + 18* e- 13.5 * e * e
	N2n = lambda e,n: 0.0
	N3n = lambda e,n: 1 - 9 * n + 13.5 * n * n
	
	N4n = lambda e,n: (9 * e * (-5 + 6 * n + 6  * e)) / 2.0
	N5n = lambda e,n: (-9 * e * (-1 + 3 * e)) / 2.0
	
	N6n = lambda e,n: (9 * e * (-1 + 3 * e)) / 2.0
	N7n = lambda e,n: (9 * e * (-1 + 6 * n)) / 2.0
	
	N8n = lambda e,n: (-9 * (1 + 9 * n * n - e + n * (-8 + 6 * e))) / 2.0
	N9n = lambda e,n: (9 * (2 + 9 * n * n - 5 * e + 3 * e * e + 2 * n * (-5 + 6 * e))) / 2.0
	
	N10n = lambda e,n: -27 * p * (-1 + 2 * n + e)


	return [N1n,N2n,N3n,N4n, N5n,N6n,N7n,N8n, N9n, N10n]


def tri_jacobian_mat_cubic(x,y):
	# defining the Jacobian matrix for the transformation from a curved edge triangle element 
	# to an isoparametric element: [ dx/de dy/de; dx/dn dy/dn]
	dNe = tri_iso_deriv_wrt_to_e_cubic()
	dNn = tri_iso_deriv_wrt_to_n_cubic()
	
	# dx / de
	J_11 = lambda e,n: ( dNe[0](e,n) * x[0] + dNe[1](e,n) * x[1] + dNe[2](e,n) * x[2] + dNe[3](e,n) * x[3] + 
						dNe[4](e,n) * x[4] + dNe[5](e,n) * x[5] + dNe[6](e,n) * x[6] + dNe[7](e,n) * x[7] +
						dNe[8](e,n) * x[8] + dNe[9](e,n) * x[9])
	# dy / de
	J_12 = lambda e,n: ( dNe[0](e,n) * y[0] + dNe[1](e,n) * y[1] + dNe[2](e,n) * y[2] + dNe[3](e,n) * y[3] + 
						dNe[4](e,n) * y[4] + dNe[5](e,n) * y[5] + dNe[6](e,n) * y[6] + dNe[7](e,n) * y[7] +
						dNe[8](e,n) * y[8] + dNe[9](e,n) * y[9] )
	# dx / dn
	J_22 = lambda e,n: ( dNn[0](e,n) * y[0] + dNn[1](e,n) * y[1] + dNn[2](e,n) * y[2] + dNn[3](e,n) * y[3] + 
						dNn[4](e,n) * y[4] + dNn[5](e,n) * y[5] + dNn[6](e,n) * y[6] + dNn[7](e,n) * y[7] +
						dNn[8](e,n) * y[8] + dNn[9](e,n) * y[9] )
	# dy / dn
	J_21 = lambda e,n: ( dNn[0](e,n) * x[0] + dNn[1](e,n) * x[1] + dNn[2](e,n) * x[2] + dNn[3](e,n) * x[3] +
						dNn[4](e,n) * x[4] + dNn[5](e,n) * x[5] + dNn[6](e,n) * x[6] + dNn[7](e,n) * x[7] +
						dNn[8](e,n) * x[8] + dNn[9](e,n) * x[9] )
	
	J = numpy.zeros((2,2))
	J = [[J_11, J_12],[J_21, J_22]]

	return J

# QUADRILATERAL Bi-Quadratic
def quad_iso_basis_fct_bi_quadratic():
	# defining basis functions for an isoparametri TRIANGLE
	N0 = lambda e,n: 1.0/4.0 * (1 - e) * (1 - n) * ( -e - n - 1 )
	N1 = lambda e,n: 1.0/2.0 * (1 - e * e) * ( 1 - n )
	N2 = lambda e,n: 1.0/4.0 * (1 + e) * (1 - n) * (e - n - 1) 
	N3 = lambda e,n: 1.0/2.0 * (1 + e) * (1 - n * n)
	N4 = lambda e,n: 1.0/4.0 * (1 + e) * (1 + n) * (e + n - 1)
	N5 = lambda e,n: 1.0/2.0 * (1 - e * e) * (1 + n)
	N6 = lambda e,n: 1.0/4.0 * (1 - e) * (1 + n) * (-e + n - 1)
	N7 = lambda e,n: 1.0/2.0 * (1 - e) * (1 - n * n)

	return [N0,N1,N2,N3,N4,N5,N6,N7]

def quad_xy_fct_bi_quadratic(x,y):
	N_e = quad_iso_basis_fct_bi_quadratic()
	x_fct = lambda e,n: ( N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2] + N_e[3](e,n) * x[3] + 
						N_e[4](e,n) * x[4] + N_e[5](e,n) * x[5] + N_e[6](e,n) * x[6] + N_e[7](e,n) * x[7] )
	y_fct = lambda e,n: (N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2] + N_e[3](e,n) * y[3] + 
						N_e[4](e,n) * y[4] + N_e[5](e,n) * y[5] + N_e[6](e,n) * y[6] + N_e[7](e,n) * y[7])
	return [x_fct,y_fct]

def quad_iso_deriv_wrt_to_e_bi_quadratic():
	# triangle
	# defining the derivatives with respect to eta of the isoparametric basis functions
	
	N0e = lambda e,n: -((-1 + n) * (n + 2 * e)) / 4.0
	N1e = lambda e,n: (-1 + n) * e	
	N2e = lambda e,n: ((-1 + n) *  (n - 2 * e)) / 4.0
	N3e = lambda e,n: (1 - n*n) / 2.0
	
	N4e = lambda e,n: ((1 + n) * (n + 2 * e)) / 4.0
	N5e = lambda e,n: -((1 + n) * e)
	
	N6e = lambda e,n: -((1 + n) * (n - 2 * e)) / 4.0
	N7e = lambda e,n: (-1 + n*n) / 2.0
	
	return [N0e,N1e,N2e,N3e,N4e, N5e,N6e,N7e]

def quad_iso_deriv_wrt_to_n_bi_quadratic():
	# defining the derivatives with respect to niu of the isoparametric basis functions

	N0n = lambda e,n: -((-1 + e) * (2 * n + e)) / 4.0
	N1n = lambda e,n: (-1 + e * e) / 2.0
	N2n = lambda e,n: ((2 * n - e) * (1 + e)) / 4.0
	N3n = lambda e,n: -(n * (1 + e))
	N4n = lambda e,n: ((1 + e) * (2 * n + e)) / 4.0
	N5n = lambda e,n: (1 - e * e) / 2.0
	N6n = lambda e,n: -((2 * n - e) * (-1 + e)) / 4.0
	N7n = lambda e,n: n * (-1 + e)

	return [N0n,N1n,N2n,N3n,N4n, N5n,N6n,N7n]

def quad_jacobian_mat_bi_quadratic(x,y):
	# defining the Jacobian matrix for the transformation from a curved edge triangle element 
	# to an isoparametric element: [ dx/de dy/de; dx/dn dy/dn]
	dNe = quad_iso_deriv_wrt_to_e_bi_quadratic()
	dNn = quad_iso_deriv_wrt_to_n_bi_quadratic()
	
	# dx / de
	J_11 = lambda e,n: ( dNe[0](e,n) * x[0] + dNe[1](e,n) * x[1] + dNe[2](e,n) * x[2] + dNe[3](e,n) * x[3] + 
						dNe[4](e,n) * x[4] + dNe[5](e,n) * x[5] + dNe[6](e,n) * x[6] + dNe[7](e,n) * x[7] )
	# dy / de
	J_12 = lambda e,n: ( dNe[0](e,n) * y[0] + dNe[1](e,n) * y[1] + dNe[2](e,n) * y[2] + dNe[3](e,n) * y[3] + 
						dNe[4](e,n) * y[4] + dNe[5](e,n) * y[5] + dNe[6](e,n) * y[6] + dNe[7](e,n) * y[7] )
	# dx / dn
	J_22 = lambda e,n: ( dNn[0](e,n) * y[0] + dNn[1](e,n) * y[1] + dNn[2](e,n) * y[2] + dNn[3](e,n) * y[3] + 
						dNn[4](e,n) * y[4] + dNn[5](e,n) * y[5] + dNn[6](e,n) * y[6] + dNn[7](e,n) * y[7] )
	# dy / dn
	J_21 = lambda e,n: ( dNn[0](e,n) * x[0] + dNn[1](e,n) * x[1] + dNn[2](e,n) * x[2] + dNn[3](e,n) * x[3] +
						dNn[4](e,n) * x[4] + dNn[5](e,n) * x[5] + dNn[6](e,n) * x[6] + dNn[7](e,n) * x[7] )
	
	J = numpy.zeros((2,2))
	J = [[J_11, J_12],[J_21, J_22]]

	return J


# QUADRILATERAL Bi-Cubic
def quad_iso_basis_fct_bi_cubic():
	# defining basis functions for an isoparametri TRIANGLE
	N1 = lambda e,n: 1.0/32.0 * (1 - e) * (1 - n) * ( -10 + 9 * ( e * e + n * n) )
	N4 = lambda e,n: 1.0/32.0 * (1 + e) * (1 - n) * ( -10 + 9 * ( e * e + n * n) )
	N7 = lambda e,n: 1.0/32.0 * (1 + e) * (1 + n) * ( -10 + 9 * ( e * e + n * n) )
	N10 = lambda e,n: 1.0/32.0 * (1 - e) * (1 + n) * ( -10 + 9 * ( e * e + n * n) )
	
	N2 = lambda e,n: 9.0/32.0 * (1 - n) * (1 - e*e) * (1 - 3 * e)
	N3 = lambda e,n: 9.0/32.0 * (1 - n) * (1 - e*e) * (1 + 3 * e)
	
	N5 = lambda e,n: 9.0/32.0 * (1 + e) * (1 - n*n) * (1 - 3 * n)
	N6 = lambda e,n: 9.0/32.0 * (1 + e) * (1 - n*n) * (1 + 3 * n)
	
	N8 = lambda e,n: 9.0/32.0 * (1 + n) * (1 - e*e) * (1 + 3 * e)
	N9 = lambda e,n: 9.0/32.0 * (1 + n) * (1 - e*e) * (1 - 3 * e)

	N11 = lambda e,n: 9.0/32.0 * (1 - e) * (1 - n*n) * (1 + 3 * n)
	N12 = lambda e,n: 9.0/32.0 * (1 - e) * (1 - n*n) * (1 - 3 * n)
	
	return [N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12]

def quad_xy_fct_bi_cubic(x,y):
	N_e = quad_iso_basis_fct_bi_cubic()
	x_fct = lambda e,n: ( N_e[0](e,n) * x[0] + N_e[1](e,n) * x[1] + N_e[2](e,n) * x[2] + N_e[3](e,n) * x[3] + 
						N_e[4](e,n) * x[4] + N_e[5](e,n) * x[5] + N_e[6](e,n) * x[6] + N_e[7](e,n) * x[7] +
						N_e[8](e,n) * x[8] + N_e[9](e,n) * x[9] + N_e[10](e,n) * x[10] + N_e[11](e,n) * x[11])
	y_fct = lambda e,n: (
						N_e[0](e,n) * y[0] + N_e[1](e,n) * y[1] + N_e[2](e,n) * y[2] + N_e[3](e,n) * y[3] + 
						N_e[4](e,n) * y[4] + N_e[5](e,n) * y[5] + N_e[6](e,n) * y[6] + N_e[7](e,n) * y[7] +
						N_e[8](e,n) * y[8] + N_e[9](e,n) * y[9] + N_e[10](e,n) * y[10] + N_e[11](e,n) * y[11]
						)
	return [x_fct,y_fct]

def quad_iso_deriv_wrt_to_e_bi_cubic():
	# triangle
	# defining the derivatives with respect to eta of the isoparametric basis functions
	
	N0e = lambda e,n: 0.28125 * (-1 + n) * (-1.11111 + 1 * n * n - 2 * e + 3 * e * e)
	N3e = lambda e,n: -0.28125 * (-1 + n) * (-1.11111 + 1 * n * n + 2 * e + 3 * e * e)
	N6e = lambda e,n: 0.28125 * (1 + n) * (-1.11111 + 1 * n * n+ 2 * e + 3 *  e * e)
	N9e = lambda e,n: -0.28125 * (1 + n) * (-1.11111 + 1 * n * n - 2 * e + 3 * e * e)
	
	N1e = lambda e,n: -0.84375 - 0.5625 * e + 2.53125 * e * e + n * (0.84375 + 0.5625 * e - 2.53125 * e * e)
	N2e = lambda e,n: 0.84375 - 0.5625 * e - 2.53125 * e * e + n * (-0.84375 + 0.5625 * e + 2.53125 * e * e)
	
	N4e = lambda e,n: 0.84375 * (-0.333333 + 1 * n) * (-1 + n * n)
	N5e = lambda e,n: -0.84375 * (0.333333 + 1 * n) * (-1 + n * n)
	
	N7e = lambda e,n: -2.53125 * (1 + n) * (-0.333333 + 0.222222 * e + 1  * e * e)
	N8e = lambda e,n: 2.53125 * (1 + n) * (-0.333333 - 0.222222 * e + 1 * e * e)

	N10e = lambda e,n: 0.84375 * (0.333333 + 1 * n) * (-1 + n * n)
	N11e = lambda e,n: -0.84375 * (-0.333333 + 1 * n) * (-1 + n * n)
	
	return [N0e,N1e,N2e,N3e,N4e, N5e,N6e,N7e, N8e, N9e, N10e, N11e]

def quad_iso_deriv_wrt_to_n_bi_cubic():
	# defining the derivatives with respect to niu of the isoparametric basis functions

	N0n = lambda e,n: 0.28125 * (-1 + e) * (-1.11111 - 2 * n + 3 * n * n + 1 * e * e)
	N3n = lambda e,n: -0.28125 * (1 + e) * (-1.11111 - 2 * n + 3* n * n + 1 * e * e)
	N6n = lambda e,n: 0.28125 * (1 + e) * (-1.11111 + 2 * n + 3 * n * n + 1 * e * e)
	N9n = lambda e,n: -0.28125 * (-1 + e) * (-1.11111 + 2 * n + 3 * n * n + 1 * e * e)
	
	N1n = lambda e,n: -0.84375 * (-0.333333 + 1 * e) * (-1 + e * e)
	N2n = lambda e,n: 0.84375 * (0.333333 + 1 * e) * (-1 + e * e)
	
	N4n = lambda e,n: 2.53125 * (-0.333333 - 0.222222 * n + 1 * n * n) * (1 + e)
	N5n = lambda e,n: -2.53125 * (-0.333333 + 0.222222 * n + 1 * n * n) * (1 + e)
	
	N7n = lambda e,n: -0.84375 * (0.333333 + 1 * e) * (-1 + e * e)
	N8n = lambda e,n: 0.84375 * (-0.333333 + 1 * e) * (-1 + e * e)

	N10n = lambda e,n: 0.84375 + n * (-0.5625 + 0.5625 * e) - 0.84375 * e + n * n * (-2.53125 + 2.53125 * e)
	N11n = lambda e,n: -0.84375 + n * n *(2.53125 - 2.53125 * e) + n * (-0.5625 + 0.5625 * e) + 0.84375 * e

	return [N0n,N1n,N2n,N3n,N4n, N5n,N6n,N7n, N8n, N9n, N10n, N11n]

def quad_jacobian_mat_bi_cubic(x,y):
	# defining the Jacobian matrix for the transformation from a curved edge triangle element 
	# to an isoparametric element: [ dx/de dy/de; dx/dn dy/dn]
	dNe = quad_iso_deriv_wrt_to_e_bi_cubic()
	dNn = quad_iso_deriv_wrt_to_n_bi_cubic()
	
	# dx / de
	J_11 = lambda e,n: ( dNe[0](e,n) * x[0] + dNe[1](e,n) * x[1] + dNe[2](e,n) * x[2] + dNe[3](e,n) * x[3] + 
						dNe[4](e,n) * x[4] + dNe[5](e,n) * x[5] + dNe[6](e,n) * x[6] + dNe[7](e,n) * x[7] +
						dNe[8](e,n) * x[8] + dNe[9](e,n) * x[9] + dNe[10](e,n) * x[10] + dNe[11](e,n) * x[11])
	# dy / de
	J_12 = lambda e,n: ( dNe[0](e,n) * y[0] + dNe[1](e,n) * y[1] + dNe[2](e,n) * y[2] + dNe[3](e,n) * y[3] + 
						dNe[4](e,n) * y[4] + dNe[5](e,n) * y[5] + dNe[6](e,n) * y[6] + dNe[7](e,n) * y[7] +
						dNe[8](e,n) * y[8] + dNe[9](e,n) * y[9] + dNe[10](e,n) * y[10] + dNe[11](e,n) * y[11])
	# dx / dn
	J_22 = lambda e,n: ( dNn[0](e,n) * y[0] + dNn[1](e,n) * y[1] + dNn[2](e,n) * y[2] + dNn[3](e,n) * y[3] + 
						dNn[4](e,n) * y[4] + dNn[5](e,n) * y[5] + dNn[6](e,n) * y[6] + dNn[7](e,n) * y[7] +
						dNn[8](e,n) * y[0] + dNn[9](e,n) * y[9] + dNn[10](e,n) * y[10] + dNn[11](e,n) * y[11])
	# dy / dn
	J_21 = lambda e,n: ( dNn[0](e,n) * x[0] + dNn[1](e,n) * x[1] + dNn[2](e,n) * x[2] + dNn[3](e,n) * x[3] +
						dNn[4](e,n) * x[4] + dNn[5](e,n) * x[5] + dNn[6](e,n) * x[6] + dNn[7](e,n) * x[7] +
						dNn[8](e,n) * x[8] + dNn[9](e,n) * x[9] + dNn[10](e,n) * x[10] + dNn[11](e,n) * x[11] )
	
	J = numpy.zeros((2,2))
	J = [[J_11, J_12],[J_21, J_22]]

	return J