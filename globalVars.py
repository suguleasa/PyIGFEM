#binBnd = [0,180,256] # for real slice
binBnd = [0,110,256] # for fake images
#binBnd = [0,50,256] # for doughnut image
AREA_INCLUSION = 4.0 # or a 4x4 small inclusion of pixels
PROB = 0.3

DIV_F = 1000.0

#real slice
#MAX_SIZE = 85
#MIN_SIZE = 20#6

# fake data
MAX_SIZE_X = 33
MAX_SIZE_Y = 33
MIN_SIZE = 33
#MAX_SIZE_X = 17
#MAX_SIZE_Y = 17
#MIN_SIZE =  17 #17, 33, 64, 127, 256

ALT_MIN_SIZE = 8
TOL_LINEARS = 2.0
TOL_QUADRATICS = 2.0 # cancel the cubic approximation by making a negative number
TOL_CUBICS = 2.0 # cancel the cubic approximation by making a negative number

TOL_NURBS = 3.0

STRESS_MIN = 2 # minimum number of elements between interfaces
TOL_error = 3.0 # tolerance error for pixel distances

EPS_FACTOR = 1e-2


# in elements where no linear, quadratic or cubic polynomial could approximate
# choose the polynomial 
# 0 for linears
# 1 for quadratics
# 2 for cubics
# 3 for NURBS 
POL_APPROX = 2

# activate norm computation: 1, else 0
NORM_COMP = 1

# fake circles radii
Rs = 1.0/3.0
R1 = 1.0/6.0
R2 = 1.0/6.0
