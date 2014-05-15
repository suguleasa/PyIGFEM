#binBnd = [0,180,256] # for real slice
binBnd = [0,110,256] # for fake images
AREA_INCLUSION = 16.0 # or a 4x4 small inclusion of pixels
PROB = 0.05

#real slice
#MAX_SIZE = 85
#MIN_SIZE = 20#6

# fake data
MAX_SIZE_X = 256
MAX_SIZE_Y = 256    
MIN_SIZE = 256

ALT_MIN_SIZE = 16
TOL_LINEARS = 2.0
TOL_QUADRATICS = 2.0 # cancel the cubic approximation by making a negative number
TOL_CUBICS = 2.0 # cancel the cubic approximation by making a negative number

TOL_NURBS = 3.0

STRESS_MIN = 4 # minimum number of elements between interfaces
TOL_error = 3.0 # tolerance error for pixel distances

EPS_FACTOR = 1e-2


# in elements where no linear, quadratic or cubic polynomial could approximate
# choose the polynomial 
# 0 for linears
# 1 for quadratics
# 2 for cubics
# 3 for NURBS - not yet implemented
POL_APPROX = 0