#binBnd = [0,180,256] # for real slice
binBnd = [0,110,256] # for fake images
AREA_INCLUSION = 16.0 # or a 4x4 small inclusion of pixels
PROB = 0.05

#real slice
#MAX_SIZE = 85
#MIN_SIZE = 20#6

# fake data
MAX_SIZE = 512     
MIN_SIZE = 128

ALT_MIN_SIZE = 64
TOL_LINEARS = -1.0
TOL_QUADRATICS = 120.0 # cancel the cubic approximation by making a negative number
TOL_CUBICS = 120.0 # cancel the cubic approximation by making a negative number
STRESS_MIN = 4 # minimum number of elements between interfaces
TOL_error = 3