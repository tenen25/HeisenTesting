
#Needed libraries:
import math
import numpy as np
import time
from Aux_0 import *
from Aux_Check import *



# Needed paramaters from the user are: 
# the stretch factor 'lamb',
# the iteration number 'N',
# the known testing domain 'known_set',
# the suspected testing domain 'init_set'


###################### User Input ###########################

# Set the underlying stretch factor
lamb = 4

# Set the iteration number for the algorithm
N = 2

            

#################### Auxiliary computations #####################

# V equals V(0,e) intersected with the lattice
V = [np.array([0, 0, 0])]



# Coordinates for the XY projection of the boxes
x_start = -2
x_end = 0

# The following sets are sets which we hope are smaller
# testing domains
T_fin = GenerateFaces(-2, 0,x_start,x_end)


#The following are sets which we can show are testing domains
T_1 = GenerateFaces(-6*2, 6*2,-2,0)

T_2 = GenerateFaces(-3*2, 3*2,-2,0)

T_3 = GenerateFaces(-3*2, 2*2,-2,0)
# GenerateFaces is a function improted from Aux_0 to generate box-like sets

# K is a proven testing domain for 'lamb' greater than 4
K = DilatEdge(lamb)

# K3 is a proven testing domain for 'lamb' 3
K3 = FundemFaceIter(3, 2, V)


#################### Suggested Run Options  #########################

### We include several suggestsions to run with the function 'InclCheck'
### Uncomment only ONE of the following options to run the function with the desired parameters 


# Option 1 #
## result = InclCheck(lamb, 1, T_1, K)

# Option 2 #
## result = InclCheck(lamb, 1, T_2, T_1)

# Option 3 #
## result = InclCheck(3, 2, K, K3)

# Option 4 #
## result = InclCheck(3, 2, T_2, K)

# Option 5 #
## result = InclCheck(3, 2, T_3, T_2)

# Option 6 #
## result = InclCheck(lamb, N, T_fin, T_2)


###################### List of gammas for x ########################
   
gamma_lst = result[1]             
#the list of tuples of (x,gamma) satisfying the set inclusion condition
