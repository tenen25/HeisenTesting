
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

# We define the dilation vector
D = np.array([lamb, lamb, lamb**2])
# V equals V(0,e) intersected with the lattice
V = [np.array([0, 0, 0])]
# T will be the suspected testing domain


# Coordinates for the XY projection of the boxes
x_start = -2
x_end = 0

# The following sets are sets which we hope are smaller
# testing domains
T_fin = GenerateFaces(-2, 0,x_start,x_end)


#The following are sets which we can show are testing domains
T_2 = GenerateFaces(-2*2, 2,x_start,x_end)

T_1 = GenerateFaces(-3*2, 3*2,x_start,x_end)

# K is a proven testing domain for 'lamb' greater than 4
K = DilatEdge(lamb)

# K3 is a proven testing domain for 'lamb' 3
K3 = FundemFaceIter(3, 2, V)


#################### Suggested Run Options  #########################

### We include several suggestsions to run with the function 'InclCheck'
### Uncomment only ONE of the following options to run the function with the desired parameters 


# Option 1 #
## result = InclCheck(lamb, 2, T_1, K)

# Option 2 #
## result = InclCheck(3, 2, K, K3)

# Option 3 #
## result = InclCheck(3, 2, T_1, K)

# Option 4 #
## result = InclCheck(lamb, N, T_2, T_1)

# Option 5 #
## result = InclCheck(lamb, N, T_fin, T_2)


###################### List of gammas for x ########################
   
gamma_lst = result[1]             
#the list of tuples of (x,gamma) satisfying the set inclusion condition
