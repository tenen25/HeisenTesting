
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


# Set the underlying stretch factor
lamb = 3



            

#################### Auxiliary computations #####################

# We define the dilation vector
D = np.array([lamb, lamb, lamb**2])
# V equals V(0,e) intersected with the lattice
V = [np.array([0, 0, 0])]
# T will be the suspected testing domain


# Coordinates for the XY projection of the boxes
x_start = -2
x_end = 0

# The following two sets are sets which we hope are smaller
# testing domains
T_fin = GenerateFaces(-2, 0,x_start,x_end)


#The following are sets which we can show are testing domains
T_2 = GenerateFaces(-2*2, 2,x_start,x_end)

T_1 = GenerateFaces(-3*2, 3*2,x_start,x_end)

# K is known theoretic testing domain for 'lamb' greater than 4
K = DilatEdge(lamb)

# K3 is known theoretic testing domain for 'lamb' 3
K3 = FundemFaceIter(3, 2, V)

#########################################################


#################### User input #########################

# Set the iteration number for the algorithm
N = 1

#Plug in the suspected testing domain
init_set = T_3                   # set which we see whether is testing domain


#Plug in the known testing domain
known_set = T_2               # set which is 'known' to be a testing domain



#Check the inclusion condition
gamma_lst = InclCheck(lamb, N, init_set, known_set)

###################### List of gammas for x ########################
   
gamma_lst = gamma_lst[1]             #the list of tuples of (x,gamma)
                                     #satisfying the set inclusion condition
