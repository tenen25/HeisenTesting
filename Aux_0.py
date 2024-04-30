import math
import numpy as np
from itertools import chain
import time

# Calculates product (a,b,c)(x, y, z) in the Heisenberg Group
def MH(x,y):
    return np.array([x[0]+y[0],x[1]+y[1],x[2]+y[2]+1/2*x[0]*y[1]+2*x[1]*y[0]])

# Calculates inverse of (a,b,c) in the Heisenberg Group
def Inv(x):
    return np.array([-x[0],-x[1],-x[2]])

# Calculates Koranye norm of (a,b,c) in the Heisenberg Group
def Norm(x):
    return ((x[0]**2 + x[1]**2)**2 + x[2]**2)**(1/4)

def Shift(A,B):
    A = np.array(A)
    B = np.array(B)
    # We add the jth element of A to each element of B, then we crate a new axis and do the same with the (j+1)th
    # element. This is done by the parameter "None". The third component just says that we take the whole array in the
    # last dimension. Since we defined before A = np.array(A), this is the vector of length 3
    AB = A[:, None, :] + B[None, :, :]
    # The code c -= a means c = c - a. So we substract in the following 
    # something from the second component
    AB[:, :, 2] -= 2 * (A[:, None, 0] * B[None, :, 1] - A[:, None, 1] * B[None, :, 0])
    return list(chain.from_iterable(AB))



def AltDilatV(lamb):                    # Compute V(1,e) intersected with
                                        # the lattice H_3(2\mathbb{Z}) for stretch factor
                                        # lamb
    D = np.array([lamb, lamb, lamb**2])
    ran_arr = np.array([[0,0,0],[0,0,0]])
    for i in range(3):
        # cases for even or odd underlying stretch factor
        if lamb % 2 == 0:
            ran_arr[0][i] = -D[i]
            ran_arr[1][i] = D[i]-2
        else:
            ran_arr[0][i] = -D[i]+1
            ran_arr[1][i] = D[i]-1
    DV1 = set()
    # Run over previous ranges of values and add lattice points to the set
    for x in range(ran_arr[0][0], ran_arr[1][0]+1,2):
        for y in range(ran_arr[0][1], ran_arr[1][1]+1,2):
            for z in range(ran_arr[0][2], ran_arr[1][2]+1,2):
                DV1.add((x, y, z))
    DV1 = sorted(DV1)
    return DV1


  
# return V(N,init_set) intersected with the lattice for stretch factor 'lamb' and 'init_set' as standard set
def FundemDomIter(lamb, N, init_set):
    # D is the dilation vector
    D = np.array([lamb, lamb, lamb ** 2])
    # W1 is the set V(1) intersected with the lattice as a standard set 
    W1= AltDilatV(lamb)
    W_N = []
    M = init_set
    # For each n from 1 to N we compute V(n,e) from V(n-1,e)
    for n in range(1, N + 1):
        # We dilate all lattice points in V by the dilation vector D which is given by V*D. 
        #Then we add the lattice points
        # of one dilated cell to each point in V*D
        W_N = Shift(M * D, W1)
        # Now we have computed from the previous support the next support of 
        #the substitution. So we set V = V_N
        M = W_N
    return W_N

def N_Dilat(lamb, N):       # return the lattice points in the N dilations of 
                            # the fundamental domain
    # D is the dilation vector
    D = [lamb,lamb,lamb**2]
    ran_arr = np.array([[0, 0, 0], [0, 0, 0]])
    for i in range(3):
        # cases for even or odd underlying stretch factor
        if lamb % 2 == 0:
            ran_arr[0][i] = -(D[i] ** N)
            ran_arr[1][i] = D[i] ** N - 2
        else:
            ran_arr[0][i] = -(D[i] ** N) + 1
            ran_arr[1][i] = D[i] ** N - 1
    DV_N = []
    # Run over previous ranges of values and add lattice points to the set
    for x in range(ran_arr[0][0], ran_arr[1][0] + 1, 2):
        for y in range(ran_arr[0][1], ran_arr[1][1] + 1, 2):
            for z in range(ran_arr[0][2], ran_arr[1][2] + 1, 2):
                DV_N.append(np.array([x, y, z]))
    DV_N = set(map(tuple, DV_N))
    return DV_N


def XY_Lattice_Approximant(x, lamb, N):         # returns set of all closest 
                                                # xy-coordinates of valid 
                                                # points in the lattice D^N of
                                                # the lattice Gamma
    gamma_xy = np.array( [x[0], x[1]] )
    # D_N is a factor we divide against to see whether a point is in D^N[V]
    D_N = np.array([2 * lamb ** N, 2 * lamb ** N])
    gamma_init = np.divide(np.asarray(gamma_xy), D_N)
    ranges = [ [] for i in range(2)]
    for i in range(2):
        if gamma_init[i] % 1 == 0:
            ranges[i] = range( int(gamma_init[i]) - 1, int(gamma_init[i]) + 2)
        else:
            ranges[i] = range(math.floor(gamma_init[i]), math.ceil(gamma_init[i]) + 1)
    approx_set = set()
    # we run over nearest XY pairs permissible in D^N[V] and save them to approx_set
    for x in ranges[0]:
        for y in ranges[1]:
            temp_gamma = np.array([x, y])
            temp_gamma = temp_gamma * D_N
            approx_set.add(tuple(temp_gamma))
    return approx_set


def ShiftFaces(F, vect):                    # returns a shift of the 'z'-faces
                                            # of a box 'F' by a shift of 'vect'
    F_shift = [0]*len(F)
    # Generate list of list to work on upper and lower 'z'-values separately
    F_temp = [ [], [] ]
    for  direc in range(2):
        # Shift each face separately
        F_temp[direc] = Shift([vect], F[direc:len(F):2])
        # Insert two faces to one list at odd or even indices
        F_shift[direc:len(F):2] = F_temp[direc]
    return F_shift



def DilatEdge(lamb):                        #return list of tuples corresponding 
                                            # to the 'z'-faces of V(1) as function
                                            # of stretch factor 'lamb'
    D = np.array([lamb, lamb**2])
    ran_arr = np.array([[0,0],[0,0]])       #ran_arr will correspond to 
                                            # range of coordinates in V(1)
                                            # ran_arr[0] is XY range value
                                            # ran_arr[1] is 'z' range value
    for i in range(2):
        # cases for even or odd underlying stretch factor
        if lamb % 2 == 0:
            ran_arr[i][0] = -D[i]
            ran_arr[i][1] = D[i]-2
        else:
            ran_arr[i][0] = -D[i]+1
            ran_arr[i][1] = D[i]-1
    DB = []
    for x in range(ran_arr[0][0], ran_arr[0][1]+1,2):
        for y in range(ran_arr[0][0], ran_arr[0][1]+1,2):
            for direc in range(2):          
                DB.append ( (x,y,ran_arr[1][direc]) )   
    return DB


def FundemFaceIter(lamb, N, M):     # return list of numpy arrays corresponding 
                                    # to the 'z'-faces of V(N,M)
    # D is the dilation vector
    D = np.array([lamb, lamb, lamb ** 2])
    # W is the set V(1) intersected with the lattice saved as 'z'-intervals
    W = DilatEdge(lamb)
    F_N = []
    # Generate list of lists to work on upper and lower 'z'-values separately
    M_temp = [ [], [] ]
    W_temp = [ [], [] ]
    F_temp = [ [], [] ]
    # For each n from 1 to N we compute V(n,e) from V(n-1,e), all as 'z'-intervals
    for n in range(1, N + 1):
        leng = int( (len(M)*len(W))/2 )
        F= [0] * leng               # a blank list of length 'leng'
        for direc in range(2):      # direc corresponds to odd or even index
                                    # indicating upper or lower 'z'-value
            M_temp[direc] = M[direc:len(M):2] 
            W_temp[direc] = W[direc:len(W):2] 
            # We shift the dilated face of 'M' by the face of W
            F_temp[direc] = Shift(M_temp[direc] * D, W_temp[direc])
            # We put both upper and lower faces in one list
            F[direc:leng:2] = F_temp[direc]
        if n == N:
            F_N = F
        else:
            M = F           # Move on to compute the next step
    return F_N

def GenerateFaces(HeightDown, HeightUp,xystart, xyend):
                                    # returns a box with XY coordinates ranging
                                    # from 'xystart' to 'xyend' with extreme
                                    # faces at 'HeightUp' and 'HeightDown'
    Height = [HeightDown, HeightUp]
    Faces = []
    XY_Slice = [  val for val in range(xystart, xyend+2, 2) ]
    for x in XY_Slice:
        for y in XY_Slice:
            for direc in range(2):
                Faces.append ( (x, y, Height[direc]) )
    return Faces


############################# Other proposed functions #########################


def DecompLatt(x, lamb, N):                     # returns the decomposition
                                                # of 'x' as a product  of 
                                                # elements in the dilated 
                                                # lattice and fundemantal domain
    gamma = (0,0,0)
    v = (0,0,0)
    for i in range(2):
        up = np.ceil(x[i]/lamb)
        down = np.floor(x[i]/lamb)
        if x[i]/lamb-down>1:
            gamma[i] = down
        else: gamma[i] = up
        v[i] = x[i]-lamb*gamma[i]
    temp = x[2]/(lamb**2)-1/(2*lamb)*(gamma[0]*v[1]-gamma[1]*v[0])
    up = np.ceil(temp)
    down = np.floor(temp)
    if temp-down > 1:
        gamma[2]= down
    else: gamma[2]=up
    v[2] = int(temp*lamb**2 -gamma[2]*lamb**2)
    return[ gamma, v ]




def ListOfTup(lst):             # Convert list of arrays to list of tuples
    return list(map(tuple, lst))

def DiffReturn(lst):            #Return pairs with different z-coordinate
    aux_lst = []
    for elem in lst:
        if elem[0][2]!=elem[1][2]:
            aux_lst.append(elem)
    return aux_lst

def LatticeCheck(gamma, lamb, N):       #checks whether gamma is in the N dilation of the lattice
    check = True
    D_N = np.array([lamb ** N, lamb ** N, lamb ** (2 * N)])
    eta = np.divide(gamma, 2*D_N)
    for x in eta:
        if x%1!= 0: check = False
    return check


def Project_to_Axis(var_set, axis_num):    #Project a set onto an axis to see values
    projected_set = set()
    for element in var_set:
        projected_set.add(element[axis_num])
    projected_list = sorted(projected_set)
    return projected_list


# A program to save pairs of gamma and x values satisfying the desired set inclusions
def Save_Arrays_to_Text(arr, N, known_name, check_name):    

    file = open(check_name+"_to_"+known_name, "w+")
    start_str = "The condition is satisfied for " + check_name + " with respect to " \
    +  known_name +"and " + str(N) + "with the following shift information:\n"
    file.write(start_str)
    init_str1 = "For x="
    init_str2 = ",  the corresponding gamma is "
    for element in arr:
        curr_str = init_str1 + str(element[0]) + init_str2 + str(element[1]) + ".\n"
        file.write(curr_str)
    file.close()
    return
