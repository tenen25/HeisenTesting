import math
import numpy as np
import time
from Aux_0 import *


def Check_x(x, lamb, N, shifted, fixed):         
   # for a fixed 'x' check whether shifting 'shifted' by 'x' can be contained by a shift of
   # of 'fixed' by some 'gamma', in the N-dilated lattice
   # returns list with truth value for the condition and suitable 'gamma' if found  
    start_time = time.time()
    check = False
    H = np.array([2 * lamb ** N, 2 * lamb ** N , 2 * lamb ** (2*N)])
    # H corresponds to possible jumps between adjacent element on the N-dilated lattice
    # //The next line implements line 11 in Algorithm 1//
    rad_search =  (lamb**2) * (fixed[1][2]-fixed[0][2]+2)
    # 'rad_search' is equal to '4* lamb**(2*N)', which is how much we should shift the 'z' intervals of
    # of 'fixed'
    print("The radius for searching gamma in the z-direction is " + str(rad_search))
    # //The following lines implement lines 14 to 16 in Algorithm 1// 
    # Generate the suspected XY values for 'gamma' which are closest to 'x'
    search_set = XY_Lattice_Approximant(x, lamb, N)
    # 'z_center' is the closest 'z'-values in the N-dilated lattice below x[2]
    z_center = int(x[2]/ H[2])
    # we search possible 'gamma[2]' values concenterically around 'z_center'
    # in jumps corresponding to H[2]
    for r in range(0, int(rad_search), H[2] ):
        for sign in range(2):
           # The next line is to avoid redundancy when r==0
           if r == 0 and sign==1: continue
            z_temp = (-1)**sign * r + z_center 
            for xy_gamma in search_set:
                    # Generate 'gamma' with gamma[2] = z_temp ranging on closest
                    # possible XY values
                    gamma = (int(xy_gamma[0]), int( xy_gamma[1]), int(z_temp))
                    # \\The following line implements line 17 in Algorithm 1\\
                    Faces_shift = ShiftFaces(fixed, gamma)   # Shifts 'fixed' by the generated 'gamma'
                    # \\The following line implements line 18 in Algorithm 1\\
                    check = CheckContain(Faces_shift, shifted)   # check the containment of shifted in Faces_shift
                    
                    if check:
                        gamma_return = gamma    #The correct gamma for said x
                        end_time = time.time()
                        print("Run time to check at x=" + str(x) + " is " + str(end_time - start_time))
                        # \\The following implements line 19 in Algorithm 1\\
                        return [True, gamma_return]     # return value if a
                                                        # gamma is found
   # In case no suitable 'gamma' is found \\Implements line 22 in Algorithm 1\\
   if check == False: print("For " + str(x) +" with N="+str(N)+ " and lamb="+ str(lamb)+
                             ", there is no corresponding gamma.")
                                                # In case no gamma is found
                                                # prints message
    end_time = time.time()
    print("Run time to check at x="+ str(x)+ " is "  + str(end_time - start_time))
    return [check, None]       # return value if no gamma is found


def InclCheck(lamb, N, init, known):         
  # checks whether the set inclusion hold for all x in the set D^N[V] intersected with Gamma
    start_time = time.time()
    D = np.array([lamb, lamb, lamb**2])         #The dilation in each coordinate
    # //The next line implements line 9 in Alogrithm 1//
    DV_N = N_Dilat(lamb, N)                     #DV_N is is D^N[V] intersected with Gamma
    # //The next line implements line 8 in Algorithm 1//
    W = FundemFaceIter(lamb, N, init)           #The 'z'-faces of the V(N,init) set
    end_time = time.time()
    print("Finished preparatory computations after "+ str(end_time-start_time))
    start_time = time.time()

    #This is the list of gamma corresponding to 'x'-s satisfying the necessary inclusion
    gammaTox_lst = [ ]
    # //The next loop implements line 12 in algorithm 1//
    for x in DV_N:
    # We loop on all x in D^N[V] to check the condition
        tru_val = False

        # //The next line implements line 13 in Algorithm 1//
        K_shift = ShiftFaces(known, x)     #  shift the list 'known' by a lattice vector 'x' 
        #  checks whether there is a 'gamma' satisfying the desired set inclusion 
        # for the current 'x'
        x_output = Check_x(x, lamb, N, K_shift, W)
        tru_val = x_output[0]
        if tru_val == False: break      # Break loop in case there exists 'x' has no corresponding 'gamma'
        # Append x and corresponding gamma satisfying inclusion to the list 'gammaTox_lst'   

        gammaTox_lst.append( [x, x_output[1] ] )
    print("The statement is " + str(tru_val) +" for N="+str(N)+" and lamb="+str(lamb)+".")
    
    end_time = time.time()
    print("Run time for the program in this construction is " + str(end_time-start_time)+".")
    return [tru_val, gammaTox_lst]



def XY_Range(lst):                      # for a list of triples corresponding 
                                        # 'z'-intervals returns all 'xy' values
    arr = np.asarray(lst)
    temp_lst = np.delete(arr, 2,1)
    xy_lst = temp_lst[0::2]
    return list(xy_lst)
    
            
def CheckContain(lst_set1, lst_set2):    
   #determines whether 'z'-intervals defined by lst_set1 by is superset
   # \\ This program implements line 18 in Algorithm 1, with respect to how we saved the sets\\
   # of lst_set2  
    time_start_tot = time.time()
   # saves the possible XY values of elements in lst_set2
    xy_lst = XY_Range(lst_set2)
   # Generate list of indices for XY values in lst_set1
    ind_lst = [ FindXY(lst_set1, xy_tup) for xy_tup in xy_lst ]
    for ind in range(0, len(lst_set2), 2):
        # saves the XY values of 'lst_set2[ind]' as 'xy_tup'

        xy_tup =  ( lst_set2[ind][0], lst_set2[ind][1]  ) 
        j = FindXY(lst_set1, xy_tup)
        # if there is no interval in lst_set1 with xy_tup entries
        if j == None: 
            tru_val=False
            break
        # check whether the intervals in lst_set1 contain the intervals in lst_set2
        tru_val = lst_set1[j][2] <= lst_set2[ind][2] and lst_set1[j+1][2] >= lst_set2[ind+1][2]
        if tru_val == False: break
    time_end_tot = time.time()
    print("Run time for the check function is " + str(time_end_tot - time_start_tot))
    #tru_val is TRUE if all interval inlusions hold, FALSE if one inclusion fails
    return tru_val

def FindXY(lst, tup):                       # find first index of element with
                                            # 'xy' values equal to tup
    for ind in range(len(lst) ):
        if lst[ind][0]== tup[0] and lst[ind][1]==tup[1]:
            return ind
    return None

##################################### Other proposed functions ################

def FailReturn(lst):                        # return list of elements for which
                                            # no gamma values were found
    fail_lst = []
    for elem in lst:
        if elem[1] == None:
            fail_lst.append(elem[0])
    return fail_lst

def FailIncCheck(lamb, N, init, known):         # given the condition does not
                                                # hold for 'N', returns list
                                                # of gammas for which the set 
                                                #inclusions hold
    start_time = time.time()
    D = np.array([lamb, lamb, lamb**2]) 
    DV_N = N_Dilat(lamb, N)
    W = FundemFaceIter(lamb, N, init)
    end_time = time.time()
    print("Finished preparatory computations after "+ str(end_time-start_time))
    start_time = time.time()
    gammaTox_lst = [ ]
    
    for x in DV_N:
        tru_val = False
        K_shift = ShiftFaces(known, x)
        x_output = Check_x(x, lamb, N, K_shift, W)
        tru_val = x_output[0]
        gammaTox_lst.append( [x, x_output[1] ] )
    end_time = time.time()
    print("Run time for the program in this construction is " + str(end_time-start_time)+".")
    return [tru_val, gammaTox_lst]
