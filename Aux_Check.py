import math
import numpy as np
import time
from Aux_0 import *


def NewCheck_x(x, lamb, N, shifted, fixed):         # for a fixed 'x' checks
                                                    # whether the required set 
                                                    # inclusion holds
    start_time = time.time()
    check = False
    H = np.array([2 * lamb ** N, 2 * lamb ** N , 2 * lamb ** (2*N)])
                    # gaps between on the dilated lattice
    rad_search =  (lamb**2) * (fixed[1][2]-fixed[0][2]+2)
    #rad_search = 4* lamb**(2*N)
    print("The radius for searching gamma in the z-direction is " + str(rad_search))
    search_set = XY_Lattice_Approximant(x, lamb, N)
    z_center = int(x[2]/ H[2])
    for r in range(0, int(rad_search), H[2] ):
        for sign in range(2):
            z_temp = (-1)**sign * r + z_center 
            #else:
            for xy_gamma in search_set:
                    gamma = (int(xy_gamma[0]), int( xy_gamma[1]), int(z_temp))
                    Faces_shift = ShiftFaces(fixed, gamma)
                    check = CheckContain(Faces_shift, shifted)
                    
                    if check:
                        #print("For " + str(x) + " the corresponding gamma is " + str(gamma) + 
                        #      ", with z-distance " + str(  gamma[2] - x[2]) +".")
                        gamma_return = gamma    #The correct gamma for said x
                        end_time = time.time()
                        print("Run time to check at x=" + str(x) + " is " + str(end_time - start_time))
                        return [True, gamma_return]     # return value if a
                                                        # gamma is found
        if r == 0: continue
    if check == False: print("For " + str(x) +" with N="+str(N)+ " and lamb="+ str(lamb)+
                             ", there is no corresponding gamma.")
                                                # In case no gamma is found
                                                # prints message
    end_time = time.time()
    print("Run time to check at x="+ str(x)+ " is "  + str(end_time - start_time))
    return [check, None]       # return value if no gamma is found


def NewInclCheck(lamb, N, init, known):         # checks whether the set inclusions
                                                # hold for all x in the set
                                                # D^N[V] intersected with Gamma
    start_time = time.time()
    D = np.array([lamb, lamb, lamb**2])         #The dilation in each coordinate
    DV_N = N_Dilat(lamb, N)
    W = FundemFaceIter(lamb, N, init)           #The 'z'-faces of the V(N,init)
                                                # set
    end_time = time.time()
    print("Finished preparatory computations after "+ str(end_time-start_time))
    start_time = time.time()
    gammaTox_lst = [ ]
    
    for x in DV_N:
        tru_val = False
        K_shift = ShiftFaces(known, x)
        x_output = NewCheck_x(x, lamb, N, K_shift, W)
        tru_val = x_output[0]
        if tru_val == False: break
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
    
            
def CheckContain(lst_set1, lst_set2):    #determines whether intervals defined
                                         # by lst_set1 by is superset of 
                                         # lst_set2
    
    time_start_tot = time.time()
    xy_lst = XY_Range(lst_set2)
    ind_lst = [ FindXY(lst_set1, xy_tup) for xy_tup in xy_lst ]
    for ind in range(0, len(lst_set2), 2):
        xy_tup =  ( lst_set2[ind][0], lst_set2[ind][1]  ) 
        j = FindXY(lst_set1, xy_tup)
        if j == None: 
            tru_val=False
            break
        tru_val = lst_set1[j][2] <= lst_set2[ind][2] and lst_set1[j+1][2] >= lst_set2[ind+1][2]
        if tru_val == False: break
    time_end_tot = time.time()
    print("Run time for the check function is " + str(time_end_tot - time_start_tot))
    return tru_val

def FindXY(lst, tup):                       # find first index of element with
                                            # 'xy' values equal to tup
    for ind in range(len(lst) ):
        if lst[ind][0]== tup[0] and lst[ind][1]==tup[1]:
            return ind
    return None

#####################################

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
        x_output = NewCheck_x(x, lamb, N, K_shift, W)
        tru_val = x_output[0]
        gammaTox_lst.append( [x, x_output[1] ] )
    end_time = time.time()
    print("Run time for the program in this construction is " + str(end_time-start_time)+".")
    return [tru_val, gammaTox_lst]