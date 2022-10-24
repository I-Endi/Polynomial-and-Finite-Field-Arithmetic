##
# 2WF90 Algebra for Security -- Software Assignment 2
# Polynomial and Finite Field Arithmetic
# solve.py
#
#
# Group number:
# group_number 
#
# Author names and student IDs:
# Endi Isuf (author_student_ID_1) 
# Dea Llazo (author_student_ID_2)
# Ilesh Yadav (author_student_ID_3)
# author_name_4 (author_student_ID_4)
##

# Import built-in json library for handling input/output 
import json



def solve_exercise(exercise_location : str, answer_location : str):
    """
    solves an exercise specified in the file located at exercise_location and
    writes the answer to a file at answer_location. Note: the file at
    answer_location might not exist yet and, hence, might still need to be created.
    """
    
    # Open file at exercise_location for reading.
    with open(exercise_location, "r") as exercise_file:
        # Deserialize JSON exercise data present in exercise_file to corresponding Python exercise data 
        exercise = json.load(exercise_file)
        

    ### Parse and solve ###

    # Check type of exercise
    if exercise["type"] == "polynomial_arithmetic":
        # Check what task within the polynomial arithmetic tasks we need to perform
        if exercise["task"] == "addition":
            # Solve polynomial arithmetic addition exercise
            pass
        elif exercise["task"] == "subtraction":
            # Solve polynomial arithmetic subtraction exercise
            pass
        elif exercise["task"] == "multiplication":
            # Solve polynomial arithmetic addition exercise
            pass
        elif exercise["task"] == "long_division":
            # Solve polynomial arithmetic subtraction exercise
            pass
        elif exercise["task"] == "extended_euclidean_algorithm":
            # Solve polynomial arithmetic addition exercise
            pass
        elif exercise["task"] == "irreducibility_check":
            # Solve polynomial arithmetic subtraction exercise
            pass
        elif exercise["task"] == "irreducibility_check":
            # Solve polynomial arithmetic subtraction exercise
            pass
    else: # exercise["type"] == "finite_field_arithmetic"
        # Check what task within the finite field arithmetic tasks we need to perform
        if exercise["task"] == "addition":
            # Solve finite field arithmetic addition exercise
            pass
        elif exercise["task"] == "subtraction":
            # Solve finite field arithmetic addition exercise
            pass
        elif exercise["task"] == "multiplication":
            # Solve finite field arithmetic addition exercise
            pass
        elif exercise["task"] == "division":
            # Solve finite field arithmetic addition exercise
            pass
        elif exercise["task"] == "inversion":
            # Solve finite field arithmetic addition exercise
            pass
        elif exercise["task"] == "primitivity_check":
            # Solve finite field arithmetic addition exercise
            pass
        elif exercise["task"] == "primitive_element_generation":
            # Solve finite field arithmetic addition exercise
            pass


    # Open file at answer_location for writing, creating the file if it does not exist yet
    # (and overwriting it if it does already exist).
    with open(answer_location, "w") as answer_file:
        # Serialize Python answer data (stored in answer) to JSON answer data and write it to answer_file
        # json.dump(answer, answer_file, indent=4)
        pass

### Control Function ###

##

### Main Functions ###
    
## Polynomial Arithmetic ##

#Addition
def poly_addition(poly_f, poly_g, mod):

    poly_f, poly_g = adjust_degree(poly_f, poly_g)
    result = []

    for i in range(len(poly_f)):
        coef = poly_f[i] + poly_g[i]
        coef = coef % mod
        result.append(coef)

    return poly_clean(result)
    


### Helper Functions ###

def adjust_degree(poly_f, poly_g):

    deg_f, deg_g = len(poly_f), len(poly_g)

    while deg_f != deg_g:
        if deg_f > deg_g:
            poly_g.append(0)
        else:
            poly_f.append(0)

        deg_f, deg_g = len(poly_f), len(poly_g)

    return poly_f, poly_g
    
def poly_clean(poly_f):

    while len(poly_f) > 0:

        if poly_f[-1] != 0:
            return poly_f
        else: 
            del poly_f[-1]

    return poly_f


f = [
        2,
        1
    ]
g = [
        0,
        0,
        1
    ]
modulus = 3
print(poly_addition(f, g, 3))