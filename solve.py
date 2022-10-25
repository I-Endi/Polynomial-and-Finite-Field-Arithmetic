##
# 2WF90 Algebra for Security -- Software Assignment 2
# Polynomial and Finite Field Arithmetic
# solve.py
#
#
# Group number:
# 34
#
# Author names and student IDs:
# Endi Isuf (1542591) 
# Dea Llazo (1589857)
##

# Import built-in json library for handling input/output 
import json
from random import randint



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
            # if poly_g == [] then ? else div
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


def poly_addition(poly_f, poly_g, mod):
    """
    Do a polynomial subtraction on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """

    poly_f, poly_g = adjust_degree(poly_f, poly_g)
    result = []

    for i in range(len(poly_f)):
        coef = poly_f[i] + poly_g[i]
        coef = coef % mod
        result.append(coef)

    poly_f, poly_g = poly_clean(poly_f), poly_clean(poly_g)
    return poly_clean(result)


def poly_subtraction(poly_f, poly_g, mod):
    """
    Do a polynomial subtraction on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """

    poly_f, poly_g = adjust_degree(poly_f, poly_g)
    result = []

    for i in range(len(poly_f)):
        coef = poly_f[i] - poly_g[i]
        coef = coef % mod
        result.append(coef)

    poly_f, poly_g = poly_clean(poly_f), poly_clean(poly_g)
    return poly_clean(result) 


def poly_multiplication(poly_f, poly_g, mod):
    """
    Do a polynomial subtraction on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """

    deg_f, deg_g = poly_double_deg(poly_f, poly_g)
    result = [0] * (deg_f + deg_g + 1)
    poly_f, poly_g = adjust_degree(poly_f, poly_g)


    for f in range(deg_f + 1):
        for g in range(deg_g + 1):

            prod = poly_f[f]  * poly_g[g]
            result[f + g] += prod

    poly_f, poly_g = poly_clean(poly_f), poly_clean(poly_g)
    return poly_clean(poly_reduce(result, mod))


def poly_division(poly_f, poly_g, mod):
    """
    Do a polynomial long division on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    deg_f, deg_g = poly_double_deg(poly_f, poly_g)
    q, r = [0] * (deg_f + 1), poly_f
    deg_r = len(r) - 1

    while deg_r >= deg_g:
        takeaway = (r[-1] * inverse(poly_g[-1], mod))
        q[deg_r - deg_g] = (q[deg_r - deg_g] + takeaway) % mod
        takeaway = ([0] * (deg_r - deg_g)) + [takeaway]
        takeaway = poly_multiplication(takeaway, poly_g, mod)
        r = poly_subtraction(r, takeaway, mod)
        r = poly_clean(r)
        deg_r = len(r) - 1

    poly_f, poly_g = poly_clean(poly_f), poly_clean(poly_g)
    return poly_clean(q), poly_clean(r)


def poly_EEA(poly_f, poly_g, mod):
    """
    Do a polynomial version of extended euclidian algorithm on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field and return a,b,gcd such that a*poly_f + b*poly_g = gcd(poly_f, poly_g).
    The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    x, v, y, u = [1], [1], [], []

    while poly_g != []:
        q, r = poly_division(poly_f, poly_g, mod)
        poly_f,poly_g = poly_g, r
        x1, y1 = x, y
        x, y = u, v
        u = poly_subtraction(x1, poly_multiplication(q, u, mod), mod)
        v = poly_subtraction(y1, poly_multiplication(q, v, mod), mod)

    x, y = poly_clean(x), poly_clean(y)
    x, y = poly_multiplication(x, [inverse(poly_f[-1], mod)], mod), poly_multiplication(y, [inverse(poly_f[-1], mod)], mod)
    gcd = poly_multiplication(poly_f, [inverse(poly_f[-1], mod)], mod)
    return x, y, gcd


def poly_irr_check(poly_f, mod):
    """
    Check if the given polynomial 'poly_f' is irreducable or not where 'mod' is the prime modulus that
    denotes the coefficient field. The polynomial is represented as a list of coefficients up to the degree of the polynomial.
    """
    deg_f = len(poly_f) - 1
    if deg_f < 2:
        return True
    
    t = 1
    poly_g = [0, mod - 1] + ([0]* (pow(mod, t) - 2)) + [1]

    while poly_gcd(poly_f, poly_g, mod) == [1]:
        t += 1
        poly_g = [0, mod - 1] + ([0]* (pow(mod, t) - 2)) + [1]

    if t == deg_f:
        return True
    else:
        return False

def poly_irr_gen(deg_f, mod):
    """
    Generate an irreducible polynomial given it's degree as 'deg_f' and the prime modulus that denotes the coefficient 
    field as 'mod'. The polynomial is represented as a list of coefficients up to the degree of the polynomial.
    """    
    poly_f = poly_clean(poly_random_gen(deg_f, mod))
    while not poly_irr_check(poly_f, mod):
        poly_f = poly_clean(poly_random_gen(deg_f, mod))

    return poly_f


### Finite Field Arithmetic ###

def ff_addition(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field addition on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    poly_sum = poly_addition(poly_f, poly_g, mod)
    q, r = poly_division(poly_sum, poly_h, mod)
    return r

def ff_subtraction(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field subtraction on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    poly_diff = poly_subtraction(poly_f, poly_g, mod)
    q, r = poly_division(poly_diff, poly_h, mod)
    return r

def ff_multiplication(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field multiplication on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    poly_prod = poly_multiplication(poly_f, poly_g, mod)
    q, r = poly_division(poly_prod, poly_h, mod)
    return r

def ff_division(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field division on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # poly_q, poly_r = poly_division(poly_f, poly_g, mod)
    
    # while poly_r != []:
    #     print(poly_r[0])
    #     poly_q, poly_r = poly_division(poly_addition(poly_r, poly_h, mod), poly_g, mod)
    
    # poly_q, poly_r = poly_division(poly_q, poly_h, mod)
    # return poly_clean(poly_r)
    poly_inv = poly_clean(ff_inverse(poly_g, mod, poly_h))
    result = ff_multiplication(poly_f, poly_inv, mod, poly_h)
    return result


def ff_inverse(poly_f, mod, poly_h):
    """
    Find the finite field inverse of the given polynomial, 'poly_f', on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    x, y, gcd = poly_EEA(poly_f, poly_h, mod)
    if gcd == [1]:
        q, r = poly_division(x, poly_h, mod)
        return r
    else:
        return

def ff_primitivity_check(poly_f, mod, poly_h):
    """
    Check whether the given polynomial, 'poly_f' is a primitive polynomial or not in the finite field defined by 'mod' and 'poly_h'.
    The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    prime_divisors = find_prime_divisors(mod - 1)
    order = pow(mod, len(poly_f) - 1)
    for prime in prime_divisors:
        if ff_exp(poly_f, (order - 1) / prime, mod, poly_h) == [1]:
            return False
    
    return True


def ff_primitivity_gen(mod, poly_h):
    """
    Given 'mod' and 'poly_f' which define a finite field, generate a primitive element for this field.
    The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    poly_f = poly_clean(poly_random_gen(len(poly_h) - 1, mod))
    while not ff_primitivity_check(poly_f, mod, poly_h):
        poly_f = poly_clean(poly_random_gen(len(poly_h) - 1, mod))

    print(poly_irr_check(poly_f, mod))
    return poly_f


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

def poly_reduce(poly_f, mod):
    
    for i in range(len(poly_f)):
        poly_f[i] %= mod
    
    return poly_f

def poly_double_deg(poly_f, poly_g):
    
    return len(poly_f) - 1, len(poly_g) - 1

def inverse(x, mod):
    return pow(x, -1, mod)

def poly_gcd(poly_f, poly_g, mod):

    while poly_g != []:
        q, r = poly_division(poly_f, poly_g, mod)
        poly_f, poly_g = poly_g, r

    gcd = poly_multiplication(poly_f, [inverse(poly_f[-1], mod)], mod)
    return gcd

def poly_random_gen(deg_f, mod):

    poly_f = []
    for i in range(deg_f + 2):
        poly_f = poly_f + [randint(0, mod - 1)]

    return poly_f

def find_prime_divisors(a):

    prime = 2
    result = []
    while prime <= a:
        if a % prime == 0:
            result.append(prime)
            while a % prime == 0:
                a = a / prime
        prime += 1

    return result

def ff_exp(poly_f, exp, mod, poly_h):
    
    if exp == 0:
        return [1]

    result = poly_f
    while exp > 1:
        if exp % 2 == 0:
            result = poly_multiplication(result, result, mod)
            exp /= 2
        else:
            result = poly_multiplication(result, poly_f, mod)
            exp -= 1
        result = poly_division(result, poly_h, mod)[1]

    return result


    




## Testing correctness and time, check special values and write comments

f = [
        0,
        1
    ]
g = []
h = [
        5,
        2,
        1
    ]
modulus = 7
# print(ff_primitivity_gen(13, h))
print(ff_primitivity_gen(13, h))
