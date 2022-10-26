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
# Ilesh Yadav (1540025)
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
            answer = poly_addition(exercise['f'], exercise['g'], exercise['integer_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "subtraction":
            # Solve polynomial arithmetic subtraction exercise
            answer = poly_subtraction(exercise['f'], exercise['g'], exercise['integer_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "multiplication":
            # Solve polynomial arithmetic addition exercise
            answer = poly_multiplication(exercise['f'], exercise['g'], exercise['integer_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "long_division":
            # Solve polynomial arithmetic subtraction exercise
            answer_q, answer_r = poly_division(exercise['f'], exercise['g'], exercise['integer_modulus'])
            answer = {"answer-q" : answer_q, "answer-r" : answer_r}

        elif exercise["task"] == "extended_euclidean_algorithm":
            # Solve polynomial arithmetic addition exercise
            answer_x, answer_y, answer_gcd = poly_EEA(exercise['f'], exercise['g'], exercise['integer_modulus'])
            answer = {"answer-a" : answer_x, "answer-b" : answer_y, "answer-gcd" : answer_gcd}

        elif exercise["task"] == "irreducibility_check":
            # Solve polynomial arithmetic subtraction exercise
            answer = poly_irr_check(exercise['f'], exercise['integer_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "irreducible_element_generation":
            # Solve polynomial arithmetic subtraction exercise
            answer = poly_irr_gen(exercise['degree'], exercise['integer_modulus'])
            answer = {"answer" : answer}

    else: # exercise["type"] == "finite_field_arithmetic"
        # Check what task within the finite field arithmetic tasks we need to perform
        if exercise["task"] == "addition":
            # Solve finite field arithmetic addition exercise
            answer = ff_addition(exercise['f'], exercise['g'], exercise['integer_modulus'], exercise['polynomial_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "subtraction":
            # Solve finite field arithmetic addition exercise
            answer = ff_subtraction(exercise['f'], exercise['g'], exercise['integer_modulus'], exercise['polynomial_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "multiplication":
            # Solve finite field arithmetic addition exercise
            answer = ff_multiplication(exercise['f'], exercise['g'], exercise['integer_modulus'], exercise['polynomial_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "division":
            # Solve finite field arithmetic addition exercise
            if exercise['g'] == [0]:
                answer = None
            else:
                answer = ff_division(exercise['f'], exercise['g'], exercise['integer_modulus'], exercise["polynomial_modulus"])
            answer = {"answer" : answer}

        elif exercise["task"] == "inversion":
            # Solve finite field arithmetic addition exercise
            answer = ff_inverse(exercise['f'], exercise['integer_modulus'], exercise['polynomial_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "primitivity_check":
            # Solve finite field arithmetic addition exercise
            answer = ff_primitivity_check(exercise['f'], exercise['integer_modulus'], exercise['polynomial_modulus'])
            answer = {"answer" : answer}

        elif exercise["task"] == "primitive_element_generation":
            # Solve finite field arithmetic addition exercise
            answer = ff_primitivity_gen(exercise['integer_modulus'], exercise['polynomial_modulus'])
            answer = {"answer" : answer}


    # Open file at answer_location for writing, creating the file if it does not exist yet
    # (and overwriting it if it does already exist).
    with open(answer_location, "w") as answer_file:
        # Serialize Python answer data (stored in answer) to JSON answer data and write it to answer_file
        json.dump(answer, answer_file, indent=4)

### Control Function ###

##

### Main Functions ###
    
## Polynomial Arithmetic ##


def poly_addition(poly_f, poly_g, mod):
    """
    Do a polynomial subtraction on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """

    # Make g and f arrays of equal length and initiate the result
    poly_f, poly_g = adjust_degree(poly_f, poly_g)
    result = []

    # Add the coefficients next to the same power together starting from the coefficient next to power of zero
    for i in range(len(poly_f)):
        coef = poly_f[i] + poly_g[i]
        # Reduce the sum
        coef = coef % mod
        result.append(coef)

    # Remove the leading zeros and return
    return poly_clean(result)


def poly_subtraction(poly_f, poly_g, mod):
    """
    Do a polynomial subtraction on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """

    # Make g and f arrays of equal length and initiate the result
    poly_f, poly_g = adjust_degree(poly_f, poly_g)
    result = []

    # Subtract the coefficients next to the same power together starting from the coefficient next to power of zero
    for i in range(len(poly_f)):
        coef = poly_f[i] - poly_g[i]
        # Reduce the result
        coef = coef % mod
        result.append(coef)

    # Remove the leading zeros and return
    return poly_clean(result) 


def poly_multiplication(poly_f, poly_g, mod):
    """
    Do a polynomial subtraction on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """

    # Get the degrees of the polynomials
    deg_f, deg_g = poly_double_deg(poly_f, poly_g)
    # Initiate a polynomial with all zero coefficients of the calculated degree
    result = [0] * (deg_f + deg_g + 1)
    # REMOVE     # Make g and f arrays of equal length
    # REMOVE    poly_f, poly_g = adjust_degree(poly_f, poly_g)

    # Double for loop through the coefs of the polynomials
    for f in range(deg_f + 1):
        for g in range(deg_g + 1):

            # Multiply coef next to power "f" in poly_f and coef next to power "g" in poly_g
            # and add to the coef of the result next to power of "f + g"
            prod = poly_f[f] * poly_g[g]
            # Reduce the sum
            result[f + g] = (result[f + g] + prod) % mod

    # Remove the leading zeros and return
    return poly_clean(result)


def poly_division(poly_f, poly_g, mod):
    """
    Do a polynomial long division on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # If the divisor is 0 then we return undefined
    if poly_g == [0]:
        return None, None
    elif poly_f == [0]:
        return [0], [0]
    # Find the degrees of the polynomials
    deg_f, deg_g = poly_double_deg(poly_f, poly_g)
    # Initiate quotient and reminder where quotient is initially 0
    # and reminder is initially the dividend
    q, r, deg_r = [0] + [0] * (deg_f - deg_g), poly_f, deg_f

    # While degree of remainder >= defree of the divisor that means that we can still divide
    while deg_r >= deg_g and r != [0]:
        # Multiply leading coef of remainder with the inverse of the leading coef of the divisor
        takeaway = (r[-1] * inverse(poly_g[-1], mod)) % mod
        # Add the product to the coef of the quotient at degree of remainder minus
        # degree of divisor since we are calcualting the inverse of the divisor
        q[deg_r - deg_g] = (q[deg_r - deg_g] + takeaway) % mod
        # Make the takeaway a list and multiply it with the divisor
        takeaway = ([0] * (deg_r - deg_g)) + [takeaway]
        prod = poly_multiplication(takeaway, poly_g, mod)
        # Subtract the takeaway from the remainder and update
        r = poly_subtraction(r, prod, mod)
        r = poly_clean(r)
        deg_r = len(r) - 1

    # Return result
    return poly_clean(q), poly_clean(r)


def poly_EEA(poly_f, poly_g, mod):
    """
    Do a polynomial version of extended euclidian algorithm on 'poly_f' and 'poly_g' with the given integer modulus 'mod' which is the prime modulus that
    denotes the coefficient field and return a,b,gcd such that a*poly_f + b*poly_g = gcd(poly_f, poly_g).
    The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Assign the proper values to x, v, y, u
    x, v, y, u = [1], [1], [0], [0]

    # Iterate untill 
    while poly_g != [0]:
        # Divide f by g
        q, r = poly_division(poly_f, poly_g, mod)
        # Assign g to the new f and the remainder to the new g
        poly_f,poly_g = poly_g, r
        # Assign x,y to x_1 y_1 and u, v to x, y
        x1, y1 = x, y
        x, y = u, v
        # u = x1 - q*u, v = y1 - q*v
        u = poly_subtraction(x1, poly_multiplication(q, u, mod), mod)
        v = poly_subtraction(y1, poly_multiplication(q, v, mod), mod)

    # x = x*inverse(lc(g)), y = y*inverse(lc(g))
    x, y = poly_clean(x), poly_clean(y)
    x, y = poly_multiplication(x, [inverse(poly_f[-1], mod)], mod), poly_multiplication(y, [inverse(poly_f[-1], mod)], mod)
    # Find gcd and return
    gcd = poly_multiplication(poly_f, [inverse(poly_f[-1], mod)], mod)
    return x, y, gcd


def poly_irr_check(poly_f, mod):
    """
    Check if the given polynomial 'poly_f' is irreducable or not where 'mod' is the prime modulus that
    denotes the coefficient field. The polynomial is represented as a list of coefficients up to the degree of the polynomial.
    """
    # Assign and check degree
    deg_f = len(poly_f) - 1
    if deg_f < 2:
        return True
    
    # Create t and polynomial g
    t = 1
    poly_g = [0, - 1] + ([0]* (pow(mod, t) - 2)) + [1]

    # Check gcd between f and g
    while poly_gcd(poly_f, poly_g, mod) == [1]:
        # if gcd == 1 then increment t by one and update g
        t += 1
        poly_g = [0, - 1] + ([0]* (pow(mod, t) - 2)) + [1]

    # After the iteration if t is equal to 
    # the degree of f then true else false
    if t == deg_f:
        return True
    else:
        return False

def poly_irr_gen(deg_f, mod):
    """
    Generate an irreducible polynomial given it's degree as 'deg_f' and the prime modulus that denotes the coefficient 
    field as 'mod'. The polynomial is represented as a list of coefficients up to the degree of the polynomial.
    """    
    # Generate random polynomial
    poly_f = poly_clean(poly_random_gen(deg_f, mod))
    # While not irreducible generate another one
    while not poly_irr_check(poly_f, mod):
        poly_f = poly_clean(poly_random_gen(deg_f, mod))

    # Return reducible polynomial
    return poly_f


### Finite Field Arithmetic ###

def ff_addition(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field addition on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Add the polynomials and divide by the polynomial modulus
    poly_sum = poly_addition(poly_f, poly_g, mod)
    q, r = poly_division(poly_sum, poly_h, mod)
    return r

def ff_subtraction(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field subtraction on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Subtract the polynomials and divide by h
    poly_diff = poly_subtraction(poly_f, poly_g, mod)
    q, r = poly_division(poly_diff, poly_h, mod)
    return r

def ff_multiplication(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field multiplication on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Multiply the polynomals and divide by h
    poly_prod = poly_multiplication(poly_f, poly_g, mod)
    q, r = poly_division(poly_prod, poly_h, mod)
    return r

def ff_division(poly_f, poly_g, mod, poly_h):
    """
    Do a finite field division on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Inverse poylnomial g and multiply it by f
    poly_inv = poly_clean(ff_inverse(poly_g, mod, poly_h))
    result = ff_multiplication(poly_f, poly_inv, mod, poly_h)
    return result


def ff_inverse(poly_f, mod, poly_h):
    """
    Find the finite field inverse of the given polynomial, 'poly_f', on the irreducible polynomials 'poly_f' and 'poly_g' where 'mod' and 'poly_h' define 
    the arithmetic field. The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Find the gcd
    x, y, gcd = poly_EEA(poly_f, poly_h, mod)
    # If there exists an inverse then divide x by h and return
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
    # Calculate the order and find the prime divisors
    order = pow(mod, len(poly_f) - 1)
    prime_divisors = find_prime_divisors(order - 1)
    # For every prime divisor raise f to the power (order - 1/ the prime number)
    for prime in prime_divisors:
        # If the exponentiation returns 1 then it is not primitive
        if ff_exp(poly_f, (order - 1) / prime, mod, poly_h) == [1]:
            return False
    
    return True


def ff_primitivity_gen(mod, poly_h):
    """
    Given 'mod' and 'poly_f' which define a finite field, generate a primitive element for this field.
    The polynomials are represented as a list of coefficients up to the degree of the polynomial.
    """
    # Generate a random polynomial of a ceratin degree
    poly_f = poly_clean(poly_random_gen(len(poly_h) - 1, mod))
    # If not primitive generate another one
    while not ff_primitivity_check(poly_f, mod, poly_h):
        poly_f = poly_clean(poly_random_gen(len(poly_h) - 1, mod))

    # Return primitive polynomial
    return poly_f


### Helper Functions ###

# Make both lists equal length
def adjust_degree(poly_f, poly_g):

    # get the degree
    deg_f, deg_g = len(poly_f), len(poly_g)

    # Append 0 to the smaller list till they're equal
    while deg_f != deg_g:
        if deg_f > deg_g:
            poly_g.append(0)
        else:
            poly_f.append(0)

        deg_f, deg_g = len(poly_f), len(poly_g)

    return poly_f, poly_g


# Remove leading zeros in the polynomial
def poly_clean(poly_f):

    while len(poly_f) > 1:

        if poly_f[-1] != 0:
            return poly_f
        else: 
            del poly_f[-1]

    return poly_f

# Reduce the coefficients in the polynomial by mod
def poly_reduce(poly_f, mod):
    
    for i in range(len(poly_f)):
        poly_f[i] %= mod
    
    return poly_f

# Find the degree of 2 polynomials
def poly_double_deg(poly_f, poly_g):
    
    return len(poly_f) - 1, len(poly_g) - 1

# Calculate an integer inverse
def inverse(x, mod):
    return pow(x, -1, mod)

# Calculate the gcd of two polynomials without finding x and y
def poly_gcd(poly_f, poly_g, mod):

    while poly_g != [0]:
        q, r = poly_division(poly_f, poly_g, mod)
        poly_f, poly_g = poly_g, r

    gcd = poly_multiplication(poly_f, [inverse(poly_f[-1], mod)], mod)
    return gcd

# Generate a random polynomial of degree 'deg_f'
def poly_random_gen(deg_f, mod):

    poly_f = []
    for i in range(deg_f + 1):
        poly_f = poly_f + [randint(0, mod - 1)]

    return poly_f

# Find all prime divisors for number 'a'
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

# Raise a polynomial to the given exponent in finite field defined by mod and poly_h
def ff_exp(poly_f, exp, mod, poly_h):
    
    if exp == 0:
        return [1]

    result = poly_f
    while exp > 1:
        if exp % 2 == 0:
            result = poly_multiplication(result, result, mod)
            exp = exp / 2
        else:
            result = poly_multiplication(result, poly_f, mod)
            exp = exp - 1
        result = poly_division(result, poly_h, mod)[1]

    return result


    


########################################## ASSIGNMENT IS ABOVE. BELOW IS FOR TESTING PURPOSES ONLY ################################################

## Write comments while chekcing for special values

### Correctness Checking ###
import time

def compare(comp_loc, ans_loc, exc_loc):
    with open(comp_loc, "r") as computed_file:
        # Deserialize JSON exercise data present in exercise_file to corresponding Python exercise data 
        computed = json.load(computed_file)
    with open(ans_loc, "r") as ans_file:
        # Deserialize JSON exercise data present in exercise_file to corresponding Python exercise data 
        ans = json.load(ans_file)

    with open(exc_loc, "r") as exc_file:
        # Deserialize JSON exercise data present in exercise_file to corresponding Python exercise data 
        exc = json.load(exc_file)

    if "gen" in exc['task']:
        return "Generator"
    if computed == ans:
        return "Correct"
    else:
        return "Incorrect"

for i in range(18):
    diff = time.perf_counter()
    solve_exercise("Realistic/Exercises/exercise" + str(i) + ".json", "Solved/Realistic/answer" + str(i) + ".json" )
    diff = time.perf_counter() - diff
    print("Realistic " + str(i) + ": " + str(diff))

for i in range(18):
    diff = time.perf_counter()
    solve_exercise("Simple/Exercises/exercise" + str(i) + ".json", "Solved/Simple/answer" + str(i) + ".json" )
    diff = time.perf_counter() - diff
    print("Simple " + str(i) + ": " + str(diff))


for i in range(18):
    print("Realistic " + str(i) + ": " + compare("Solved/Realistic/answer" + str(i) + ".json", "Realistic/Answers/answer" + str(i) + ".json", "Realistic/Exercises/exercise" + str(i) + ".json"))
for i in range(18):
    print("Simple " + str(i) + ": " + compare("Solved/Simple/answer" + str(i) + ".json", "Simple/Answers/answer" + str(i) + ".json", "Simple/Exercises/exercise" + str(i) + ".json"))



# i=9
# diff = time.perf_counter()
# solve_exercise("Realistic/Exercises/exercise" + str(i) + ".json", "Solved/Realistic/answer" + str(i) + ".json" )
# diff = time.perf_counter() - diff
# print(str(i) + ": " + str(diff))

# print(str(i) + ": " + compare("Solved/Realistic/answer" + str(i) + ".json", "Realistic/Answers/answer" + str(i) + ".json", "Realistic/Exercises/exercise" + str(i) + ".json"))

# print(ff_primitivity_gen(5, [4,
#         1,
#         2,
#         1
#     ]))

