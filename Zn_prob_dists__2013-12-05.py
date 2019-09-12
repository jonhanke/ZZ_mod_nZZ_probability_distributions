#########################################################
##
##  Python class for two kinds of probability distributions on Z/nZ
##
##  Written by Jonathan Hanke on 2013-12-05
##
#########################################################

from copy import deepcopy
import fractions

def gcd(a, b):
    """
    My modified GCD froutine to always return a non-negative integer.
    """
    return abs(fractions.gcd(a,b))





def linear_diophantine_solution(a,b):
    """
    Returns an integer tuple (x, y) so that a * x + b * y = 1.
    
    We assume that the intgers a and b are relativiely prime.
    """
    ## SANITY CHECK: Are a and b integers?
    
    ## SANITY CHECK: Are a and b relatively prime?
    if gcd(a,b) != 1:
        raise RuntimeError, "Oops, a and b must be relatively prime!"

    ## Arrange that A and B satisfy A < B
    if a < b:
        A = b
        B = a
    else:
        A = a
        B = b

    ## Perform the repeated Euclid's algorithm
    alpha_list = []
    r = 1
    while r != 0:
        q = A / B
        r = A - B*q
        alpha_list.append(q)
        A = B
        B = r

    ## Assemble the convergents
    num_list = [0, 1]
    denom_list = [1, 0]
    for al in alpha_list:
        num_list.append(num_list[-1] * al + num_list[-2])    
        denom_list.append(denom_list[-1] * al + denom_list[-2])
        
    ## Assemble the solution
    pm = (-1)**len(alpha_list)
    if a > b:
        y = pm * -num_list[-2]
        x = pm * denom_list[-2]
    else:
        x = pm * -num_list[-2]
        y = pm * denom_list[-2]

    ## SANITY CHECK:
    if a*x + b*y != 1:
        raise RuntimeError, "ERROR: The expected solution (x,y) = " + str((x,y)) + " failed to solve a*x + b*y = 1!"

    ## Return the solution
    return (x, y)



        
class Zn:

    def __init__(self, n, a):
        """
        Initializes the number a in the ring Z/nZ, 
        where a is an integer and n in a positive integer.
        """
        ## SANITY CHECK:
        
        ## SANITY CHECK:
        
        ## SANITY CHECK:

        
        ## Store the modulus and representative
        if isinstance(a, int):
            self.__n = n
            self.__elt = a % n
        elif isinstance(a, Zn):
            if (a.modulus() % n == 0):
                self.__n = n
                self.__elt = a.rep() % n
            else:
                raise RuntimeError, "ERROR: The new modulus must be a divisor of the old modulus!"
            
            
    def __repr__(self):
        return str(self.__elt) + " (mod " + str(self.__n) + ")"

    def modulus(self):
        return deepcopy(self.__n)

    def n(self):
        return deepcopy(self.__n)

    def rep(self):
        """
        Return the reduced representative x for this equivalence class 
        satisfying 0 <= x <= self.modulus().
        """
        return deepcopy(self.__elt)


    
    def __eq__(self, other):
        if not isinstance(other, Zn):   
            return False        
        if (self.n() != other.n()):
            return False        
        if (self.rep() != other.rep()):
            return False        
        return True

    def __neq__(self, other):
        return not self.__eq__(other)


    
    
    def __add__(self, other):

        if isinstance(other, Zn):   
            ## SANITY CHECK: Are these in the same ring?
            if (self.n() != other.n()):
                raise TypeError, "The two elements must be in the same ring Z/nZ."
            
            return Zn(self.n(), self.__elt + other.__elt)

        elif isinstance(other, int):
            return Zn(self.n(), self.__elt + other)

        else:
            raise TypeError, "Something's wrong here... in __add__()"

        
    def __sub__(self, other):

        if isinstance(self, Zn) and isinstance(other, Zn):        
            ## SANITY CHECK: Are these in the same ring?
            if (self.n() != other.n()):
                raise TypeError, "The two elements must be in the same ring Z/nZ."
            
            return Zn(self.n(), self.__elt - other.__elt)

        elif isinstance(other, int):
            return Zn(self.n(), self.__elt - other)

    
    def __mul__(self, other):

        if isinstance(other, Zn):        
            ## SANITY CHECK: Are these in the same ring?
            if (self.n() != other.n()):
                raise TypeError, "The two elements must be in the same ring Z/nZ."
            
            return Zn(self.n(), self.__elt * other.__elt)

        elif isinstance(other, int):
            return Zn(self.n(), self.__elt * other)

    
    def __div__(self, other):
        if isinstance(other, Zn):        
            ## SANITY CHECK: Are these in the same ring?
            if (self.n() != other.n()):
                raise TypeError, "The two elements must be in the same ring Z/nZ."
            
            ## SANITY CHECK: Is the denominator a unit?
            if other.is_zero_divisor():
                raise RuntimeError, "Oops, we can only divide by units!"
        
            return Zn(self.n(), self.__elt * other.mult_inverse().rep())

        elif isinstance(other, int):
            return self / Zn(self.n(), other)

        
    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        if isinstance(other, int):
            return Zn(self.n(), other - self.__elt)
        else:
            raise TypeError, "ERROR: The other argument in __rsub__() must be an int."

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rdiv__(self, other):
        if isinstance(other, int):
            return Zn(self.n(), other) / self
        else:
            raise TypeError, "ERROR: The other argument in __rdiv__() must be an int."

    


    def is_unit(self):
        """
        Determines if the element is a unit (i.e. is invertible).
        """
        return (gcd(self.__elt, self.__n) == 1)

    
    def is_zero_divisor(self):
        """
        Determines if the element is a zero-divisor.
        """
        return not self.is_unit()

    
    def mult_inverse(self):
        """
        Returns the multiplicative inverse of the given element.
        """
        n = self.n()
        a = deepcopy(self.__elt)
        x, y = linear_diophantine_solution(a,n)
        return Zn(n, x)



###################################
# In [195]: two = Zn(5, 2)

# In [196]: two + two + two
# Out[196]: 1 (mod 5)

# In [197]: two + two + two + two
# Out[197]: 3 (mod 5)

# In [198]: two * two
# Out[198]: 4 (mod 5)

# In [199]: two / two
# Out[199]: 1 (mod 5)

# In [200]: two - two
# Out[200]: 0 (mod 5)

################################
# In [204]: two = Zn(5, 2)

# In [205]: Zn(5, two)
# Out[205]: 2 (mod 5)

# In [206]: two
# Out[206]: 2 (mod 5)

# In [209]: Zn(1, two)
# Out[209]: 0 (mod 1)
#################################

# In [495]: two = Zn(5, 2)

# In [496]: 1 - two
# Out[496]: 4 (mod 5)

# In [497]: two - 1
# Out[497]: 1 (mod 5)

###############################

class Zn_dist:

    def __init__(self, n, weight_list, weight_ring="Zn"):
        """
        Initializes a probability distribution 
        from a list of weights.

        EXAMPLES:
        > A = Zn_dist(5, range(5), weight_ring="real"); A
        Probability distribution on Z/5Z with density [0.0, 0.1, 0.2, 0.3, 0.4]
        """
        ## SANITY CHECK: Is n a positive integer?
        if (not isinstance(n, int)) or n <= 0:
            raise RuntimeError, "The number n must be a positive integer!"
                    
        ## SANITY CHECK: Are there the correct number of weights?
        if len(weight_list) != n:
            raise TypeError, "The list must have n elements!"

        ## SANITY CHECK: Check that all of the weights are >=0
        for w in weight_list:
            if (weight_ring == "real") and (w < 0):
                raise TypeError, "Oops, all of the weights must be non-negative!"

        ## Coerce the weights into the appropriate ring
        if weight_ring == "Zn":
            new_weight_list = [Zn(n, w)  for w in weight_list]
        elif weight_ring == "real":
            new_weight_list = [1.0 * w  for w in weight_list]
        else:
            raise TypeError, "ERROR: The weight_ring must be 'Zn' or 'real' only!"
            
        ## Normalize and store n and the weights
        self.__weight_ring_str = weight_ring
        self.__n = n
        w_sum = sum(new_weight_list)
        if weight_ring == "Zn":
            if not w_sum.is_unit():
                raise RuntimeError, "Oops, the total weight must be a unit for this to give a probability distrribution!"
        self.__prob_list = [w / w_sum  for w in weight_list]

        
    def __repr__(self):
        return "Probability distribution on Z/" + str(self.n()) + "Z with density " + str(self.density())

    def total_measure(self):
        """
        Returns the total measure of Z/nZ, which should always be one!
        """
        return sum(self.__prob_list)
    
    def density(self):
        return deepcopy(self.__prob_list)

    def weight_ring(self):
        """
        Returns a string specifying the weight ring.
        """
        return deepcopy(self.__weight_ring_str)
    
    def modulus(self):
        return deepcopy(self.__n)

    def n(self):
        return deepcopy(self.__n)

    
    def __add__(self, other):
        """
        Define the sum of two probability distributions.

        EXAMPLES:
            In [40]: A + A
            Out[40]: Probability distribution on Z/5Z with density [0.2, 0.25, 0.25, 0.20000000000000004, 0.1]
        """        
        ## SANITY CHECK: Are these distributions on the same group?
        if self.n() != other.n():
            raise TypeError, "The two distributions must be on the same group Z/nZ."
            
        ## SANITY CHECK: Do they have the same weight ring?
        if self.weight_ring() != other.weight_ring():
            raise TypeError, "The two distributions must have the same weight ring."

        ## Do the convolution
        n = self.n()
        new_dist_list = [sum([self.__prob_list[i] * other.__prob_list[(k-i) % n]  \
                              for i in range(n)])  for k in range(n)]
        return Zn_dist(n, new_dist_list, self.weight_ring())
    

    def __mul__(self, other):
        """
        Define the k-fold sum of this probability distribution with itself,
        where k is a positive integer.

        EXAMPLES:
            In [383]: A = Zn_dist(5, [1,2,3,4,2]); A
            Out[383]: Probability distribution on Z/5Z with density [2 (mod 5), 1 (mod 5), 4 (mod 5), 3 (mod 5), 1 (mod 5)]
            
            In [384]: A
            Out[384]: Probability distribution on Z/5Z with density [2 (mod 5), 1 (mod 5), 4 (mod 5), 3 (mod 5), 1 (mod 5)]
            
            In [385]: 1*A
            Out[385]: Probability distribution on Z/5Z with density [2 (mod 5), 1 (mod 5), 4 (mod 5), 3 (mod 5), 1 (mod 5)]

            In [386]: A == 1 *A
            Out[386]: True
            
            In [388]: 1*A == A
            Out[388]: True

            In [420]: 2*A == A + A
            Out[420]: True
            
            In [421]: 3*A == A + A + A
            Out[421]: True
        """
        if isinstance(other, int) and (other > 0):
            new_dist = self
            for i in range(other-1):
                new_dist = new_dist + self
            return new_dist
        else: 
            raise TypeError, "ERROR: k must be a positive integer!"

    def __rmul__(self, other):
        return self.__mul__(other)
        

    def __eq__(self, other):
        """
        Check if the two probability distributions are equal.
        """
        if self.__n != other.__n:
            return False
        if self.__prob_list != other.__prob_list:
            return False
        return True

    def __neq__(self, other):
        """
        Check if the two probability distributions are not equal.
        """
        return not self.__eq__(other)


    ## These methods may not make sense in general, since we need to have division available!
    
    def mean(self):
        """
        Return the mean (expected) value of the distribution.
        """
        ## SANITY CHECK: Are the weights in Zn?
        if self.__weight_ring_str != "Zn":
            raise TypeError, "Oops, the mean is only defined when the weight ring is Zn!"

        n = self.n()
        return sum([self.__prob_list[i] * Zn(n,i)  for i in range(n)])
        
        
    def variance(self):
        """
        Return the mean (expected) value of the distribution.
        """
        ## SANITY CHECK: Are the weights in Zn?
        if self.__weight_ring_str != "Zn":
            raise TypeError, "Oops, the variance is only defined when the weight ring is Zn!"

        n = self.n()
        m = self.mean()
        return sum([self.__prob_list[i] * Zn(n, (i - m)**2) for i in range(n)])

        
    def moment(self, k):
        """
        Return the k-th moment of the distribution.
        """
        ## SANITY CHECK: Are the weights in Zn?
        if self.__weight_ring_str != "Zn":
            raise TypeError, "Oops, the moments are only defined when the weight ring is Zn!"

        ## SANITY CHECK: Is k a non-negative integer
        if not isinstance(k, int) or (k < 0):
            raise TypeErrror, "Oops, k must be a non-negative integer!"

        n = self.n()            
        return sum([self.__prob_list[i] * Zn(n, i**k)  for i in range(n)])
        

    def moment_gen_func(self, var='t'):
        """
        Returns the moment generating function as a rational function
        """
        numerator_list 
        

    
###########################

# In [58]: D
# Out[58]: Probability distribution on Z/5Z with density [1.0, 0.0, 0.0, 0.0, 0.0]

# In [59]: D + D
# Out[59]: Probability distribution on Z/5Z with density [1.0, 0.0, 0.0, 0.0, 0.0]

# In [60]: D + D + D
# Out[60]: Probability distribution on Z/5Z with density [1.0, 0.0, 0.0, 0.0, 0.0]

# In [61]: D + D + D + D
# Out[61]: Probability distribution on Z/5Z with density [1.0, 0.0, 0.0, 0.0, 0.0]

##########################

# In [38]: A = Zn_dist(5,range(5),"real")

# In [39]: A
# Out[39]: Probability distribution on Z/5Z with density [0.0, 0.1, 0.2, 0.3, 0.4]

# In [44]: B = A

# In [45]: B = B + A; B
# Out[45]: Probability distribution on Z/5Z with density [0.2, 0.25, 0.25, 0.20000000000000004, 0.1]

# In [46]: B = B + A; B
# Out[46]: Probability distribution on Z/5Z with density [0.22499999999999995, 0.19999999999999996, 0.175, 0.175, 0.22499999999999998]

# In [47]: B = B + A; B
# Out[47]: Probability distribution on Z/5Z with density [0.18999999999999997, 0.19, 0.20249999999999996, 0.21499999999999997, 0.20249999999999996]

# In [48]: B = B + A; B
# Out[48]: Probability distribution on Z/5Z with density [0.19999999999999998, 0.205, 0.20375000000000001, 0.19625, 0.195]

# In [49]: B = B + A; B
# Out[49]: Probability distribution on Z/5Z with density [0.20187499999999997, 0.19937499999999997, 0.19749999999999995, 0.19937499999999997, 0.20187499999999994]

# In [50]: B = B + A; B
# Out[50]: Probability distribution on Z/5Z with density [0.1990625, 0.19937500000000002, 0.200625, 0.20093750000000002, 0.2]

# In [51]: B = B + A; B
# Out[51]: Probability distribution on Z/5Z with density [0.20012500000000005, 0.20043750000000002, 0.20012500000000003, 0.19965625000000004, 0.19965624999999998]

# In [52]: B = B + A; B
# Out[52]: Probability distribution on Z/5Z with density [0.20010937499999998, 0.19989062499999996, 0.19982812499999997, 0.19999999999999996, 0.200171875]

# In [53]: B = B + A; B
# Out[53]: Probability distribution on Z/5Z with density [0.19992187500000003, 0.19997656250000004, 0.2000625, 0.2000625, 0.1999765625]

# In [54]: B = B + A; B
# Out[54]: Probability distribution on Z/5Z with density [0.20001953125000002, 0.20003125, 0.2, 0.19996875000000003, 0.19998046875000003]

# In [55]: B = B + A; B
# Out[55]: Probability distribution on Z/5Z with density [0.200004296875, 0.19998867187500002, 0.19998867187500002, 0.20000429687500004, 0.2000140625]

# In [56]: B = B + A; B
# Out[56]: Probability distribution on Z/5Z with density [0.19999433593750002, 0.20000000000000004, 0.20000566406250003, 0.20000351562500002, 0.19999648437500003]

##########################

# In [62]: C = Zn_dist(5, [0,0,0,0,1], "real")

# In [63]: C + C
# Out[63]: Probability distribution on Z/5Z with density [0.0, 0.0, 0.0, 1.0, 0.0]

# In [64]: C + C + C
# Out[64]: Probability distribution on Z/5Z with density [0.0, 0.0, 1.0, 0.0, 0.0]

# In [65]: C + C + C + C
# Out[65]: Probability distribution on Z/5Z with density [0.0, 1.0, 0.0, 0.0, 0.0]

# In [66]: C + C + C + C + C
# Out[66]: Probability distribution on Z/5Z with density [1.0, 0.0, 0.0, 0.0, 0.0]

# In [67]: C + C + C + C + C + C
# Out[67]: Probability distribution on Z/5Z with density [0.0, 0.0, 0.0, 0.0, 1.0]

#############################
# In [68]: E = Zn_dist(5, [1, 0, 0, 0, 1], "real")

# In [69]: E
# Out[69]: Probability distribution on Z/5Z with density [0.5, 0.0, 0.0, 0.0, 0.5]

# In [70]: E + E
# Out[70]: Probability distribution on Z/5Z with density [0.25, 0.0, 0.0, 0.25, 0.5]

# In [71]: E + E + E
# Out[71]: Probability distribution on Z/5Z with density [0.125, 0.0, 0.125, 0.375, 0.375]

# In [72]: E + E + E + E
# Out[72]: Probability distribution on Z/5Z with density [0.0625, 0.0625, 0.25, 0.375, 0.25]

# In [73]: E + E + E + E + E
# Out[73]: Probability distribution on Z/5Z with density [0.0625, 0.15625, 0.3125, 0.3125, 0.15625]

# In [74]: E + E + E + E + E + E
# Out[74]: Probability distribution on Z/5Z with density [0.109375, 0.234375, 0.3125, 0.234375, 0.109375]

# In [75]: E + E + E + E + E + E + E
# Out[75]: Probability distribution on Z/5Z with density [0.171875, 0.2734375, 0.2734375, 0.171875, 0.109375]

# In [76]: E + E + E + E + E + E + E + E
# Out[76]: Probability distribution on Z/5Z with density [0.22265625, 0.2734375, 0.22265625, 0.140625, 0.140625]

# In [77]: E + E + E + E + E + E + E + E + E
# Out[77]: Probability distribution on Z/5Z with density [0.248046875, 0.248046875, 0.181640625, 0.140625, 0.181640625]

# In [78]: E + E + E + E + E + E + E + E + E + E
# Out[78]: Probability distribution on Z/5Z with density [0.248046875, 0.21484375, 0.1611328125, 0.1611328125, 0.21484375]

# In [79]: E + E + E + E + E + E + E + E + E + E + E
# Out[79]: Probability distribution on Z/5Z with density [0.2314453125, 0.18798828125, 0.1611328125, 0.18798828125, 0.2314453125]

# In [80]: E + E + E + E + E + E + E + E + E + E + E + E
# Out[80]: Probability distribution on Z/5Z with density [0.209716796875, 0.174560546875, 0.174560546875, 0.209716796875, 0.2314453125]

# #####################################

# In [82]: F = Zn_dist(5, [1, 1, 0, 0, 1], "real"); F
# Out[82]: Probability distribution on Z/5Z with density [0.3333333333333333, 0.3333333333333333, 0.0, 0.0, 0.3333333333333333]

# In [83]: F + F
# Out[83]: Probability distribution on Z/5Z with density [0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111, 0.2222222222222222]

# In [84]: F + F + F
# Out[84]: Probability distribution on Z/5Z with density [0.2592592592592593, 0.22222222222222224, 0.14814814814814817, 0.14814814814814817, 0.22222222222222224]

# In [85]: F + F + F + F
# Out[85]: Probability distribution on Z/5Z with density [0.2345679012345679, 0.20987654320987653, 0.1728395061728395, 0.1728395061728395, 0.20987654320987656]

# In [86]: F + F + F + F + F
# Out[86]: Probability distribution on Z/5Z with density [0.21810699588477372, 0.20576131687242802, 0.18518518518518517, 0.1851851851851852, 0.20576131687242805]

# In [87]: F + F + F + F + F + F
# Out[87]: Probability distribution on Z/5Z with density [0.2098765432098766, 0.2030178326474623, 0.19204389574759945, 0.19204389574759945, 0.20301783264746232]

# In [88]: F + F + F + F + F + F + F
# Out[88]: Probability distribution on Z/5Z with density [0.2053040695016004, 0.20164609053497945, 0.1957018747142204, 0.1957018747142204, 0.20164609053497945]

# In [89]: F + F + F + F + F + F + F + F
# Out[89]: Probability distribution on Z/5Z with density [0.2028654168571864, 0.20088401158360006, 0.19768327998780674, 0.19768327998780671, 0.2008840115836001]

# In [90]: F + F + F + F + F + F + F + F + F
# Out[90]: Probability distribution on Z/5Z with density [0.2015444800081288, 0.20047756947619771, 0.19875019051973783, 0.1987501905197378, 0.20047756947619771]

# ###################################

