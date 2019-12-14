import crlibm
import struct
import math
import secrets
import numpy as np
import gmpy2
import copy

class Snapping_Mechanism:
    """ Implementation of the Snapping Mechanism from Mironov (2012)

    This class presents an implementation of the Snapping Mechanism developed in Mironov (2012)
    (paper found at http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.366.5957&rep=rep1&type=pdf). The mechanism
    and associated proofs are described in section 5.2 of the paper.

    The mechanism was developed in order to avoid vulnerabilities in the Laplace Mechanism due to imprecise
    floating-point calculations. The 'snapping' can be thought of (loosely) as rounding the final output (input + laplace noise) such that
    potential output values are spaced slightly more than lambda apart.

    The mechanism itself is a little different than many others in that it takes the non-private statistic as an argument
    (it could theoretically be rewritten to take the raw data as an argument as well). That is, the noise is a function of
    the non-private statistic (among other things).

    Big thanks to Daniel Alabi and Ira Globus-Harris for their help. Daniel helped point me in the right direction with regards
    to rounding via IEEE representation manipulation. Ira came up with the main idea we eventually used for the IEEE
    manipulation, wrote the first implementation of it, and was very helpful in general with ideas and feedback on various
    piece of the implementation.

    Attributes:
        mechanism_input (numeric): Raw (non-private) value of statistic for which you want a private release
        sensitivity (numeric): Sensitivity of function generating mechanism input
        epsilon (numeric): DP-epsilon
        accuracy (numeric): Desired accuracy level
        alpha (numeric): Desired alpha level for accuracy calculations
        gamma (numeric): Desired upper bound on the probability that clamping bounds bind
        min_B (numeric): Maximum possible value of the non-private statistic (corresponds to minimum viable value of B)
        B (numeric): Clamping bound
        _B_scaled (numeric): Clamping bound scaled by function sensitivity (used only by class internals)
        precision (numeric): Number of bits of precision used for arithmetic operations
        epsilon_prime (numeric): Epsilon value used for parameterization of Laplace within the Snapping Mechanism
        Lambda_prime (numeric): Level at which output is rounded
        _Lambda_prime_scaled (numeric): Lambda_prime scaled by function sensitivity (used only by class internals)
        _m (numeric): Power to which 2 is raised to get _Lambda_prime_scaled (used only by class internals)

    Methods:
        get_snapped_noise: Returns epsilon-DP noise in according with snapping mechanism

    Example:
        snap_object = Snapping_Mechanism(mechanism_input = 50, sensitivity = 0.012, min_B = 200, epsilon = 0.001, accuracy = None)
        snapped_noise = snap_object.get_snapped_noise()
    """

    def __init__(self, mechanism_input, sensitivity, min_B, epsilon = None, accuracy = None, alpha = 0.05, gamma = 0.05):
        if gmpy2.get_max_precision() < 118:
            raise ValueError('Software does not have access to sufficient precision to use the Snapping Mechanism')

        self.mechanism_input = mechanism_input
        self.sensitivity = sensitivity
        self.alpha = alpha # needs to be set before accuracy
        self.gamma = gamma
        self.epsilon = epsilon
        self.accuracy = accuracy
        if self.epsilon is None and self.accuracy is None:
            raise ValueError('Either epsilon or accuracy must be provided.')
        if self.epsilon:
            k = (2 + 24*2**-52) / (self.epsilon - 2**-117)
        else:
            epsilon_lb = math.log(1/self.alpha) * (self.sensitivity / self.accuracy)
            k = (2 + 24*2**-52) / (epsilon_lb - 2**-117)
        self.B = min_B + k/2 * (1 + 2*math.log(1/self.gamma))

        # ensure that B*precision <= 2^-52
        if self.B <= 2**66:
            self.precision = 118
        else:
            # find smallest greater power of two (k is the power to which two is raised)
            _, k = self._get_smallest_greater_or_eq_power_of_two(self.B)

            # add to precision to ensure that B*precision <= 2^-52
            extra_precision = k - 66
            self.precision = 118 + extra_precision

        if self.epsilon is None:
            self.epsilon = self._get_epsilon()
        elif self.accuracy is None:
            self.accuracy = self._get_accuracy()
        self._B_scaled,  self.epsilon_prime, self.Lambda_prime, self._Lambda_prime_scaled, self._m = self._parameter_setup()
        if self.epsilon_prime < 0:
            raise ValueError('The desired accuracy/epsilon results in epsilon_prime < 0. Please choose a smaller accuracy or larger epsilon.')

    def _double_to_bin(self, x):
        """
        Converts a numeric variable (int, float, etc.) to its IEEE 64-bit representation.
        Representation consists of sign (length 1), exponent (length 11), and mantissa (length 52).

        Parameters:
            x (numeric): Number for which 64-bit representation will be returned

        Return:
            String: 64-bit representation of x
        """
        return(bin(struct.unpack('!Q', struct.pack('!d', x))[0])[2:].zfill(64))

    def _bin_to_double(self, binary):
        """
        Converts an IEEE 64-bit representation to a 64-bit float (double).
        Representation consists of sign (length 1), exponent (length 11), and matissa (length 52).

        Parameters:
            binary (String): 64-bit representation of x

        Return:
            float64/double: Number corresponding to 64-bit representation
        """
        return(struct.unpack('!d',struct.pack('!Q', int(binary, 2)))[0])

    def _clamp(self, x, B):
        """
        If abs(x) > abs(B), clamps x towards 0 such that abs(x) = abs(B).

        Parameters:
            x (numeric): Number to be clamped
            B (numeric): Bounds to which x should be clamped

        Return:
            numeric: Clamped version of x
        """
        if (x < -abs(B)):
            return(-abs(B))
        elif (x > abs(B)):
            return(abs(B))
        return(x)

    def _get_ieee_representation(self, binary):
        """
        Get IEEE representation of floating point number from 64-bit binary representation.
        """
        sign = binary[0]
        exponent = binary[1:12]
        mantissa = binary[12:]
        return(sign, exponent, mantissa)

    def _get_smallest_greater_or_eq_power_of_two(self, _lambda):
        """
        Gets closest power of two that is >= _lambda.

        Parameters:
            _lambda (numeric): Argument to laplace noise

        Return:
            numeric: Smallest power of two that is >= _lambda -- we call this Lambda
            int: m such that Lambda = 2^m
        """
        # get IEEE representation of x
        binary = self._double_to_bin(_lambda)
        sign, exponent, mantissa = self._get_ieee_representation(binary)

        # return smallest power of two >= x
        if all(bit == '0' for bit in mantissa):
            return(_lambda, int(exponent, 2)-1023)
        else:
            # create mantissa of all zeros and incremented exponent, then use these to create smallest power of two >= x
            zero_mantissa = ''.zfill(len(mantissa))
            exponent_plus_one = bin(int(exponent, base = 2) + 1)[2:].zfill(len(exponent))
            return(self._bin_to_double(sign + exponent_plus_one + zero_mantissa), int(exponent_plus_one, 2)-1023)

    def _divide_by_power_of_two(self, sign, exponent, mantissa, power):
        """
        Divide by power of two and return updated exponent.
        """
        exponent_num = int(exponent, 2) - power
        exponent = bin(exponent_num)[2:].zfill(11) # use [2:] to remove 0b from beginning of binary representation
        return(sign, exponent, mantissa)

    def _round_to_nearest_int(self, sign, exponent, mantissa):
        """
        Round the IEEE representation to the nearest integer.
        This proceeds by different cases for when the unbiased exponent >= 0 vs < 0.
        """

        # generate numeric version of unbiased exponent
        unbiased_exponent_num = int(exponent, 2) - 1023

        if unbiased_exponent_num >= 52:
            return(sign, exponent, mantissa)
        elif unbiased_exponent_num >= 0:
            '''IEEE_rep >= Lambda'''
            # get elements of mantissa that represent integers
            # (after being multiplied by 2^unbiased_exponent_num)
            mantissa_subset = mantissa[0:unbiased_exponent_num]

            # check to see if mantissa needs to be rounded up or down
            if mantissa[unbiased_exponent_num] == '1':
                '''mantissa needs to be rounded up'''
                # if integer part of mantissa is all 1s, rounding needs to be reflected
                # in the exponent instead
                if mantissa_subset == '1' * len(mantissa_subset):
                    mantissa = ''.ljust(52, '0')[0:52]
                    exponent_num = int(exponent, 2) + 1
                    exponent = bin(exponent_num)[2:].ljust(11, '0')
                else:
                    # if integer part of mantissa not all 1s, just increment mantissa
                    mantissa_subset_inc = bin(int(mantissa_subset, 2) + 1)[2:].zfill(len(mantissa_subset))
                    mantissa = mantissa_subset_inc.ljust(52, '0')[0:52]
            elif mantissa[unbiased_exponent_num] == '0':
                '''mantissa needs to be rounded down'''
                mantissa = mantissa_subset.ljust(52, '0')[0:52]
        elif unbiased_exponent_num < 0:
            '''IEEE_rep < Lambda'''
            if unbiased_exponent_num == -1:
                # round IEEE_rep to 1
                mantissa = ''.ljust(52, '0')
                exponent = '0'.ljust(11, '1')
            elif unbiased_exponent_num < -1:
                # round IEEE_rep to 0
                mantissa = ''.ljust(52, '0')
                exponent = ''.ljust(11, '0')

        return(sign, exponent, mantissa)

    def _multiply_by_power_of_two(self, sign, exponent, mantissa, power):
        """
        Multiply by power of two and return updated exponent.
        """
        if exponent == '0'*len(exponent):
            return(sign, exponent, mantissa)
        else:
            exponent_num = int(exponent, 2) + power
            exponent = bin(exponent_num)[2:].zfill(11) # use [2:] to remove 0b from beginning of binary representation
            return(sign, exponent, mantissa)

    def _get_closest_multiple_of_Lambda(self, x, m):
        """
        Rounds x to the closest multiple of Lambda, resolving ties toward positive infinity.

        Parameters:
            x (numeric): Number to be rounded to closest multiple of Lambda
            m (int): Integer such that Lambda = 2^m

        Return:
            numeric: x rounded to nearest multiple of Lambda
        """
        sign, exponent, mantissa = self._get_ieee_representation(self._double_to_bin(x))
        sign_a, exponent_a, mantissa_a = self._divide_by_power_of_two(sign, exponent, mantissa, m)
        sign_b, exponent_b, mantissa_b = self._round_to_nearest_int(sign_a, exponent_a, mantissa_a)
        sign_c, exponent_c, mantissa_c = self._multiply_by_power_of_two(sign_b, exponent_b, mantissa_b, m)
        return(self._bin_to_double(str(sign_c) + str(exponent_c) + str(mantissa_c)))

    def _sample_from_uniform(self):
        """
        Samples from uniform (0,1) as described by Mironov.
        Meant to sample floating points proportional to their unit of least precision.

        Return:
            numeric: Sample from Unif(0,1)
        """

        # generate draw from geometric distribution
        # run through entire range to protect against timing attacks
        # TODO: if geom produces no successes in 1024 attempts, we currently set it as if its last attempt was a
        #       success -- should think more about this
        geom = 0
        for i in range(1, 1025):
            if (secrets.randbits(1) == 1):
                if geom == 0:
                    geom = copy.deepcopy(i)
        if geom == 0:
            geom = 1024

        # generate draw from uniform as described in Mironov
        u_star_exponent = bin(-geom + 1023)[2:]
        u_star_mantissa = ''.join([str(secrets.randbits(1)) for i in range(52)])
        u_star_sample = self._bin_to_double('0' + str(u_star_exponent) + str(u_star_mantissa))
        return(u_star_sample)

    def _redefine_epsilon(self, epsilon, B, precision):
        """
        Calculates epsilon for the Laplace inside the Snapping Mechanism such that the
        Snapping Mechanism as a whole is (epsilon, 0)-DP.

        Note:
            This is a work in progress -- maybe could be improved
        """

        eta = gmpy2.mpfr(2**-precision)
        return((epsilon - 2*eta) / (1 + 12*B*eta))

    def _get_accuracy(self):
        """
        Get accuracy as described in notes document
        """
        accuracy = ( (1+12*self.B*2**-self.precision) / (self.epsilon-2*2**-self.precision) ) * (1 + math.log(1 / self.alpha)) * (self.sensitivity)
        return(float(accuracy))

    def _get_epsilon(self):
        """
        Get epsilon as described in notes document
        """
        epsilon = ( (1+12*self.B*2**-self.precision) / (self.accuracy) ) * (1 + math.log(1 / self.alpha)) * (self.sensitivity) + 2*2**-self.precision
        return(float(epsilon))

    def _parameter_setup(self):
        """
        Set up parameters for calculations (noise, accuracy, etc.)
        """
        gmpy2.get_context().precision = self.precision

        '''
        Scale input and bounds to sensitivity 1

        The results in the paper all rely on the function that calculates the statistic in question having sensitivity 1.
        So, we scale the input and bounds such that the function has sensitivity 1, implement the snapping mechanism, and rescale back to
        the original sensitivity.
        '''
        # scale clamping bound by sensitivity
        B_scaled = self.B / self.sensitivity
        epsilon_prime = self._redefine_epsilon(self.epsilon, B_scaled, self.precision)

        # NOTE: this Lambda is calculated relative to lambda = 1/epsilon' rather than sensitivity/epsilon' because we have already
        #       scaled by the sensitivity
        Lambda_prime_scaled, m = self._get_smallest_greater_or_eq_power_of_two(1/epsilon_prime)
        Lambda_prime = Lambda_prime_scaled * self.sensitivity

        return(B_scaled, epsilon_prime, Lambda_prime, Lambda_prime_scaled, m)

    def get_snapped_noise(self):
        """
        Generates noise for Snapping Mechanism

        See section 5.2 of Mironov (2012) for the definition of the snapping mechanism implemented below.
        """
        gmpy2.get_context().precision = self.precision

        # scale mechanism input by sensitivity
        mechanism_input_scaled = self.mechanism_input / self.sensitivity

        # generate random sign and draw from Unif(0,1)
        sign = gmpy2.mpfr(2*secrets.randbits(1) - 1) # using mpfr precision here, even though it's an integer, so that the
                                                                               # precision carries through in later steps
        u_star_sample = self._sample_from_uniform()
        inner_result = self._clamp(mechanism_input_scaled, self._B_scaled) + (sign * 1/self.epsilon_prime * crlibm.log_rn(u_star_sample))

        # perform rounding/snapping
        inner_result_rounded = self._get_closest_multiple_of_Lambda(inner_result, self._m)
        private_estimate = self._clamp(self.sensitivity * inner_result_rounded, self.B) # put private estimate back on original scale
        snapped_noise = private_estimate - self.mechanism_input

        return(snapped_noise)