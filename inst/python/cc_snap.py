import crlibm
import struct
import math
import secrets
import numpy as np
import gmpy2

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
        B (numeric): Suggested bound on the mechanism_input, corresponds to B from Mironov

    Methods:
        get_snapped_noise: Returns epsilon-DP noise in according with snapping mechanism

    Example:
        snap_object = Snapping_Mechanism(50, 0.012, 0.001, 200)
        snapped_noise = snap_object.get_snapped_noise()
    """

    def __init__(self, mechanism_input, sensitivity, epsilon, B):
        self.mechanism_input = mechanism_input
        self.sensitivity = sensitivity
        self.epsilon = epsilon
        self.B = B

    def _double_to_bin(self, x):
        """
        Converts a numeric variable (int, float, etc.) to its IEEE 64-bit representation.
        Representation consists of sign (length 1), exponent (length 11), and matissa (length 52).

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

    def _get_smallest_greater_power_of_two(self, _lambda):
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

    def _round_up_to_nearest_int(self, sign, exponent, mantissa):
        """
        Round the IEEE representation up to the nearest integer.
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

            '''round mantissa up'''
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

        elif unbiased_exponent_num < 0:
            '''IEEE_rep < Lambda'''
            # round IEEE_rep to 1
            mantissa = ''.ljust(52, '0')
            exponent = '0'.ljust(11, '1')

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

    def _get_closest_larger_multiple_of_Lambda(self, x, m):
        """
        Rounds x to the closest larger multiple of Lambda.

        Parameters:
            x (numeric): Number to be rounded to closest larger multiple of Lambda
            m (int): Integer such that Lambda = 2^m

        Return:
            numeric: x rounded to nearest larger multiple of Lambda
        """
        sign, exponent, mantissa = self._get_ieee_representation(self._double_to_bin(x))
        sign_a, exponent_a, mantissa_a = self._divide_by_power_of_two(sign, exponent, mantissa, m)
        sign_b, exponent_b, mantissa_b = self._round_up_to_nearest_int(sign_a, exponent_a, mantissa_a)
        sign_c, exponent_c, mantissa_c = self._multiply_by_power_of_two(sign_b, exponent_b, mantissa_b, m)
        return(self._bin_to_double(str(sign_c) + str(exponent_c) + str(mantissa_c)))

    def _sample_from_uniform(self):
        """
        Samples from uniform (0,1) as described by Mironov.
        Meant to sample floating points proportional to their unit of least precision.

        Return:
            numeric: Sample from Unif(0,1)
        """

        '''
        Looking for more elegant way to sample from geometric
        '''
        geom = 1
        while (secrets.randbits(1) == 0):
            geom += 1
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

    def _get_laplace_CDF(self, x, _lambda):
        """
        Gets P(Laplace(_lambda) <= x) for mean 0 Laplace

        Parameters:
            x (numeric): Number for which we want the area to the left
            _lambda (numeric): Scale parameter for Laplace distribution

        Return:
            numeric: P(Y <= x) where Y ~ Laplace(_lambda)
        """
        if x <= 0:
            return( (1/2) * math.e**(x/_lambda) )
        else:
            return(1 - (1/2) * math.e**(-x/_lambda) )

    def _F_Z_plus(self, z, _lambda):
        """
        Gets upper bound on F_Z(z), where F_Z is the CDF of |Y'| + Lambda/2.
        See section 6.2 of snapping_implementation_notes for explanation.

        Parameters:
            z (numeric): Number for which we want the area to the left
            _lambda (numeric): Scale parameter for Laplace distribution (represented by Y')

        Return:
            numeric: Upper bound on P(Z <= z) where Z ~ |Laplace(_lambda)| + Lambda/2
        """
        return(1 - math.e**(-1/_lambda * z + 1/2))

    def _F_Z_minus(self, z, _lambda):
        """
        Gets lower bound on F_Z(z), where F_Z is the CDF of |Y'| + Lambda/2.
        See section 6.2 of snapping_implementation_notes for explanation.

        Parameters:
            z (numeric): Number for which we want the area to the left
            _lambda (numeric): Scale parameter for Laplace distribution (represented by Y')

        Return:
            numeric: Lower bound on P(Z <= z) where Z ~ |Laplace(_lambda)| + Lambda/2
        """
        return(1 - math.e**(-1/_lambda * z + 1))

    def _get_accuracy(self, f_D, B, P_l_plus, P_u_minus, _lambda, alpha):
        """
        Gets upper bound on accuracy (lower accuracy number means more accurate estimate).
        See section 6.2 of snapping_implementation_notes for explanation.

        Parameters:
            f_D (numeric): Private estimate (aka mechanism output)
            B (numeric): Snapping Mechanism bound
            P_l_plus (numeric): Upper bound on probability that lower snapping bound is binding
            P_u_minus (numeric): Lower bound on probability that upper snapping bound is binding
            _lambda (numeric): Scale parameter for Laplace distribution
            alpha (numeric): Desired confidence level

        Return:
            numeric: Accuracy guarantee
        """
        if P_l_plus >= alpha:
            return(B + f_D)
        elif self._F_Z_minus(B - f_D, _lambda) >= 1 + P_l_plus - alpha:
            return(_lambda * (1 - math.log(alpha - P_l_plus)))
        elif self._F_Z_minus(B - f_D, _lambda) + P_u_minus >= 1 + P_l_plus - alpha:
            return(B - f_D)
        else:
            return(_lambda * (1 - math.log(alpha - P_l_plus + P_u_minus)))

    def get_snapped_noise(self):
        """
        Generates noise for Snapping Mechanism
        """
        # Set bits of numerical precision
        #
        # The precision is set based on one of two components of the mechanism; the epsilon redefinition or the
        # exact rounding of the natural log, whichever needs more precision.
        #
        # The redefinition of epsilon (to epsilon') relies on the numerical precision and
        # can lead to negative epsilon' values (which are not allowed as input to the laplace mechanism)
        # if epsilon < 2*eta, where eta = 2^-precision
        #
        # Exact rounding, described in section 1.1 of http://www.ens-lyon.fr/LIP/Pub/Rapports/RR/RR2005/RR2005-37.pdf,
        # is an alternative to accurate-faithful calculations (which is what most mathemematical libraries are).
        # If the real-valued number falls between two floating point numbers, accurate-faithful calculations
        # return one of those two floating point numbers and usually returns the one that is closer to the
        # real number. Exact rounding always returns the floating point number that is closer.
        #
        # I don't actually know what precision is needed for exact rounding.
        # I chose 118 because this is the number of bits necessary for exact rounding of the log (see section
        # 2.1 at http://www.ens-lyon.fr/LIP/Pub/Rapports/RR/RR2005/RR2005-37.pdf), but I should clarify this
        # with someone who might know better

        # find necessary precision to ensure epsilon_prime > 0
        eps_precision = abs(self._get_smallest_greater_power_of_two(self.epsilon)[1]) + 2
        precision = max(118, eps_precision)
        gmpy2.get_context().precision = precision

        '''
        Scale input and bounds to sensitivity 1

        The results in the paper all rely on the function that calculates the statistic in question having sensitivity 1.
        So, we scale the input and bounds such that the function has sensitivity 1, implement the snapping mechanism, and rescale back to
        the original sensitivity.
        '''
        # scale clamping bound and mechanism input by sensitivity
        mechanism_input_scaled = self.mechanism_input / self.sensitivity
        B_scaled = self.B / self.sensitivity

        '''
        Construct private estimate

        See section 5.2 of Mironov (2012) for the definition of the snapping mechanism implemented below.
        '''
        # generate random sign and draw from Unif(0,1)
        sign = gmpy2.mpfr(2*secrets.randbits(1) - 1) # using mpfr precision here, even though it's an integer, so that the
                                                                               # precision carries through in later steps
        u_star_sample = self._sample_from_uniform()
        epsilon_prime = self._redefine_epsilon(self.epsilon, B_scaled, precision)
        inner_result = self._clamp(mechanism_input_scaled, B_scaled) + (sign * 1/epsilon_prime * crlibm.log_rn(u_star_sample))

        # NOTE: this Lambda is calculated relative to lambda = 1/epsilon rather than sensitivity/epsilon because we have already
        #       scaled by the sensitivity
        Lambda, m = self._get_smallest_greater_power_of_two(1/epsilon_prime)

        inner_result_rounded = self._get_closest_multiple_of_Lambda(inner_result, m)
        private_estimate = self._clamp(self.sensitivity * inner_result_rounded, self.B) # put private estimate back on original scale
        # snapped_noise = private_estimate - self.mechanism_input

        # # calculate extra quantities needed for accuracy guarantee
        # # see snapping_implementation_notes.pdf for explanation
        # lambda_prime = self.sensitivity / epsilon_prime
        # P_l_plus = self._get_laplace_CDF(x = -self.B - private_estimate + lambda_prime, _lambda = lambda_prime)
        # P_u_minus = 1 - self._get_laplace_CDF(x = self.B - private_estimate - lambda_prime, _lambda = lambda_prime)

        # # TODO: just testing this below -- not sure if correct and will likely need to define outside
        # #       of this function to allow user to set alpha
        # accuracy = self._get_accuracy(f_D = private_estimate,
        #                               B = self.B,
        #                               P_l_plus = P_l_plus,
        #                               P_u_minus = P_u_minus,
        #                               _lambda = lambda_prime,
        #                               alpha = 0.05)


        return(private_estimate)