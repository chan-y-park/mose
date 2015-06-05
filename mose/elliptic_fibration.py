import logging
import math
import cmath
import pdb
import sympy
import numpy
import scipy

import weierstrass as wss

from cmath import pi, exp, phase, sqrt
from numpy import linalg as LA

from branch import BranchPoint, reference_period

from scipy.integrate import odeint
from scipy import interpolate

# temporarily added the following: should check if something is superfluous
from sympy import diff, N, simplify
from sympy import mpmath as mp
from sympy.utilities.lambdify import lambdify
from cmath import exp, pi
from numpy import array, linspace
from numpy.linalg import det
from operator import itemgetter
from misc import path_derivative, data_plot, complexify, id_generator, \
                    rotate_poly


NEGLIGIBLE_BOUND = 0.1**12

class EllipticFibration:
    def __init__(self, w_f, w_g, params, branch_point_charges=None):
        self.params = params

        self.sym_f = sympy.sympify(w_f)
        self.num_f = self.sym_f.subs(self.params)
        self.sym_g = sympy.sympify(w_g)
        self.num_g = self.sym_g.subs(self.params)

        u = sympy.Symbol('u')
        self.f_coeffs = map(complex, sympy.Poly(self.num_f, u).all_coeffs())
        self.g_coeffs = map(complex, sympy.Poly(self.num_g, u).all_coeffs())

        ### rotation of the u-plane
        self.u_rot_phase = 0

        # We will work with the convention that the DSZ matrix is fixed to be
        # the following. Must keep this attribute of the class, as it will be 
        # used when computing the KSWCF for new Kwalls.
        self.dsz_matrix = [[0, -1], [1, 0]]

        ### The following variable is true as long as all the
        ### branch point monodromies are computed using the 
        ### SAME rotation of the x-plane.
        ### It's important that the rotation be the same, because
        ### a rotation may interchange the positions of, say e_1,e_2
        ### resulting in a change of basis in which we express the 
        ### monodromies. As long as we have the same basis for all
        ### monodromies, we are good.
        x_rotation_consistency_check = True

        rotation_n = 10
        zeta = numpy.exp(2 * numpy.pi * 1j / rotation_n)

        for i in range(rotation_n):
            phase = zeta ** (i)
            
            r_num_f = self.sym_f.subs(self.params) / (phase ** 2)
            r_num_g = self.sym_g.subs(self.params) / (phase ** 3)

            r_f_coeffs = map(complex, sympy.Poly(r_num_f, u).all_coeffs())
            r_g_coeffs = map(complex, sympy.Poly(r_num_g, u).all_coeffs())

            self.w_model = wss.WeierstrassModelWithPaths(r_f_coeffs, \
                                                                r_g_coeffs)

            self.u_rot_phase = self.w_model.u_rot_phase
            new_r_f_coeffs = rotate_poly(r_f_coeffs, self.u_rot_phase)
            new_r_g_coeffs = rotate_poly(r_g_coeffs, self.u_rot_phase)

            logging.info('I have rotated the x-plane {} times so far.'\
                                                                    .format(i))
            logging.info('\nf = \n{}\n\ng = \n{}\n'.format(\
                                                numpy.poly1d(new_r_f_coeffs), \
                                                numpy.poly1d(new_r_g_coeffs) \
                                                ))

            branch_point_loci = list(self.w_model.get_discriminant_locus())

            branch_point_monodromies = [
                wss.monodromy_at_point_wmodel(i, self.w_model) 
                for i in range(len(branch_point_loci))
            ]
            
            x_rotation_consistency_check = \
                                    self.w_model.x_rotation_consistency_check

            ### If, after the i-th rotation, we do have consistency,
            ### we can eventually exit the cycle.
            if x_rotation_consistency_check == True:
                self.x_rotation = phase
                break
            else:
                if i == rotation_n - 1:
                    raise ValueError('Rotated a full 2pi angle, and couldnt '+\
                                    'compute the braiding!')
                else:
                    logging.info('\n\
                    *************************************\n\
                    Having trouble tracking root braiding\n\
                            Will rotate the x-plane      \n\
                    *************************************\n')
        
        ### Now we update the f,g of the fibration class:
        ### we must take into account both the x-plane rotation
        ### and the u-plane rotation.
        self.sym_f = (sympy.sympify(w_f) / (phase ** 2))
        self.num_f = self.sym_f.subs(self.params).subs(u, u * self.u_rot_phase)
        self.sym_g = sympy.sympify(w_g) / (phase ** 3)
        self.num_g = self.sym_g.subs(self.params).subs(u, u * self.u_rot_phase)
        self.f_coeffs = map(complex, sympy.Poly(self.num_f, u).all_coeffs())
        self.g_coeffs = map(complex, sympy.Poly(self.num_g, u).all_coeffs())

        branch_point_charges = [
            monodromy_eigencharge(m) for m in branch_point_monodromies
        ]

        ### Here we determine the functional expressions for the e_i
        ### these will be used both in hair evolution and in primary kwall
        ### evolution. For each branch point, this expression is thus used 
        ### 3 times, hence computing it here only once is a great boost.
        x = sympy.Symbol('x')
        eq = sympy.simplify(x ** 3 + self.num_f * x + self.num_g)
        self.sym_roots = sympy.simplify(sympy.solve(eq, x))
        

        ### The following might be useful if some trouble is encountered
        ### with automatic computation of discriminant charges from weierstrass 
        ### module.
        ###
        # if branch_point_charges is None:
        #     # Calculate branch point charges using weierstrss.py
        #     branch_point_monodromies = [
        #         wss.monodromy_at_point_wmodel(i, self.w_model) 
        #         for i in range(len(branch_point_loci))
        #     ]
            
        #     branch_point_charges = [
        #         monodromy_eigencharge(m) for m in branch_point_monodromies
        #     ]
        # ### ONCE THE WEIERSTRASS MODULE WORKS RELIABLY, NEED TO REMOVE THE 
        # ### 'IF' STRUCTURE, THE FOLLOWING PART WON'T BE NEEDED.
        # else:
        #     # THIS IS A TEMPORARY DUMMY ASSIGNMENT
        #     branch_point_monodromies = [
        #         [[1, -2],[0, 1]] for i in range(len(branch_point_loci))
        #     ]
        ###

        ### Introducing string identifiers to label branch-points.
        ### These will be used when building genealogies of intersection 
        ### points, to compare them and build MS walls accordingly.
        bp_identifiers = [id_generator() 
                          for i in range(len(branch_point_loci))]

        logging.info('\
              \nHere is the full list of branch point data:\
              \n-------------------------------------------\n')
        for i in range(len(branch_point_loci)):
            logging.info(
                '\
                \nBranch point {}\
                \n--------------------\
                \nLocus: {}\
                \nMonodromy: {}\
                \nCharge: {}\
                \nInternal label: {}\
                \n'\
                .format(
                    i,
                    branch_point_loci[i],
                    numpy.array(branch_point_monodromies[i]).tolist(),
                    branch_point_charges[i],
                    bp_identifiers[i],
                )
            )

        
        self.branch_points = [
            BranchPoint(
                locus=branch_point_loci[i],
                positive_charge=branch_point_charges[i], 
                monodromy_matrix=branch_point_monodromies[i],
                identifier=bp_identifiers[i],
                fibration = self,
            )
            for i in range(len(branch_point_loci))
        ]

        ### Now determine the "positive periods" of the holomorphic form
        ### on the pinching (gauge) cycle for each branch point.
        ### This will be used to fix the sign ambiguity on the charges of
        ### primary K walls.    
        for bp in self.branch_points:
            ### It is necessary to grow hair before computing the 
            ### reference period, because the latter will be determined 
            ### by comparison with the period at the tip of the hair.
            ### The hair growth will make use of the pf_matrix of the
            ### elliptic fibration, hence of f, g; therefore this 
            ### should be called only after no more rotations of the 
            ### x-plane are necessary.
            logging.debug('\nGrowing hair for branch point {}\n'\
                                                            .format(bp.count))
            bp.grow_hair()
            logging.debug('\nDetermine positive period for branch point {}\n'\
                                                            .format(bp.count))
            bp.determine_positive_period(reference_period(bp))
    
    def pf_matrix(self, z):
        f = numpy.poly1d(self.f_coeffs)
        g = numpy.poly1d(self.g_coeffs)
        f_p = f.deriv()
        g_p = g.deriv()
        f_p_p = f_p.deriv()
        g_p_p = g_p.deriv()

        def M10(z): 
            return (\
                18432.0 * (f(z) ** 2) * (f_p(z) ** 2) * g_p(z) \
                - 12.0 * f(z) * (1792.0 * g(z) * (f_p(z) ** 3) \
                - 2560.0 * (g_p(z) ** 3)) \
                + 8192.0 * (f(z) ** 3) \
                * (g_p(z) * f_p_p(z) - f_p(z) * g_p_p(z)) \
                + 432.0 * g(z) \
                * (128.0 * g(z) * g_p(z) * f_p_p(z) \
                + (-4.0 * f_p(z)) \
                * ( 16.0 * (g_p(z) ** 2) + 32.0 * g(z) * g_p_p(z))) \
                ) \
                / (16 * (-64.0 * (f(z) ** 3) - 432.0 * (g(z) ** 2)) \
                * (-48.0 * g(z) * f_p(z) + 32.0 * f(z) * g_p(z))) 
        def M11(z):
            return \
            (192.0 * (f(z) ** 2) * f_p(z) + 864.0 * g(z) * g_p(z)) \
            / (-64.0 * (f(z) ** 3) - 432.0 * (g(z) ** 2)) \
            + (16.0 * f_p(z) * g_p(z) + 48.0 * g(z) * f_p_p(z) \
            - 32.0 * f(z) * g_p_p(z)) \
            / (48.0 * g(z) * f_p(z) - 32.0 * f(z) * g_p(z))
        
        return [[0, 1], [M10(z), M11(z)]]



def monodromy_eigencharge(monodromy):
    ### Leave the following, might be useful in the future to speed 
    ### up the code. But might need a revision to work properly
    ###
    # # m = magnetic charge
    # # n = electric charge
    # # we work in conventions of Seiberg-Witten 2, 
    # # with monodromy matrices acting from the LEFT (hence our matrices
    # # are TRANSPOSE os those in SW!), 
    # # and given by:
    # #       1 + 2 n m    &    - m^2       \\
    # #       4 n^2        &    1 - 2 n m

    # # print "this is the monodromy matrix: \n%s" % monodromy
    # # print "this is the monodromy matrix type: %s" % type(monodromy)

    # nm = (monodromy[0,0] - monodromy[1,1]) / 4.0
    # m_temp = math.sqrt(-1.0 * monodromy[0,1])
    # n_temp = 2 * math.sqrt(monodromy[1,0] / 4.0)
    # if nm != 0:
    #     if nm > 0:
    #         n = n_temp
    #         m = m_temp
    #     else:
    #         n = -n_temp
    #         m = m_temp
    # else:
    #     m = m_temp
    #     n = n_temp
    # return (int(m), int(n))

    transpose_m = monodromy.transpose()

    eigen_syst = LA.eig(transpose_m)
    eigen_vals = eigen_syst[0]
    eigen_vects = eigen_syst[1].transpose()

    eigen_val_0 = int(numpy.rint(eigen_vals[0]))
    eigen_val_1 = int(numpy.rint(eigen_vals[1]))
    eigen_vec_0 = numpy.array(eigen_vects[0])[0]
    eigen_vec_1 = numpy.array(eigen_vects[1])[0]

    min_0 = min(abs(eigen_vec_0))
    max_0 = max(abs(eigen_vec_0))
    min_1 = min(abs(eigen_vec_1))
    max_1 = max(abs(eigen_vec_1))
    
    ### A number may be numerically close to zero
    if min_0 < 10**-8:
        min_0 = 0.0
    if max_0 < 10**-8:
        max_0 = 0.0
    if min_1 < 10**-8:
        min_1 = 0.0
    if max_1 < 10**-8:
        max_1 = 0.0

    ### A vector may be of the type (1,0) or (0,1)
    if min_0 == 0.0:
        min_0 = max_0
    if min_1 == 0.0:
        min_1 = max_1

    ### One or both vectors may be (0,0)
    if max_0 == 0.0 and max_1 == 0.0:
        ### both eigenvectors are (0,0)
        logging.info('Cannot determine eigencharges, they seem to be (0,0).')
    elif max_0 == 0.0 and max_1 != 0.0:
        ### the first eigenvector is (0,0), but the second one is not
        norm_eigen_vec_1 = [int(round(x)) for x in list(eigen_vec_1 / min_1)]
        return norm_eigen_vec_1
    elif max_0 != 0.0 and max_1 == 0.0:
        ### the second eigenvector is (0,0), but the first one is not
        norm_eigen_vec_0 = [int(round(x)) for x in list(eigen_vec_0 / min_0)]
        return norm_eigen_vec_0
    else:
        norm_eigen_vec_0 = [int(round(x)) for x in list(eigen_vec_0 / min_0)]
        norm_eigen_vec_1 = [int(round(x)) for x in list(eigen_vec_1 / min_1)]
        n_norm_eigen_vec_1 = [-1 * int(round(x)) \
                                        for x in list(eigen_vec_1 / min_1)]
        if norm_eigen_vec_0 == norm_eigen_vec_1 \
                                    or norm_eigen_vec_0 == n_norm_eigen_vec_1:
            return norm_eigen_vec_0
        else:
            raise ValueError('Cannot compute monodromy eigenvector for:\
                \n'+str(monodromy)+'\ntwo independent eigenvectors!\n')







