/** \file main.h
    \brief File containing user's guide documentation
*/
/** \mainpage 

    \section ug_section User's Guide

    This user's guide documents a code which computes semi-infinite
    nucleonic matter in 1-dimension in the Thomas-Fermi approximation.
    This documentation is generated based on commit 
    \htmlinclude rev.txt

    This code requires 
    <a href="https://isospin.roam.utk.edu/static/code/o2scl">O<sub>2</sub>scl</a>,
    <a href="http://www.gnu.org/software/gsl">GSL</a>,
    <a href="http://www.boost.org">Boost</a>,
    <a href="http://www.hdfgroup.org/HDF5/">HDF5</a>.
    
    The code appears to work marginally for Skyrme model NRAPR (see
    \ref seminf_nr) and RMF model RAPR (see \ref seminf_rel), but
    still needs quite a bit of work.

    \comment
    Eigen isn't working yet...

    This code is designed to work with generic vector or matrix types.
    By default, uBlas vector and matrix types are used and the O2scl
    dense solver is used. Better performance can be obtained with the
    <a href="http://eigen.tuxfamily.org">Eigen</a> linear algebra
    library. If USE_EIGEN is defined, then Eigen vector and matrix
    types are used and the Eigen QR solver with column pivoting is
    used.
    \endcomment

    The surface tension (surface energy per unit area) is
    \f[
    \omega \equiv \int~dz \left[ \varepsilon(z) - 
    \mu_n n_n(z) - \mu_p n_p(z) \right]
    \f]
    where \f$ \varepsilon(z) \f$ is the full energy density including
    the contributions from gradient terms and the neutron and proton
    chemical potentials are independent of \f$ z \f$. It is often
    useful to expand the surface tension in terms of the isospin
    asymmetry of the dense phase (the LHS)
    \f[
    \omega(\delta) = \omega_0 + \delta^2 \omega_{\delta}
    \f]

    This code was originally written for \ref Steiner05ia .
    
    \anchor Steiner05ia Steiner05ia:
    <a href="http://dx.doi.org/10.1016/j.physrep.2005.02.004">A.W. 
    Steiner, M. Prakash, J.M. Lattimer, and P.J. Ellis</a>,
    Phys. Rep. \b 411 (2005) 325.

*/
