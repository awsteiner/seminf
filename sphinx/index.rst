Seminf
======

A code for computing semi-infinite nuclear matter

This documentation was generated from git commit

.. include:: commit.rst
	     
This user's guide documents a code which computes semi-infinite
nucleonic matter in 1-dimension in the Thomas-Fermi approximation.
This documentation is generated based on commit 

This code requires requires the installation of `Boost
<http://www.boost.org>`_, `GSL <http://www.gnu.org/software/gsl>`_
(versions 1.16 and later), `HDF5 <http://www.hdfgroup.org>`_ (versions
1.8.14 and later), the most current version of `O2scl
<https://isospin.roam.utk.edu/static/code/o2scl>`_ .
    
The code appears to work marginally for Skyrme model NRAPR (see
\ref seminf_nr) and RMF model RAPR (see \ref seminf_rel), but
still needs quite a bit of work.

.. Eigen isn't working yet...

.. This code is designed to work with generic vector or matrix types.
.. By default, uBlas vector and matrix types are used and the O2scl
.. dense solver is used. Better performance can be obtained with the
.. <a href="http://eigen.tuxfamily.org">Eigen</a> linear algebra
.. library. If USE_EIGEN is defined, then Eigen vector and matrix
.. types are used and the Eigen QR solver with column pivoting is
.. used.

The surface tension (surface energy per unit area) is

.. math::
   \omega \equiv \int~dz \left[ \varepsilon(z) - 
   \mu_n n_n(z) - \mu_p n_p(z) \right]
   
where :math:`\varepsilon(z)` is the full energy density including
the contributions from gradient terms and the neutron and proton
chemical potentials are independent of :math:`z`. It is often
useful to expand the surface tension in terms of the isospin
asymmetry of the dense phase (the LHS)

.. math:: \omega(\delta) = \omega_0 + \delta^2 \omega_{\delta}

This code was originally written for [Steiner05ia]_

.. toctree::
   :maxdepth: 2

   classdoc
   bib
   
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

   
  
