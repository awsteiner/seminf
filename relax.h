#ifndef RELAX_H
#define RELAX_H

#include <iostream>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/err_hnd.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::vector<int> ubvector_int;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/** \brief ODE solution by relaxation 

    1/13/02 - Error codes: \n
    1: Singular matrix - row all 0, in pinvs \n
    2: Singular matrix in routine pinvs \n
    3: Memory error (irow<0) in pinvs \n
    10: Iterations exceeded maximum \n

*/
class relax {

 public:

  /// default 1.0e-9
  double eps;

  /// Number of equations
  int ne;

  /// Number of LHS boundary conditions
  int nb;

  /// Grid size
  int ngrid;

  /// Rearrangement
  ubvector_int indexv;

  /// Maximum number of iterations (default 40)
  int itmax;

  /// Scaling of functions
  ubvector scalv;

  /// coordinate grid
  ubvector *x;

  /// function grid
  ubmatrix *y;

#ifndef DOXYGENP

  ubmatrix *s;

  ubmatrix **c;

#endif

  /// Create an object with ne='tne', nb='tnb', ngrid='tngrid'.
  relax(int tne, int tnb, int tngrid);

  virtual ~relax();

  /// Change uniform (possible non-uniform) grid of size ngrid to a 
  /// uniform grid of size newgrid (sets ngrid=newgrid afterwards). 
  /// Uses simple linear interpolation.
  void regrid(int newne, int newnb, int newgrid);
    
  /// Expand left use linear extrapolation to increase
  /// the grid on the LHS or the RHS. Maybe should consider adding 
  /// versions which simply put the boundary y-values in for all of the
  /// new y-values.
  //    void expandleft(int newpts);

  /// Expand right use linear extrapolation to increase
  /// the grid on the LHS or the RHS. Maybe should consider adding 
  /// versions which simply put the boundary y-values in for all of the
  /// new y-values.
  //    void expandright(int newpts);
  
  /// Calculate derivatives
  int calcderiv(int k, int k1, int k2, int jsf, int is1, int isf);

  /// Solve
  int solve(double conv, double slowc);

  /*
    solve(inte ne, void *pa, int (*difeq)(int, int, int, void *));
    int set_iter(int (*iter)(int, double, double, int *, double *, void *));
  */
  
  /// Provide differential equations
  virtual int difeq(int k, int k1, int k2, int jsf, int is1, int isf)=0;
  
  /// Code to execute every iteration
  virtual int iter(int k, double err, double fac, ubvector_int &kmax,
		   ubvector &ermax);

  void bksub(int jf, int k1, int k2);

  void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, 
	   int jmf, int ic1, int jc1, int jcf, int kc);

  int pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k);
  
};
 
#endif
