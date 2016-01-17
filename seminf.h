/*
  -------------------------------------------------------------------
  
  Copyright (C) 2002-2016, Andrew W. Steiner
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include <o2scl/ode_it_solve.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

/** \brief A local version of the ODE solver \ref o2scl::ode_it_solve
    
    \future This should be replaced with exact derivatives.
*/
class ode_it_solve2 :
public o2scl::ode_it_solve<o2scl::ode_it_funct11,ubvector,ubmatrix,
  ubmatrix_row,ubvector,ubmatrix> {

protected:
  
  /** \brief The relative stepsize
   */
  double eps;
  
public:
  
  /** \brief Create a new ODE solver
   */
 ode_it_solve2() : o2scl::ode_it_solve<o2scl::ode_it_funct11,ubvector,ubmatrix,
    ubmatrix_row,ubvector,ubmatrix>() {
    eps=1.0e-9;
  }
  
protected:
  
  /** \brief Compute the derivatives of the LHS boundary conditions

      This function computes \f$ \partial f_{left,\mathrm{ieq}} /
      \partial y_{\mathrm{ivar}} \f$
  */
  virtual double fd_left(size_t ieq, size_t ivar, double x,
			 ubmatrix_row &y) {
    
    double ret, dydx;
    
    h=eps*fabs(y[ivar]);
    if (fabs(h)<1.e-15) h=eps;
    
    y[ivar]+=h;
    ret=(*fl)(ieq,x,y);
    
    y[ivar]-=h;
    ret-=(*fl)(ieq,x,y);
    
    ret/=h;
    return ret;
  }
  
  /** \brief Compute the derivatives of the RHS boundary conditions
	
      This function computes \f$ \partial f_{right,\mathrm{ieq}} /
      \partial y_{\mathrm{ivar}} \f$
  */
  virtual double fd_right(size_t ieq, size_t ivar, double x,
			  ubmatrix_row &y) {

    double ret, dydx;
    
    h=eps*fabs(y[ivar]);
    if (fabs(h)<1.e-15) h=eps;

    y[ivar]+=h;
    ret=(*fr)(ieq,x,y);
    
    y[ivar]-=h;
    ret-=(*fr)(ieq,x,y);
    
    ret/=h;
    return ret;
  }
  
  /** \brief Compute the finite-differenced part of the 
      differential equations

      This function computes \f$ \partial f_{\mathrm{ieq}} / \partial
      y_{\mathrm{ivar}} \f$
  */
  virtual double fd_derivs(size_t ieq, size_t ivar, double x,
			   ubmatrix_row &y) {

    double ret, dydx;
    
    h=eps*fabs(y[ivar]);
    if (fabs(h)<1.e-15) h=eps;

    y[ivar]+=h;
    ret=(*fd)(ieq,x,y);
    
    y[ivar]-=h;
    ret-=(*fd)(ieq,x,y);
    
    ret/=h;
    
    return ret;
  }

};

