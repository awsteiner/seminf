/*
  -------------------------------------------------------------------
  
  Copyright (C) 2002-2020, Andrew W. Steiner
  
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
#include <iostream>
#include <cmath>
#include <functional>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#ifdef USE_EIGEN
#include <Eigen/Dense>
#endif

#include <o2scl/table.h>
#include <o2scl/constants.h>
#include <o2scl/part.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/lib_settings.h>
#include <o2scl/ode_it_solve.h>
#include <o2scl/linear_solver.h>
#include <o2scl/interp.h>

#ifdef USE_EIGEN
typedef Eigen::VectorXd si_vector_t;
typedef Eigen::MatrixXd si_matrix_t;
typedef Eigen::MatrixXd::RowXpr si_matrix_row_t;
typedef Eigen::MatrixXd si_sp_matrix_t;
#else
typedef boost::numeric::ublas::vector<double> si_vector_t;
typedef boost::numeric::ublas::matrix<double> si_matrix_t;
typedef boost::numeric::ublas::matrix_row<si_matrix_t> si_matrix_row_t;
typedef boost::numeric::ublas::matrix<double> si_sp_matrix_t;
#endif

typedef std::function<
  int(size_t,si_vector_t &,size_t,si_vector_t &,si_matrix_t &)> jac_funct;
typedef std::function
<int(size_t,const si_vector_t &,si_vector_t &)> mm_funct;
typedef std::function
<double(size_t,double,si_matrix_row_t &)> ode_it_funct;
    
