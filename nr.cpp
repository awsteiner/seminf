/*
  -------------------------------------------------------------------
  
  Copyright (C) 2002-2018, Andrew W. Steiner
  
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
#include "seminf.h"

#include <o2scl/eos_had_apr.h>
#include <o2scl/eos_had_schematic.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/eos_had_potential.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/cli.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

/** \brief Semi-infinite nuclear matter for nonrelativistic
    models in the Thomas-Fermi approximation

    Given a energy density functional separated into
    bulk and gradient terms
    \f[
    \varepsilon(n_n,n_p,\nabla n_n,\nabla n_p ) = 
    \varepsilon_{\mathrm{bulk}}(n_n,n_p) +
    \frac{Q_{nn}}{2} \left( \nabla n_n \right)^2 +
    \frac{Q_{pp}}{2} \left( \nabla n_p \right)^2 +
    Q_{np} \nabla n_n \nabla n_p 
    \f]
    one can minimize the thermodynamic potential
    \f[
    \Omega = \int \left(-P\right) dV
    \f]
    to obtain the neutron and proton density profiles. 
    The Euler-Lagrange equations are 
    \f[
    Q_{nn} \nabla^2 n_n + Q_{np} \nabla^2 n_p -
    \frac{1}{2} \left[ 
    \frac{\partial Q_{nn}}{\partial n_n} \left(\nabla n_n\right)^2 +
    \frac{\partial Q_{np}}{\partial n_n} \nabla n_n \nabla n_p +
    \frac{\partial Q_{pp}}{\partial n_n} \left(\nabla n_p\right)^2 +
    \right] = \frac{f_{\mathrm{bulk}}}{n_n} - \mu_{n,0}
    \f]
    and a second equation with \f$ n \leftrightarrow p \f$,
    where \f$ f_{\mathrm{bulk}} = \varepsilon_{\mathrm{bulk}} - 
    T s_{\mathrm{bulk}} \f$ and \f$ \mu_{n,0} \f$ is a constant.
    
    In the
    semi-infinite matter approximation with a single coordinate,
    \f$ x \f$,
    the E-L equations become
    \f{eqnarray}
    Q_{nn} n_n^{\prime \prime} + Q_{np} n_p^{\prime \prime} = 
    \mu_n - \mu_{n,0} \nonumber \\
    Q_{np} n_n^{\prime \prime} + Q_{pp} n_p^{\prime \prime} = 
    \mu_p - \mu_{p,0}
    \f}

    The surface tension is 
    \f{eqnarray}
    \omega &=& \int_{-\infty}^{\infty} 
    \left( f - \mu_{n,0} n_n - \mu_{p,0} n_p \right) dx
    \nonumber \\
    \omega &=& 2 \int_{-\infty}^{\infty} 
    \left( f_{\mathrm{bulk}} - \mu_{n,0} n_n - 
    \mu_{p,0} n_p \right) dx
    \f]
    with the second equality because the bulk and gradient
    contributions to the surface tension are exactly the same.

    The neutron skin thickness in semi-infinite matter is
    \f[
    t \equiv \int_{-\infty}^{\infty} \left( \frac{n_n(z)}{n_{n,\mathrm{LHS}}} - 
    \frac{n_p(z)}{n_{p,\mathrm{LHS}}} \right)~dz
    \f]
    and converting to a spherical geometry gives a factor of 
    \f$ \sqrt{3/5} \f$ so the neutron skin thickness is
    \f$ R_n - R_p = t \sqrt{3/5} \f$ .

    Derive:
    \f[
    \omega(\delta=0) &=& \sqrt{Q_{nn} + Q_{np}} 
    \int_{-\infty}^{\infty} 
    \left[ f_{\mathrm{bulk}}(n_B,\alpha=0)+n_B B \right]^{1/2} dn_B
    \f]

    In neutron-rich matter, the proton density vanishes while the
    neutron density extends out to \f$ x \rightarrow \infty \f$. In
    order to handle this, the algorithm separates the solution into
    two intervals: the first where the proton density is non-zero, and
    the second where the proton density is zero. These two solutions
    are then pasted together to create the final table.
    
    \note The first proton fraction cannot be 0.5, as the 
    algorithm cannot handle later values different from 0.5
    
    12/15/03 - Implemented a new solution method for Qnn!=Qnp for when
    nn_drip becomes > 0. It seems that in this case, it is better not
    to automatically adjust the scale of the x-axis, but use a gradual
    broadening similar to what we have done using the relativistic
    code. This allows calculation down to really low proton fractions!
    We need to see if this is consistent with the analytic code for
    both Qnn<Qnp (hard) and Qnn>Qnp (probably works).

*/
class seminf_nr {
  
protected:

  /// Neutron chemical potential
  double mun;
  
  /// Proton chemical potential
  double mup;
  
  /** \brief Neutron number density on the LHS
   */
  double nn_left;
  
  /** \brief Proton number density on the LHS
   */
  double np_left;
  
  /// Isoscalar gradient term for neutrons
  double qnn;
  
  /// Isoscalar gradient term for protons
  double qpp;
  
  /// Isovector gradient term
  double qnp;
  
  /// Saturation density of isospin-symmetric matter
  double n0half;
  
  /** \brief Desc
   */
  bool show_relax;
  
  /** \brief Desc
   */
  double barn;
  
  /** \brief Desc
   */
  double dmundn;
  
  /** \brief Desc
   */
  double dmundp;
  
  /// Saturation density of matter with current proton fraction
  double nsat;
  
  /// Current proton fraction (in high-density region)
  double protfrac;
  
  /// Proton fraction index
  int pf_index;
  
  /// EOS pointer
  eos_had_base *eos;
  
  /// Neutron
  fermion neutron;
  
  /// Proton
  fermion proton;
  
  /// Thermodynamic functions
  thermo hb;

  /// If true, output relaxation iterations (default false)
  bool relax_file;
  
  /// If true, \f$ Q_{nn}=Q_{np} \f$
  bool qnn_equals_qnp;
  
  /// If true, match to exponential decay on the RHS
  bool low_dens_part;
  
  /// Minimum err
  double min_err;
  
  /// Desc
  double monfact;
  
  /// Desc
  double expo;
  
  /// Neutron drip density
  double nn_drip;
  
  /// Proton drip density
  double np_drip;
  
  /// Desc (default false)
  bool flatden;

  /// Nonlinear equation solver
  mroot_hybrids<mm_funct,si_vector_t,si_matrix_t,jac_funct> nd;
  
  /// Neutron drip density
  double nnrhs;
  
  /// Proton drip density
  double nprhs;
  
  /// Size of the low-density surface region
  double rhs_length;

  /// If true, we have good boundary conditions
  bool goodbound;
  
  /// Initial derivative for densities on LHS (typically negative)
  double first_deriv;

  /** \brief Desc (default 3.0)
   */
  double big;
  
  /** \brief The grid size (default 100)
   */
  int ngrid;

  /// \name Storage for the solution
  //@{
  si_vector_t xstor;
  si_matrix_t ystor;
  //@}
  
  /// \name Store initial guesses for next proton fraction
  //@{
  si_vector_t xg1;
  si_vector_t xg2;
  si_matrix_t yg1;
  si_matrix_t yg2;
  //@}

  /** \brief Desc
   */
  double relax_converge;

  /** \brief Desc
   */
  double rhs_min;

  /** \brief The step size factor for constructing the initial 
      guess (default 7.0)
  */
  double initialstep;

  /// Integrator to compute integrals of final solution
  inte_qagiu_gsl<> gl2;

  /// Store the profiles
  table_units<> tab_prof;
  
  /// Store results across several profiles
  table_units<> tab_summ;

  /// \name Nucleon-nucleon interactions
  //@{
  eos_had_apr eosa;
  eos_had_potential eosg;
  eos_had_skyrme eoss;
  eos_had_schematic eosp;
  //@}

  /// To compute derivatives (currently unused)
  deriv_gsl<> df;
  
  /** \brief The differential equations to be solved
   */
  int derivs(double sx, const si_vector_t &sy, si_vector_t &dydx) {
    double rhsn, rhsp=0.0, det;
    double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

    neutron.n=sy[0];
    proton.n=sy[1];
    if (model=="APR") {
      eosa.gradient_qij2(neutron.n,proton.n,qnn,qnp,qpp,
			 dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
    }
    if (model=="gp") {
      thermo thth;
      if (false) {
	eosg.gradient_qij(neutron,proton,thth,qnn,qnp,qpp,
			  dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
      } else {
	qnn=100.0/hc_mev_fm;
	qpp=100.0/hc_mev_fm;
	qnp=90.0/hc_mev_fm;
	dqnndnn=0.0;
	dqnndnp=0.0;
	dqnpdnn=0.0;
	dqnpdnp=0.0;
	dqppdnn=0.0;
	dqppdnp=0.0;
      }
    }

    eos->calc_e(neutron,proton,hb);

    if (neutron.n<=0.0) {
      if (proton.n<=0.0) {
	dydx[0]=0.0;
	dydx[1]=0.0;
	dydx[2]=0.0;
	dydx[3]=0.0;
      } else {
	dydx[0]=0.0;
	dydx[1]=sy[3];
	dydx[2]=0.0;
	if (model!="APR" && model!="gp") {
	  dydx[3]=(proton.mu-mup)/ qpp;
	} else {
	  dydx[3]=(proton.mu-mup+0.5*dqppdnp*sy[3]*sy[3])/ qpp;
	}
      }
    } else if (proton.n<=0.0) {
      dydx[0]=sy[2];
      dydx[1]=0.0;
      if (model!="APR" && model!="gp") {
	dydx[2]=(neutron.mu-mun)/ qnn;
      } else {
	dydx[2]=(neutron.mu-mun+0.5*dqnndnn*sy[2]*sy[2])/ qnn;
      }
      dydx[3]=0.0;
    } else {
      if (low_dens_part==false) {
	if (model!="APR" && model!="gp") {
	  det=qnn*qpp-qnp*qnp;
	  rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	  rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	  rhsn/=det;
	  rhsp/=det;
	} else {
	  det=qnn*qpp-qnp*qnp;
	  rhsn=(neutron.mu-mun-0.5*dqnndnn*sy[2]*sy[2]+dqnndnp*sy[2]*sy[3]-
		dqnpdnp*sy[3]*sy[3]+0.5*dqppdnn*sy[3]*sy[3])*qpp-
	    (proton.mu-mup+0.5*dqnndnp*sy[2]*sy[2]-dqnpdnn*sy[2]*sy[2]-
	     dqppdnn*sy[2]*sy[3]-0.5*dqppdnp*sy[3]*sy[3])*qnp;
	  rhsp=(proton.mu-mup+0.5*dqnndnp*sy[2]*sy[2]-dqnpdnn*sy[2]*sy[2]-
		dqppdnn*sy[2]*sy[3]-0.5*dqppdnp*sy[3]*sy[3])*qnn-
	    (neutron.mu-mun-0.5*dqnndnn*sy[2]*sy[2]+dqnndnp*sy[2]*sy[3]-
	     dqnpdnp*sy[3]*sy[3]+0.5*dqppdnn*sy[3]*sy[3])*qnp;
	  rhsn/=det;
	  rhsp/=det;
	}
      } else {
	if (model!="APR" && model!="gp") {
	  rhsn=(neutron.mu-mun)/qnn;
	  rhsp=0.0;
	} else {
	  rhsn=(neutron.mu-mun+0.5*dqnndnn*sy[2]*sy[2])/qnn;
	  rhsp=0.0;
	}
      }
    
      dydx[0]=sy[2];
      dydx[1]=sy[3];
      dydx[2]=rhsn;
      dydx[3]=rhsp;
    }
  
    //--------------------------------------------
    // Return sensible results 

    if (!std::isfinite(dydx[2])) {
      cout << "3 not finite." << endl;
      dydx[2]=0.0;
    }
    if (!std::isfinite(dydx[3])) {
      cout << "4 not finite." << endl;
      dydx[3]=0.0;
    }
  
    return 0;
  }

  /** \brief Solve the ODEs when \f$ Q_{nn} = Q_{np} \f$ 
   */
  double solve_qnn_neq_qnp(double lmonfact) {
    bool guessdone;
    int debug=0;
    double dx, xrhs, delta, epsi;
    si_vector_t y(5), dydx(5);
    int ilast=0, interpi;
    int n_eq=5, n_b;
    si_vector_t ox(ngrid);
    si_matrix_t oy(ngrid,n_eq);

    monfact=lmonfact;

    //----------------------------------------------
    // Create object

    if (flatden) {
      n_b=3;
    } else {
      n_b=4;
    }

    ode_it_funct f_derivs=std::bind
      (std::mem_fn<double(size_t,double,si_matrix_row_t &)>
       (&seminf_nr::difeq),this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct f_left=std::bind
      (std::mem_fn<double(size_t,double,si_matrix_row_t &)>
       (&seminf_nr::left),this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct f_right=std::bind
      (std::mem_fn<double(size_t,double,si_matrix_row_t &)>
       (&seminf_nr::right),this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);   
    ode_it_solve<ode_it_funct,si_vector_t,si_matrix_t,si_matrix_row_t,
		 si_vector_t,si_sp_matrix_t> oit;
    oit.niter=60;
    
#ifdef USE_EIGEN
    o2scl_linalg::linear_solver_eigen_colQR
      <Eigen::VectorXd,Eigen::MatrixXd> sol;
    oit.set_solver(sol);
#endif

    si_sp_matrix_t A(ngrid*n_eq,ngrid*n_eq);
    si_vector_t rhs(ngrid*n_eq),dy(ngrid*n_eq);
    
    if (pf_index==1) {

      //----------------------------------------------
      // Construct a guess by shooting

      xstor[1]=0.0;
      ystor(1,1)=nn_left;
      ystor(2,1)=np_left;
      ystor(3,1)=first_deriv;
      ystor(4,1)=first_deriv;
      dx=initialstep/((double)ngrid);
    
      if (debug>0) {
	cout.width(3);
	cout << 1 << " " 
	     << xstor[1] << " " << ystor(1,1) << " " << ystor(2,1) << " "
	     << ystor(3,1) << " " << ystor(4,1) << endl;
      }
    
      guessdone=false;
      for(int i=2;guessdone==false && i<=ngrid;i++) {
	xstor[i]=xstor[i-1]+dx;
      
	for(int j=0;j<4;j++) {
	  y[j]=ystor(j+1,i-1);
	}
	derivs(xstor[i-1],y,dydx);
	for(int j=1;j<=4;j++) {
	  ystor(j,i)=ystor(j,i-1)+dx*dydx[j-1];
	}
      
	if (ystor(1,i)<nn_drip || ystor(2,i)<np_drip) {
	  guessdone=true;
	  ilast=i-1;
	  i=ngrid+10;
	} else if (debug>0) {
	  cout.width(3);
	  cout << i << " " 
	       << xstor[i] << " " << ystor(1,i) << " " << ystor(2,i) << " "
	       << ystor(3,i) << " " << ystor(4,i) << endl;
	}
      }

      // If ilast wasn't set, then nn and np never became smaller than
      // nn_drip or np_drip. In that case, just use the entire guess.
      if (ilast==0) ilast=100;
    
      //----------------------------------------------
      // Arrange guess into relaxation arrays and
      // stretch the solution to grid size=ngrid

      ox[0]=xstor[1];
      ox[ngrid-1]=xstor[ilast];
      for(int i=1;i<ngrid-1;i++) {
	ox[i]=ox[0]+((double)(i))/((double)(ngrid-1))*
	  (ox[ngrid-1]-ox[0]);
      }
      interpi=1;
      for(int i=0;i<ngrid;i++) {
	while(ox[i]>xstor[interpi+1] && interpi<ilast-1) interpi++;
	while(ox[i]<xstor[interpi] && interpi>1) interpi--;
	for(int j=0;j<4;j++) {
	  oy(i,j)=ystor(j+1,interpi)+
	    (ystor(j+1,interpi+1)-ystor(j+1,interpi))*
	    (ox[i]-xstor[interpi])/(xstor[interpi+1]-xstor[interpi]);
	}
	oy(i,4)=ox[ngrid-1];
	ox[i]=((double)(i))/((double)(ngrid-1));
      }
    
    } else {

      //----------------------------------------------
      // Use last solution for guess
    
      for(int i=0;i<ngrid;i++) {
	ox[i]=xg1[i+1];
	for(int j=0;j<n_eq;j++) {
	  oy(i,j)=yg1(j+1,i+1);
	}
      }

    }

    //----------------------------------------------
    // Try ode_it_solve

    low_dens_part=false;

    if (debug>0) {
      for(int i=0;i<ngrid;i+=9) {
	cout.precision(4);
	cout.width(3);
	cout << i << " " << ox[i] << " ";
	for(int j=0;j<5;j++) {
	  cout << oy(i,j) << " ";
	}
	cout << endl;
	cout.precision(6);
      }
    }

    //----------------------------------------------
    // Solve

    oit.tol_rel=relax_converge;
    oit.verbose=1;
    oit.solve(ngrid,n_eq,n_b,ox,oy,f_derivs,f_left,f_right,
	      A,rhs,dy);
    
    if (debug>0) {
      for(int i=0;i<ngrid;i+=9) {
	cout.precision(4);
	cout.width(3);
	cout << i << " " << ox[i] << " ";
	for(int j=0;j<5;j++) {
	  cout << oy(i,j) << " ";
	}
	cout << endl;
	cout.precision(6);
      }
    }
    

    //----------------------------------------------
    // Copy solution for next guess

    for(int i=1;i<=ngrid;i++) {
      xg1[i]=ox(i-1);
      for(int j=1;j<=n_eq;j++) {
	yg1(j,i)=oy(i-1,j-1);
      }
    }
  
    //----------------------------------------------
    // Rescale solution

    for(int i=0;i<ngrid;i++) {
      ox[i]=ox[i]*oy(i,4);
    }

    //----------------------------------------------
    // Store quantites at rhs for later use

    xrhs=ox[ngrid-1];
    if (np_drip==0.0) {
      nnrhs=oy(ngrid-1,0);
      cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
      if (nnrhs<rhs_min) {
	nnrhs=rhs_min;
	cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
      }
    } else {
      nprhs=oy(ngrid-1,1);
      cout << "RHS, nprhs=" << nprhs << endl;
      if (nprhs<rhs_min) {
	nprhs=rhs_min;
	cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
      }
    }

    //----------------------------------------------
    // Store solution and delete seminf_si_relax

    tab_prof.set_nlines(ngrid);
    for(int i=0;i<ngrid;i++) {
      tab_prof.set(0,i,ox[i]);
      for(int j=0;j<5;j++) {
	tab_prof.set(j+1,i,oy(i,j));
      }
    }

    //----------------------------------------------
    // RHS
    
    if (fabs(protfrac-0.5)>1.0e-4) {

      tab_prof.set_nlines(ngrid*2-1);

      if (pf_index==1) {

	//----------------------------------------------
	// Construct a simple linear guess
      
	if (np_drip==0.0) {
	  for(int i=0;i<ngrid;i++) {
	    ox[i]=((double)i)/((double)(ngrid-1));
	    oy(i,0)=(nnrhs-nn_drip)*(1.0-ox[i])+nn_drip;
	    oy(i,1)=-0.08*(1.0-ox[i]);
	    oy(i,2)=0.001;
	  }
	} else {
	  for(int i=0;i<ngrid;i++) {
	    ox[i]=((double)i)/((double)(ngrid-1));
	    oy(i,0)=(nprhs-np_drip)*(1.0-ox[i])+np_drip;
	    oy(i,1)=-0.08*(1.0-ox[i]);
	    oy(i,2)=0.001;
	  }
	}

      } else {
	
	//----------------------------------------------
	// Use last solution for guess
      
	for(int i=0;i<ngrid;i++) {
	  ox[i]=xg2[i+1];
	  for(int j=0;j<3;j++) {
	    oy(i,j)=yg2(j+1,i+1);
	  }
	}

      }

      //----------------------------------------------
      // Solve
  
      n_eq=3;
      n_b=1;

      if (debug>0) {
	for(int i=0;i<ngrid;i+=9) {
	  cout.precision(4);
	  cout.width(3);
	  cout << i << " " << ox[i] << " ";
	  for(int j=0;j<3;j++) {
	    cout << oy(i,j) << " ";
	  }
	  cout << endl;
	  cout.precision(6);
	}
      }
      
      low_dens_part=true;
      cout << "second solve, nn_drip: " << nn_drip << endl;
      
      oit.tol_rel=relax_converge;
      oit.verbose=1;
      oit.solve(ngrid,n_eq,n_b,ox,oy,f_derivs,f_left,f_right,
		A,rhs,dy);

      if (debug>0) {
	for(int i=0;i<ngrid;i+=9) {
	  cout.precision(4);
	  cout.width(3);
	  cout << i << " " << ox[i] << " ";
	  for(int j=0;j<3;j++) {
	    cout << oy(i,j) << " ";
	  }
	  cout << endl;
	  cout.precision(6);
	}
      }

      //----------------------------------------------
      // Copy solution for next guess

      for(int i=1;i<=ngrid;i++) {
	xg2[i]=ox(i-1);
	for(int j=1;j<=n_eq;j++) {
	  yg2(j,i)=oy(i-1,j-1);
	}
      }
  
      //----------------------------------------------
      // Store solution and delete seminf_si_relax
      
      for(int i=0;i<ngrid;i++) {
	tab_prof.set(0,i+ngrid-1,ox(i));
	for(int j=0;j<3;j++) {
	  tab_prof.set(j+1,i+ngrid-1,oy(i,j));
	}
      }
      //delete rel;

      //----------------------------------------------
      // Rearrangment and rescaling

      if (np_drip==0.0) {
	for(int i=ngrid-1;i<2*ngrid-1;i++) {
	  tab_prof.set(5,i,tab_prof.get(3,i));
	  tab_prof.set(4,i,0.0);
	  tab_prof.set(3,i,tab_prof.get(2,i));
	  tab_prof.set(2,i,0.0);
	  tab_prof.set(0,i,tab_prof.get(0,i)*tab_prof.get(5,i)+xrhs);
	}
      } else {
	for(int i=1;i<=ngrid;i++) {
	  tab_prof.set(5,i,tab_prof.get(3,i));
	  tab_prof.set(4,i,tab_prof.get(2,i));
	  tab_prof.set(3,i,0.0);
	  tab_prof.set(2,i,tab_prof.get(1,i));
	  tab_prof.set(1,i,0.0);
	  tab_prof.set(0,i,tab_prof.get(0,i)*tab_prof.get(5,i)+xrhs);
	}
      }
      //    tab_prof.set_nlines(2*ngrid-1);
    }

    if (debug>0) {
      for(int i=0;i<2*ngrid-1;i+=9) {
	if (i==108) i-=8;
	cout.precision(4);
	cout.width(3);
	cout << i << " ";
	for(int j=0;j<6;j++) {
	  cout << tab_prof.get(j,i) << " ";
	}
	cout << endl;
	cout.precision(6);
	if (i==190) i--;
      }
    }

    return 0.0;
  }
  
#ifdef NEVER_DEFINED
  /** \brief Desc
   */
  double solve_qnn_equal_qnp(double lmonfact) {
    bool guessdone, debug=false;
    double dx, y[5], dydx[5], xrhs;
    int ilast=0, interpi;
    si_vector_t sx(3), sy(3);
    int n_eq=3, n_b=2;

    monfact=lmonfact;

    //----------------------------------------------
    // Create object

    rel=new seminf_si_relax(3,2,ngrid);
    rel->itmax=200;
    rel->rp=&rp;

    //----------------------------------------------
    // Construct a guess by shooting

    if (pf_index==1) {
      xstor[1]=0.0;
      ystor(1,1)=nn_left+np_left;
      ystor(2,1)=first_deriv;
      dx=initialstep/((double)ngrid);
    
      guessdone=false;
      for(int i=2;guessdone==false && i<=ngrid;i++) {
	xstor[i]=xstor[i-1]+dx;
      
	sx[1]=(1.0-protfrac)*ystor(1,i-1);
	sx[2]=protfrac*ystor(1,i-1);
	barn=ystor(1,i-1);

	mm_funct qqf=std::bind
	  (std::mem_fn<int(size_t,const si_vector_t &,si_vector_t &)>
	   (&seminf_nr::qnnqnpfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	nd.msolve(2,sx,qqf);
      
	ystor(1,i)=ystor(1,i-1)+dx*ystor(2,i-1);
	ystor(2,i)=ystor(2,i-1)+dx*(neutron.mu-mun)/(qnn);
      
	if (sx[2]<nn_drip+np_drip+1.0e-6 || ystor(2,i)>0.0) {
	  guessdone=true;
	  ilast=i-1;
	  i=ngrid+10;
	} else if (debug) {
	  cout << i << " " << xstor[i] << " " << ystor(1,i) << " " 
	       << ystor(2,i) << endl;
	}
      }
    
      //----------------------------------------------
      // Stretch the solution to the grid ngrid
    
      rel->x[1]=xstor[1];
      rel->x[ngrid]=xstor[ilast];
      for(int i=2;i<ngrid;i++) {
	rel->x[i]=rel->x[1]+((double)(i-1))/((double)(ngrid-1))*
	  (rel->x[ngrid]-rel->x[1]);
      }
      interpi=1;
      for(int i=1;i<=ngrid;i++) {
	while(rel->x[i]>xstor[interpi+1] && interpi<ilast-1) interpi++;
	while(rel->x[i]<xstor[interpi] && interpi>1) interpi--;
	for(int j=1;j<=2;j++) {
	  rel->y(j,i)=ystor(j,interpi)+
	    (ystor(j,interpi+1)-ystor(j,interpi))*
	    (rel->x[i]-xstor[interpi])/(xstor[interpi+1]-xstor[interpi]);
	}
	rel->y(3,i)=(rel->x[ngrid]);
	rel->x[i]=((double)(i-1))/((double)(ngrid-1));
      }

    } else {

      //----------------------------------------------
      // Use last solution for guess
    
      for(int i=1;i<=ngrid;i++) {
	rel->x[i]=xg1[i];
	for(int j=1;j<=n_eq;j++) {
	  rel->y(j,i)=yg1(j,i);
	}
      }
    }

    //----------------------------------------------
    // Solve

    low_dens_part=false;
    int relret1=rel->solve(relax_converge,1.0);
  
    //----------------------------------------------
    // Copy solution for next guess

    for(int i=1;i<=ngrid;i++) {
      xg1[i]=rel->x[i];
      for(int j=1;j<=n_eq;j++) {
	yg1(j,i)=rel->y(j,i);
      }
    }
  
    //----------------------------------------------
    // Rescale solution and set neutron and proton
    // densities:

    tab_prof.set_nlines(ngrid);
    for(int i=1;i<=rel->ngrid;i++) {

      sx[1]=(1.0-protfrac)*rel->y(1,i);
      sx[2]=protfrac*rel->y(1,i);

      barn=rel->y(1,i);

      mm_funct qqf=std::bind
	(std::mem_fn<int(size_t,const si_vector_t &,si_vector_t &)>
	 (&seminf_nr::qnnqnpfun),this,std::placeholders::_1,
	 std::placeholders::_2,std::placeholders::_3);
      nd.msolve(2,sx,qqf);

      tab_prof.set(0,i-1,rel->x[i]*rel->y(3,i));
      tab_prof.set(1,i-1,sx[1]);
      tab_prof.set(2,i-1,sx[2]);
      tab_prof.set(3,i-1,rel->y(2,i));
      tab_prof.set(4,i-1,0.0);
      tab_prof.set(5,i-1,rel->y(3,i));
    }

    //----------------------------------------------
    // Store quantites at rhs for later use:
  
    xrhs=tab_prof[0][ngrid-1];
    if (np_drip==0.0) {
      nnrhs=tab_prof[1][ngrid-1];
      cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
      if (nnrhs<rhs_min) {
	nnrhs=rhs_min;
	cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
      }
    } else {
      nprhs=tab_prof[2][ngrid-1];
      cout << "RHS, nprhs=" << nprhs << endl;
      if (nprhs<rhs_min) {
	nprhs=rhs_min;
	cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
      }
    }

    //----------------------------------------------
    // Delete seminf_si_relax

    delete rel;

    //---------------------------------------------
    // RHS
    
    if (fabs(protfrac-0.5)>1.0e-4) {
      rel=new seminf_si_relax(3,1,ngrid);
      rel->rp=&rp;
    
      if (pf_index==1) {

	//----------------------------------------------
	// Construct a simple linear guess
      
	if (np_drip==0.0) {
	  for(int i=1;i<=ngrid;i++) {
	    rel->x[i]=((double)(i-1))/((double)(ngrid-1));
	    rel->y(1,i)=(nnrhs-nn_drip)*(1.0-rel->x[i])+nn_drip;
	    rel->y(2,i)=-0.08*(1.0-rel->x[i]);
	    rel->y(3,i)=0.001;
	  }
	} else {
	  for(int i=1;i<=ngrid;i++) {
	    rel->x[i]=((double)(i-1))/((double)(ngrid-1));
	    rel->y(1,i)=(nprhs-np_drip)*(1.0-rel->x[i])+np_drip;
	    rel->y(2,i)=-0.08*(1.0-rel->x[i]);
	    rel->y(3,i)=0.001;
	  }
	}
      } else {

	//----------------------------------------------
	// Use last solution for guess
      
	for(int i=1;i<=ngrid;i++) {
	  rel->x[i]=xg2[i];
	  for(int j=1;j<=n_eq;j++) {
	    rel->y(j,i)=yg2(j,i);
	  }
	}
      
      }
    
      //----------------------------------------------
      // Solve
    
      low_dens_part=true;
      int relret2=rel->solve(relax_converge,1.0);
    
      //----------------------------------------------
      // Copy solution for next guess
      
      for(int i=1;i<=ngrid;i++) {
	xg2[i]=rel->x[i];
	for(int j=1;j<=n_eq;j++) {
	  yg2(j,i)=rel->y(j,i);
	}
      }
    
      //----------------------------------------------
      // Store solution, rescale, and delete seminf_si_relax
    
      tab_prof.set_nlines(2*ngrid-1);
      for(int i=1;i<=rel->ngrid;i++) {
	tab_prof.set(5,i-2+ngrid,rel->y(3,i));
	tab_prof.set(3,i-2+ngrid,rel->y(2,i));
	tab_prof.set(4,i-2+ngrid,0.0);
	tab_prof.set(2,i-2+ngrid,0.0);
	tab_prof.set(1,i-2+ngrid,rel->y(1,i));
	tab_prof.set(0,i-2+ngrid,rel->x[i]*rel->y(3,i)+xrhs);
      }
      delete rel;
    }
  
    return 0.0;

  }
#endif
  
  /** \brief Desc
   */
  int qnnqnpfun(size_t sn, const si_vector_t &sx, si_vector_t &sy) {

    neutron.n=sx[1];
    proton.n=sx[2];
    if (sx[1]<0.0 || sx[2]<0.0) return 1;
  
    eos->calc_e(neutron,proton,hb);

    sy[1]=neutron.n+proton.n-barn;
    sy[2]=neutron.mu-proton.mu-mun+mup;

    return 0;
  }

  /** \brief Compute the boundary densities in the case of 
      a neutron drip
  */
  int ndripfun(size_t sn, const si_vector_t &sx, si_vector_t &sy) {
    double pleft, pright, munleft, munright;

    if (sx[0]<0.0 || sx[1]<0.0 || sx[2]<0.0) return 1;

    neutron.n=sx[0];
    proton.n=sx[1];
    sy[0]=proton.n-protfrac*(proton.n+neutron.n);
  
    eos->calc_e(neutron,proton,hb);
    pleft=hb.pr;
    munleft=neutron.mu;

    neutron.n=sx[2];
    proton.n=0.0;
  
    eos->calc_e(neutron,proton,hb);
    pright=hb.pr;
    munright=neutron.mu;

    sy[1]=pleft-pright;
    sy[2]=munleft-munright;

    return 0;
  }

  /** \brief Compute the boundary densities in the case of 
      a proton drip
  */
  int pdripfun(size_t sn, const si_vector_t &sx, si_vector_t &sy) {
    double pleft, pright, mupleft, mupright;

    if (sx[0]<0.0 || sx[1]<0.0 || sx[2]<0.0) return 1;

    neutron.n=sx[0];
    proton.n=sx[1];
    sy[0]=proton.n-protfrac*(proton.n+neutron.n);
  
    eos->calc_e(neutron,proton,hb);
    pleft=hb.pr;
    mupleft=proton.mu;

    neutron.n=0.0;
    proton.n=sx[2];
  
    eos->calc_e(neutron,proton,hb);
    pright=hb.pr;
    mupright=proton.mu;

    sy[1]=pleft-pright;
    sy[2]=mupleft-mupright;

    return 0;
  }

  /** \brief Reformat differential equations for \ref o2scl::ode_it_solve
   */
  double difeq(size_t ieq, double x, si_matrix_row_t &y) {
    if (low_dens_part==false) {
      si_vector_t y2=y, dydx(5);
      derivs(x,y2,dydx);
      if (ieq==0) {
	return dydx[0]*y[4];
      } else if (ieq==1) {
	return dydx[1]*y[4];
      } else if (ieq==2) {
	return dydx[2]*y[4];
      } else if (ieq==3) {
	return dydx[3]*y[4];
      }
      return 0.0;
    } else {
      si_vector_t y2(5), dydx(5);
      y2[0]=y[0];
      y2[1]=0.0;
      y2[2]=y[1];
      y2[3]=0.0;
      y2[4]=y[2];
      derivs(x,y2,dydx);
      if (ieq==0) {
	return dydx[0]*y[2];
      } else if (ieq==1) {
	return dydx[2]*y[2];
      }
      return 0.0;
    }
    O2SCL_ERR("Sanity in difeq().",exc_esanity);
    return 0.0;
  }

  /** \brief LHS boundary conditions for \ref o2scl::ode_it_solve
   */
  double left(size_t ieq, double x, si_matrix_row_t &y) {

    if (low_dens_part==false) {
      
      double delta=nn_left*monfact/exp(big*expo);
      double epsi=delta*(-dmundn+qnn*expo*expo)/
	(dmundp-qnp*expo*expo);
      
      if (epsi*delta<0.0) {
	O2SCL_ERR("Variables 'epsi' and 'delta' have different signs.",
		  o2scl::exc_efailed);
      }
      
      double ret;
      if (ieq==0) {
	ret=y[0]-(nn_left-delta*exp(big*expo));
      } else if (ieq==1) {
	ret=y[1]-(np_left-epsi*exp(big*expo));
      } else if (ieq==2) {
	ret=y[2]+delta*expo*exp(big*expo);
      } else {
	ret=y[3]+epsi*expo*exp(big*expo);
      }
      return ret;
      
    }

    // If low_dens_part is true, then the only boundary condition is
    // that the neutron density on the LHS match that from the RHS of
    // the high-density part of the solution
    return y[0]-nnrhs;
  }

  /** \brief RHS boundary conditions for \ref o2scl::ode_it_solve
   */
  double right(size_t ieq, double x, si_matrix_row_t &y) {
    // If low_dens_part is false, then the right hand boundary condition is
    // that the proton density vanishes.
    if (low_dens_part==false) {
      return y[1];
    }
    // If low_dens_part is true, then the two boundary conditions are that
    // the neutron density goes to the drip density and the derivative
    // of the neutron density vanishes.
    if (ieq==0) {
      return y[0]-nn_drip*(1.0+rhs_adjust);
    }
    return y[1];
  }
  
public:
  
  seminf_nr() {
    rhs_length=2.0;
    big=3.0;
    ngrid=100;
    xg1.resize(ngrid+1);
    xg2.resize(ngrid+1);
    yg1.resize(6,ngrid+1);
    yg2.resize(6,ngrid+1);
    rhs_adjust=0.0;
    model="NRAPR";
    out_file="nr.o2";
    min_err=1.0e-6;
    relax_converge=1.0e-9;
    show_relax=true;
    relax_file=false;
    first_deriv=-2.0e-4;
    rhs_min=1.0e-7;
    monfact=0.002;
  }

  /// Adjustment for RHS boundary
  double rhs_adjust;

  /// Type of EOS model in use
  string model;
  
  /// Name of output file
  string out_file;

  /** \brief Calculate nucleon profiles for a list of proton 
      fractions
  */
  int calc(std::vector<std::string> &sv, bool itive_com) {
    
    double lex1;
    double lex2;
    double rex1;
    double rex2;
    
    double dndnl;
    double dndnr;
    double dndpl;
    double dndpr;
    double dpdnl;
    double dpdnr;
    double dpdpl;
    double dpdpr;
    
    double wd=0.0, wd2=0.0, dtemp;
    double rhsn, rhsp, det;
    double surf, sbulk, sgrad, sssv_drop=0.0;
    vector<double> pflist;
    si_vector_t sx(5), sy(5);
    double sssv_jim=0.0, hns, delta, w0jl=0.0, wdjl=0.0, den, w0;

    // Locations of fixed relative densities
    double xn[3], xp[3];
    
    double thick;
    int npf;
    double lon, lop, hin, hip, sqt;
    double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

    int verbose=1;
    
    // Not used ATM but probably useful later
    //si_vector_t rho, alpha, ebulk, egrad, esurf, thickint;
    //si_vector_t wdint, wd2int;
    //si_vector_t vqnn, vqnp, vqpp;
    
    nd.tol_abs=1.0e-11;
    nd.tol_rel=1.0e-8;
    nd.ntrial=100;
    nd.verbose=0;

    initialstep=7.0;

    //--------------------------------------------
    // Get model parameters

    if (model==((string)"schematic")) {
      eos=&eosp;
    } else if (model==((string)"APR")) {
      eos=&eosa;
    } else if (model==((string)"gp")) {
      eos=&eosg;
    } else {

      eos=&eoss;
      skyrme_load(eoss,model);
      
      // Determine values of Q_{nn}, Q_{pp}, and Q_{np}
      qnn=0.1875*(eoss.t1*(1.0-eoss.x1)-eoss.t2*(1.0+eoss.x2));
      qpp=qnn;
      qnp=0.125*(3.0*eoss.t1*(1.0+eoss.x1/2.0)-
		 eoss.t2*(1.0+eoss.x2/2.0));

    }
  
    //--------------------------------------------
    // Test if Qnn == Qnp

    if (model!="gp" && model!="APR" && fabs(qnn-qnp)<1.0e-7) {
      qnn_equals_qnp=true;
      cout << "Qnn=Qnp" << endl;
      if (true) {
	qnn/=1.01;
	qpp/=1.01;
	cout << "Qnn hacked. Qnn!=Qnp" << endl;
	qnn_equals_qnp=false;
      }
    } else {
      qnn_equals_qnp=false;
      cout << "Qnn!=Qnp" << endl;
    }
    cout << endl;
  
    //--------------------------------------------

    flatden=false;

    // Initialize tables
    tab_prof.set_nlines(ngrid*2);
    tab_prof.line_of_names(((string)"x nn np nnp npp scale rho alpha ")+
			   "ebulk egrad esurf thickint wdint wd2int");
    tab_summ.line_of_names("pf surf thick s_bulk s_grad wd wd2 nn_drip");

    // Allocate storage
    xstor.resize(ngrid+1);
    ystor.resize(5+1,ngrid+1);

    // Initialize particles and EOS
    neutron.init(o2scl_settings.get_convert_units().convert
		 ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    proton.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_proton),2.0);
    neutron.non_interacting=false;
    proton.non_interacting=false;
    eos->set_n_and_p(neutron,proton);
    eos->set_thermo(hb);

    // Compute saturation density and r0 in symmetric matter
    n0half=eos->fn0(0.0,dtemp);
    double r0=cbrt(0.75/pi/n0half);

    // Get list of proton fractions to compute
    pflist.resize(sv.size()-1);
    for(size_t k=0;k<pflist.size();k++) {
      pflist[k]=o2scl::stod(sv[k+1]);
    }
    npf=pflist.size();

    hdf_file hf;
    
    // Loop through proton fractions
    for(pf_index=1;pf_index<=npf;pf_index++) {
      
      protfrac=pflist[pf_index-1];

      //------------------------------------------------------
      // Calculate properties of saturated nuclear matter with
      // appropriate value of proton fraction, including the
      // possibility of drip particles if they exist.
      
      nsat=eos->fn0(1.0-2.0*protfrac,dtemp);
      np_left=nsat*protfrac;
      nn_left=nsat-np_left;
      neutron.n=nn_left;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      mun=neutron.mu;
      mup=proton.mu;
      nn_drip=0.0;
      np_drip=0.0;

      if (mun>neutron.m) {
	cout << "Neutron drip" << endl;

	sx[0]=neutron.n;
	sx[1]=proton.n;
	sx[2]=0.01;

	mm_funct ndf=std::bind
	  (std::mem_fn<int(size_t,const si_vector_t &,si_vector_t &)>
	   (&seminf_nr::ndripfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	
	nd.msolve(3,sx,ndf);
	ndripfun(3,sx,sy);
	if (fabs(sy[0])>1.0e-4 || fabs(sy[1])>1.0e-4 || fabs(sy[2])>1.0e-4) {
	  O2SCL_ERR("Neutron drip solution failed.",o2scl::exc_efailed);
	}

	nn_left=sx[0];
	np_left=sx[1];
	neutron.n=nn_left;
	proton.n=np_left;
	eos->calc_e(neutron,proton,hb);
	mun=neutron.mu;
	mup=proton.mu;
	nn_drip=sx[2];
	np_drip=0.0;

      }

      if (mup>proton.m) {
	cout << "Proton drip" << endl;

	sx[0]=neutron.n;
	sx[1]=proton.n;
	sx[2]=0.01;

	mm_funct pdf=std::bind
	  (std::mem_fn<int(size_t,const si_vector_t &,si_vector_t &)>
	   (&seminf_nr::pdripfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	nd.msolve(3,sx,pdf);
	pdripfun(3,sx,sy);
	if (fabs(sy[0])>1.0e-4 || fabs(sy[1])>1.0e-4 || fabs(sy[2])>1.0e-4) {
	  O2SCL_ERR("Proton drip solution failed.",o2scl::exc_efailed);
	}

	nn_left=sx[0];
	np_left=sx[1];
	eos->calc_e(neutron,proton,hb);
	mun=neutron.mu;
	mup=proton.mu;
	nn_drip=0.0;
	np_drip=sx[2];

      }

      nsat=nn_left+np_left;
      cout << "Saturation density at x=" <<  protfrac << ": " <<  nsat 
	   << endl;
      cout << "nn_left: " <<  nn_left << " np_left: " <<  np_left << endl;
      cout << "nn_drip: " << nn_drip << " np_drip: " << np_drip << endl;

      //--------------------------------------------
      // Calculate derivatives useful for the 
      // matching to the analytical solutions

      // dmundn
   
      neutron.n=nn_left;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      hin=neutron.mu;

      neutron.n=nn_left-1.0e-4;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      lon=neutron.mu;
    
      dmundn=1.0e4*(hin-lon);
  
      // dmundp
  
      neutron.n=nn_left;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      hin=neutron.mu;

      neutron.n=nn_left;
      proton.n=np_left-1.0e-4;
      eos->calc_e(neutron,proton,hb);
      lon=neutron.mu;
    
      dmundp=1.0e4*(hin-lon);
  
      //--------------------------------------------
      // Calculate derivatives of the RHS's numerically
  
      // At the left hand side
      neutron.n=nn_left;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      if (model=="APR") {
	eosa.gradient_qij2(neutron.n, proton.n, qnn, qnp, qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(neutron,proton,thth, qnn, qnp, qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=qnn*qpp-qnp*qnp;
      if (qnn_equals_qnp==true || (model=="APR" && fabs(det)<1.0e-5)) {
	rhsn=(neutron.mu-mun)/ qpp/2.0;
	rhsp=(proton.mu-mup)/ qnn/2.0;
      } else {
	rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	rhsn/=det;
	rhsp/=det;
      }
      hin=rhsn;
      hip=rhsp;

      neutron.n=nn_left-1.0e-4;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      if (model=="APR") {
	eosa.gradient_qij2(neutron.n, proton.n, qnn, qnp, qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(neutron,proton,thth, qnn, qnp, qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=qnn*qpp-qnp*qnp;
      if (qnn_equals_qnp==true || (model=="APR" && fabs(det)<1.0e-5)) {
	rhsn=(neutron.mu-mun)/ qpp/2.0;
	rhsp=(proton.mu-mup)/ qnn/2.0;
      } else {
	rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;

      dndnl=1.0e4*(hin-lon);
      dpdnl=1.0e4*(hip-lop);
  
      neutron.n=nn_left;
      proton.n=np_left-1.0e-4;
      eos->calc_e(neutron,proton,hb);
      if (model=="APR") {
	eosa.gradient_qij2(neutron.n, proton.n, qnn, qnp, qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(neutron,proton,thth, qnn, qnp, qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=qnn*qpp-qnp*qnp;
      if (qnn_equals_qnp==true || (model=="APR" && fabs(det)<1.0e-5)) {
	rhsn=(neutron.mu-mun)/ qpp/2.0;
	rhsp=(proton.mu-mup)/ qnn/2.0;
      } else {
	rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;
  
      dndpl=1.0e4*(hin-lon);
      dpdpl=1.0e4*(hip-lop);
  
      // At the right hand side
      neutron.n=1.0e-4;
      proton.n=1.0e-4;
      eos->calc_e(neutron,proton,hb);
      if (model=="APR") {
	eosa.gradient_qij2(neutron.n, proton.n, qnn, qnp, qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(neutron,proton,thth, qnn, qnp, qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=qnn*qpp-qnp*qnp;
      if (qnn_equals_qnp==true || (model=="APR" && fabs(det)<1.0e-5)) {
	rhsn=(neutron.mu-mun)/ qpp/2.0;
	rhsp=(proton.mu-mup)/ qnn/2.0;
      } else {
	rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	rhsn/=det;
	rhsp/=det;
      }
      hin=rhsn;
      hip=rhsp;
  
      neutron.n=1.0e-8;
      proton.n=1.0e-4;
      eos->calc_e(neutron,proton,hb);
      if (model=="APR") {
	eosa.gradient_qij2(neutron.n, proton.n, qnn, qnp, qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(neutron,proton,thth, qnn, qnp, qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=qnn*qpp-qnp*qnp;
      if (qnn_equals_qnp==true || (model=="APR" && fabs(det)<1.0e-5)) {
	rhsn=(neutron.mu-mun)/ qpp/2.0;
	rhsp=(proton.mu-mup)/ qnn/2.0;
      } else {
	rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;
  
      dndnr=1.0e4*(hin-lon);
      dpdnr=1.0e4*(hip-lop);
  
      neutron.n=1.0e-4;
      proton.n=1.0e-8;
      eos->calc_e(neutron,proton,hb);
      if (model=="APR") {
	eosa.gradient_qij2(neutron.n, proton.n, qnn, qnp, qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(neutron,proton,thth, qnn, qnp, qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=qnn*qpp-qnp*qnp;
      if (qnn_equals_qnp==true || (model=="APR" && fabs(det)<1.0e-5)) {
	rhsn=(neutron.mu-mun)/ qpp/2.0;
	rhsp=(proton.mu-mup)/ qnn/2.0;
      } else {
	rhsn=(neutron.mu-mun)*qpp-(proton.mu-mup)*qnp;
	rhsp=(proton.mu-mup)*qnn-(neutron.mu-mun)*qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;
  
      dndpr=1.0e4*(hin-lon);
      dpdpr=1.0e4*(hip-lop);
  
      goodbound=true;

      sqt=sqrt(dndnl*dndnl-2.0*dpdpl*dndnl+dpdpl*dpdpl+
	       4.0*dndpl*dpdnl);
      lex1=sqrt(2.0*(dndnl+dpdpl)+2.0*sqt)/2.0;
      lex2=sqrt(2.0*(dndnl+dpdpl)-2.0*sqt)/2.0;

      if (!std::isfinite(lex1)) {
	goodbound=false;
	lex1=0.0;
      }
      if (!std::isfinite(lex2)) {
	expo=lex1;
	lex2=0.0;
	goodbound=false;
      } else if (lex2< lex1) {
	expo=lex2;
      } else {
	expo=lex1;
      }

      sqt=sqrt(dndnr*dndnr-2.0*dpdpr*dndnr+dpdpr*dpdpr+
	       4.0*dndpr*dpdnr);
      rex1=-sqrt(2.0*(dndnr+dpdpr)+2.0*sqt)/2.0;
      rex2=-sqrt(2.0*(dndnr+dpdpr)-2.0*sqt)/2.0;
  
      if (!std::isfinite(rex1)) {
	goodbound=false;
	rex1=0.0;
      }
      if (!std::isfinite(rex2)) {
	goodbound=false;
	rex2=0.0;
      }

      //--------------------------------------------
      // Solve
      
      if (qnn_equals_qnp) {
	//solve_qnn_equal_qnp(monfact);
      } else {
	solve_qnn_neq_qnp(monfact);
      }

      //--------------------------------------------
      // Calculate the value of x for nn*0.9, nn*0.5, 
      // nn*0.1, etc.

      verbose=2;
      if (verbose>1) cout << "Lookups. " << endl;

      o2scl::interp<std::vector<double>,std::vector<double> > it(itp_linear);
      
      for(int i=0;i<3;i++) {
	xn[i]=it.eval(nn_drip+(nn_left-nn_drip)/10.0*((double)(i*4+1)),
		      tab_prof.get_nlines(),tab_prof.get_column("nn"),
		      tab_prof.get_column("x"));
	xp[i]=it.eval(np_drip+(np_left-np_drip)/10.0*((double)(i*4+1)),
		      tab_prof.get_nlines(),tab_prof.get_column("np"),
		      tab_prof.get_column("x"));
      }
      
      delta=1.0-2.0*protfrac;
      barn=nsat;
      hns=eos->fesym(barn);

      if (verbose>1) cout << "Integrands." << endl;
      for(int i=0;i<((int)tab_prof.get_nlines());i++) {
	// Fix negative values
	if (tab_prof[1][i]<0.0) tab_prof.set(1,i,0.0);
	if (tab_prof[2][i]<0.0) tab_prof.set(2,i,0.0);

	neutron.n=tab_prof[1][i];
	proton.n=tab_prof.get(2,i);

	eos->calc_e(neutron,proton,hb);

	if (model=="APR") {
	  eosa.gradient_qij2(neutron.n,proton.n,qnn,qnp,qpp,
			     dqnndnn,dqnndnp,dqnpdnn,
			     dqnpdnp,dqppdnn,dqppdnp);
	  tab_prof.set("qnn",i,qnn);
	  tab_prof.set("qnp",i,qnp);
	  tab_prof.set("qpp",i,qpp);
	}
	if (model=="gp") {
	  thermo thth;
	  eosg.gradient_qij(neutron,proton,thth,qnn,qnp,qpp,
			    dqnndnn,dqnndnp,dqnpdnn,
			    dqnpdnp,dqppdnn,dqppdnp);
	  tab_prof.set("qnn",i,qnn);
	  tab_prof.set("qnp",i,qnp);
	  tab_prof.set("qpp",i,qpp);
	} 

	tab_prof.set("rho",i,neutron.n+proton.n);
	tab_prof.set("alpha",i,neutron.n-proton.n);
	tab_prof.set("ebulk",i,hb.ed-mun*tab_prof.get(1,i)-
		     mup*tab_prof.get(2,i));
	tab_prof.set("egrad",i,0.5*qnn*tab_prof.get(3,i)*tab_prof.get(3,i)+
		     0.5*qpp*tab_prof.get(4,i)*tab_prof.get(4,i)+
		     qnp*tab_prof.get(3,i)*tab_prof.get(4,i));
	tab_prof.set("esurf",i,tab_prof.get("ebulk",i)+tab_prof.get("egrad",i));
	tab_prof.set("thickint",i,tab_prof.get(1,i)/nn_left-
		     tab_prof.get(2,i)/np_left);
	if (fabs(protfrac-0.5)>1.0e-8) {
	  tab_prof.set("wdint",i,tab_prof.get("alpha",i)/delta-
		       tab_prof.get("rho",i));
	  if (tab_prof.get("rho",i)>=1.0e-5) {
	    barn=tab_prof.get("rho",i);
	    tab_prof.set("wd2int",i,tab_prof.get("rho",i)*
			 (pow(tab_prof.get("alpha",i)/delta/
			      tab_prof.get("rho",i),2.0)*
			  eos->fesym(barn)-hns));
	  }
	}
      }
    
      if (verbose>1) cout << "Integrals." << endl;

      interp<std::vector<double>,std::vector<double> > gi;
      surf=gi.integ(tab_prof.get(0,0),tab_prof[0][tab_prof.get_nlines()-1],
		    tab_prof.get_nlines(),
		    tab_prof[0],(tab_prof.get_column("esurf")));
      sbulk=gi.integ(tab_prof.get(0,0),tab_prof[0][tab_prof.get_nlines()-1],
		     tab_prof.get_nlines(),
		     tab_prof[0],(tab_prof.get_column("ebulk")));
      sgrad=gi.integ(tab_prof.get(0,0),tab_prof[0][tab_prof.get_nlines()-1],
		     tab_prof.get_nlines(),
		     tab_prof[0],(tab_prof.get_column("egrad")));
      thick=gi.integ(tab_prof.get(0,0),tab_prof[0][tab_prof.get_nlines()-1],
		     tab_prof.get_nlines(),
		     tab_prof[0],(tab_prof.get_column("thickint")));
      if (verbose>1) {
	cout << "surf: " << surf << endl;
	cout << "thick: " << thick << endl;
	cout << "sbulk: " << sbulk << endl;
	cout << "sgrad: " << sgrad << endl;
      }
      if (fabs(protfrac-0.5)>1.0e-8) {
	wd=gi.integ(tab_prof.get(0,0),tab_prof[0][tab_prof.get_nlines()-1],
		    tab_prof.get_nlines(),
		    tab_prof[0],(tab_prof.get_column("wdint")));
	wd*=hns;
      
	wd2=gi.integ(tab_prof.get(0,0),tab_prof[0][tab_prof.get_nlines()-1],
		     tab_prof.get_nlines(),
		     tab_prof[0],(tab_prof.get_column("wd2int")));
      
	sssv_drop=4*pi*r0*r0*wd/hns;
      
	w0jl=0.0;
	wdjl=0.0;
	sssv_jim=0.0;
	if (verbose>1) {
	  cout << "wd: " << wd << endl;
	  cout << "wd2: " << wd2 << endl;
	  cout << "sssv_drop: " << sssv_drop << endl;
	}
      
      }

      double xline[8]={protfrac,surf,thick,sbulk,sgrad,wd,wd2,
		       nn_drip};
      tab_summ.line_of_data(8,xline);
    
      //--------------------------------------------
      // Add columns to check table:
    
      if (true && qnn_equals_qnp==false && 
	  model!="schematic" && model!="APR" && model!="gp") {
	double tx;
	si_vector_t ty(5), tdy(5);
	if (!tab_prof.is_column("nnpp")) tab_prof.new_column("nnpp");
	if (!tab_prof.is_column("nppp")) tab_prof.new_column("nppp");
	if (!tab_prof.is_column("rhsn")) tab_prof.new_column("rhsn");
	if (!tab_prof.is_column("rhsp")) tab_prof.new_column("rhsp");
	if (!tab_prof.is_column("lhs_n")) tab_prof.new_column("lhs_n");
	if (!tab_prof.is_column("lhs_a")) tab_prof.new_column("lhs_a");
	if (!tab_prof.is_column("rhs_n")) tab_prof.new_column("rhs_n");
	if (!tab_prof.is_column("rhs_a")) tab_prof.new_column("rhs_a");
	tab_prof.set("nnpp",0,0.0);
	tab_prof.set("nppp",0,0.0);
	tab_prof.set("rhsn",0,0.0);
	tab_prof.set("rhsp",0,0.0);
	tab_prof.set("lhs_n",0,0.0);
	tab_prof.set("lhs_a",0,0.0);
	tab_prof.set("rhs_n",0,0.0);
	tab_prof.set("rhs_a",0,0.0);

	low_dens_part=false;

	for(int i=1;i<((int)tab_prof.get_nlines());i++) {
	  
	  tab_prof.set("nnpp",i,(tab_prof.get("nnp",i)-
				 tab_prof.get("nnp",i-1))/
		       (tab_prof.get("x",i)-tab_prof.get("x",i-1)));
	  tab_prof.set("nppp",i,(tab_prof.get("npp",i)-
				 tab_prof.get("npp",i-1))/
		       (tab_prof.get("x",i)-tab_prof.get("x",i-1)));
	  
	  tx=(tab_prof.get("x",i)+tab_prof.get("x",i-1))/2.0;
	  ty[0]=(tab_prof.get("nn",i)+tab_prof.get("nn",i-1))/2.0;
	  ty[1]=(tab_prof.get("np",i)+tab_prof.get("np",i-1))/2.0;
	  ty[2]=(tab_prof.get("nnp",i)+tab_prof.get("nnp",i-1))/2.0;
	  ty[3]=(tab_prof.get("npp",i)+tab_prof.get("npp",i-1))/2.0;

	  if (ty[1]<=0.0) low_dens_part=true;
	  derivs(tx,ty,tdy);
	  
	  tab_prof.set("rhsn",i,tdy[2]);
	  tab_prof.set("rhsp",i,tdy[3]);
	  
	  if (ty[1]>0.0) {
	    tab_prof.set("lhs_n",i,(neutron.mu-mun+proton.mu-mup)/2.0);
	    tab_prof.set("lhs_a",i,(neutron.mu-mun-proton.mu+mup)/2.0);
	    tab_prof.set("rhs_n",i,(tab_prof.get("nnpp",i)+
				    tab_prof.get("nppp",i))/2.0*
			 (qnn+qnp));
	    tab_prof.set("rhs_a",i,(tab_prof.get("nnpp",i)-
				    tab_prof.get("nppp",i))/2.0*
			 (qnn-qnp));
	  } else {
	    tab_prof.set("lhs_n",i,neutron.mu-mun);
	    tab_prof.set("rhs_n",i,tab_prof.get("nnpp",i)*(qnn));
	    tab_prof.set("lhs_a",i,0.0);
	    tab_prof.set("rhs_a",i,0.0);
	  }
	}
      }

      string tablename=((string)"nr")+std::to_string(pf_index);
      hf.open_or_create(out_file);
      hdf_output(hf,tab_prof,tablename);
      hf.close();
      cout << "Wrote solution to file " << out_file << " ." << endl;
      cout << endl;

      // Loop for next proton fraction
    }
    
    hf.open_or_create(out_file);
    hdf_output(hf,tab_summ,"summary");
    hf.close();
    cout << "Wrote summary to file " << out_file << " ." << endl;
    cout << endl;
    
    return 0;
  }

  /** \brief Compute standard results and compare with stored output
   */
  int check(std::vector<std::string> &sv, bool itive_com) {

    vector<string> sv2={"calc","0.45","0.46"};
    calc(sv2,itive_com);
    
    test_mgr t;
    t.set_output_level(2);
    
    hdf_file hf;
    string name;
    table_units<> tab, tab_expected;

    name="nr1";

    hf.open("nr.o2");
    hdf_input(hf,tab,name);
    hf.close();

    hf.open("nr_save.o2");
    hdf_input(hf,tab_expected,name);
    hf.close();

    t.test_rel_nonzero_table(tab,tab_expected,1.0e-8,1.0e-8,"table 1");

    name="nr2";

    hf.open("nr.o2");
    hdf_input(hf,tab,name);
    hf.close();

    hf.open("nr_save.o2");
    hdf_input(hf,tab_expected,name);
    hf.close();

    t.test_rel_nonzero_table(tab,tab_expected,1.0e-12,1.0e-8,"table 2");
    
    if (!t.report()) {
      exit(-1);
    }
    
    return 0;
  }
  
};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  seminf_nr sn;
  cli cl;
  
  static const int nopt=2;
  o2scl::comm_option_s options[nopt]={
    {0,"calc","",-1,-1,"<proton fraction 1> [pf2] [pf3] ...","",
     new o2scl::comm_option_mfptr<seminf_nr>
     (&sn,&seminf_nr::calc),o2scl::cli::comm_option_both},
    {0,"check","",0,0,"","",
     new o2scl::comm_option_mfptr<seminf_nr>
     (&sn,&seminf_nr::check),o2scl::cli::comm_option_both}
  };

  cl.set_comm_option_vec(nopt,options);
  cl.gnu_intro=false;
  
  o2scl::cli::parameter_double p_rhs_adjust;
  p_rhs_adjust.d=&sn.rhs_adjust;
  p_rhs_adjust.help="RHS adjustment (default 0.0)";
  cl.par_list.insert(make_pair("rhs_adjust",&p_rhs_adjust));

  o2scl::cli::parameter_string p_model;
  p_model.str=&sn.model;
  p_model.help="Model (default \"NRAPR\")";
  cl.par_list.insert(make_pair("model",&p_model));

  o2scl::cli::parameter_string p_out_file;
  p_out_file.str=&sn.out_file;
  p_out_file.help="Out_File (default \"nr.o2\")";
  cl.par_list.insert(make_pair("out_file",&p_out_file));

  cl.run_auto(argc,argv);
  
  return 0;
}
