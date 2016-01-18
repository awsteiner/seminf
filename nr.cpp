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
#include <iostream>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/eos_had_apr.h>
#include <o2scl/table.h>
#include <o2scl/constants.h>
#include <o2scl/part.h>
#include <o2scl/eos_had_schematic.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/eos_had_potential.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/lib_settings.h>
#include <o2scl/ode_it_solve.h>

#include "seminf.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

//--------------------------------------------
// Global variables

static const double big=3.0;
static const int ngrid=100;

/** \brief Semi-infinite nuclear matter for nonrelativistic
    models in the Thomas-Fermi approximation
    
    \warning Does not work with dripped nucleons yet.
    
    \note The first proton fraction cannot be 0.5, as the 
    algorithm cannot handle later values different from 0.5
    
    12/15/03 - Implemented a new solution method for Qnn!=Qnp for when
    nndrip becomes > 0. It seems that in this case, it is better not
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
  bool showrelax;
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
  bool relaxfile;
  /// If true, \f$ Q_{nn}=Q_{np} \f$
  bool qnn_equals_qnp;
  /// If true, match to exponential decay on the RHS
  bool rhsmode;
  /// Minimum err
  double minerr;
  /// Desc
  double monfact;
  /// Desc
  double expo;
  /// Neutron drip density
  double nndrip;
  /// Proton drip density
  double npdrip;
  /// Desc (default false)
  bool flatden;
  /// Nonlinear equation solver
  mroot_hybrids<> nd;
  /// Desc
  double nnrhs;
  /// Desc
  double nprhs;
  /// Desc
  double rhslength;

  /// If true, we have good boundary conditions
  bool goodbound;
  /// Initial derivative for densities on LHS (typically negative)
  double firstderiv;

  /// \name Desc
  //@{
  ubvector xstor;
  ubmatrix ystor;
  //@}
  
  /// \name Desc
  //@{
  double xg1[ngrid+1];
  double xg2[ngrid+1];
  double yg1[6][ngrid+1];
  double yg2[6][ngrid+1];
  //@}

  /** \brief Desc
   */
  double relaxconverge;
  /** \brief Desc
   */
  double rhsmin;
  /** \brief The step size factor for constructing the initial 
      guess (default 7.0)
  */
  double initialstep;
  /// Type of EOS model in use
  string model;
  /// Integrator to compute integrals of final solution
  inte_qagiu_gsl<> gl2;
  /// Store the solution
  table<> at;
  /// \name Nucleon-nucleon interactions
  //@{
  eos_had_apr eosa;
  eos_had_potential eosg;
  eos_had_skyrme eoss;
  eos_had_schematic eosp;
  //@}
  /// To compute derivatives
  deriv_gsl<> df;
  
  /** \brief Desc
   */
  int derivs(double sx, const ubvector &sy, ubvector &dydx) {
    double rhsn, rhsp=0.0, det;
    double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

    neutron.n=sy[0];
    proton.n=sy[1];
    if (model=="apr") {
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
	if (model!="apr" && model!="gp") {
	  dydx[3]=(proton.mu-mup)/ qpp;
	} else {
	  dydx[3]=(proton.mu-mup+0.5*dqppdnp*sy[3]*sy[3])/ qpp;
	}
      }
    } else if (proton.n<=0.0) {
      dydx[0]=sy[2];
      dydx[1]=0.0;
      if (model!="apr" && model!="gp") {
	dydx[2]=(neutron.mu-mun)/ qnn;
      } else {
	dydx[2]=(neutron.mu-mun+0.5*dqnndnn*sy[2]*sy[2])/ qnn;
      }
      dydx[3]=0.0;
    } else {
      if (rhsmode==false) {
	if (model!="apr" && model!="gp") {
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
	if (model!="apr" && model!="gp") {
	  rhsn=(neutron.mu-mun)/ qnn;
	  rhsp=0.0;
	} else {
	  rhsn=(neutron.mu-mun+0.5*dqnndnn*sy[2]*sy[2])/ qnn;
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
  double solve_qnn_neq_qnp(double lmonfact, int argc) {
    bool guessdone;
    int debug=0;
    double dx, xrhs, delta, epsi;
    ubvector y(5), dydx(5);
    int ilast=0, interpi;
    int n_eq=5, n_b;
    ubvector ox(ngrid);
    ubmatrix oy(ngrid,n_eq);

    monfact=lmonfact;

    //----------------------------------------------
    // Create object

    if (flatden) {
      n_b=3;
    } else {
      n_b=4;
    }
    
    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&seminf_nr::difeq),this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&seminf_nr::left),this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&seminf_nr::right),this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);   
    ode_it_solve2 oit;
    ubmatrix A(ngrid*n_eq,ngrid*n_eq);
    ubvector rhs(ngrid*n_eq),dy(ngrid*n_eq);
    
    if (pf_index==1) {

      //----------------------------------------------
      // Construct a guess by shooting

      xstor[1]=0.0;
      ystor(1,1)=nn_left;
      ystor(2,1)=np_left;
      ystor(3,1)=firstderiv;
      ystor(4,1)=firstderiv;
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
      
	if (ystor(1,i)<nndrip || ystor(2,i)<npdrip) {
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
      // nndrip or npdrip. In that case, just use the entire guess.
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
	  oy(i,j)=yg1[j+1][i+1];
	}
      }

    }

    //----------------------------------------------
    // Try ode_it_solve

    rhsmode=false;

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

    oit.tol_rel=relaxconverge;
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
	yg1[j][i]=oy(i-1,j-1);
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
    if (npdrip==0.0) {
      nnrhs=oy(ngrid-1,0);
      cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
      if (nnrhs<rhsmin) {
	nnrhs=rhsmin;
	cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
      }
    } else {
      nprhs=oy(ngrid-1,1);
      cout << "RHS, nprhs=" << nprhs << endl;
      if (nprhs<rhsmin) {
	nprhs=rhsmin;
	cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
      }
    }

    //----------------------------------------------
    // Store solution and delete seminf_nr_relax

    at.set_nlines(ngrid);
    for(int i=0;i<ngrid;i++) {
      at.set(0,i,ox[i]);
      for(int j=0;j<5;j++) {
	at.set(j+1,i,oy(i,j));
      }
    }

    //----------------------------------------------
    // RHS
    
    if (fabs(protfrac-0.5)>1.0e-4) {

      at.set_nlines(ngrid*2-1);

      if (pf_index==1) {

	//----------------------------------------------
	// Construct a simple linear guess
      
	if (npdrip==0.0) {
	  for(int i=0;i<ngrid;i++) {
	    ox[i]=((double)i)/((double)(ngrid-1));
	    oy(i,0)=(nnrhs-nndrip)*(1.0-ox[i])+nndrip;
	    oy(i,1)=-0.08*(1.0-ox[i]);
	    oy(i,2)=0.001;
	  }
	} else {
	  for(int i=0;i<ngrid;i++) {
	    ox[i]=((double)i)/((double)(ngrid-1));
	    oy(i,0)=(nprhs-npdrip)*(1.0-ox[i])+npdrip;
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
	    oy(i,j)=yg2[j+1][i+1];
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
      
      rhsmode=true;
      cout << "second solve: " << nndrip << endl;
      
      if (nndrip>0.0) {
	O2SCL_ERR("Neutron drip not yet supported.",exc_eunimpl);
      } else {
	oit.tol_rel=relaxconverge;
	oit.verbose=1;
	oit.solve(ngrid,n_eq,n_b,ox,oy,f_derivs,f_left,f_right,
		  A,rhs,dy);
      }

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
	  yg2[j][i]=oy(i-1,j-1);
	}
      }
  
      //----------------------------------------------
      // Store solution and delete seminf_nr_relax
      
      for(int i=0;i<ngrid;i++) {
	at.set(0,i+ngrid-1,ox(i));
	for(int j=0;j<3;j++) {
	  at.set(j+1,i+ngrid-1,oy(i,j));
	}
      }
      //delete rel;

      //----------------------------------------------
      // Rearrangment and rescaling

      if (npdrip==0.0) {
	for(int i=ngrid-1;i<2*ngrid-1;i++) {
	  at.set(5,i,at.get(3,i));
	  at.set(4,i,0.0);
	  at.set(3,i,at.get(2,i));
	  at.set(2,i,0.0);
	  at.set(0,i,at.get(0,i)*at.get(5,i)+xrhs);
	}
      } else {
	for(int i=1;i<=ngrid;i++) {
	  at.set(5,i,at.get(3,i));
	  at.set(4,i,at.get(2,i));
	  at.set(3,i,0.0);
	  at.set(2,i,at.get(1,i));
	  at.set(1,i,0.0);
	  at.set(0,i,at.get(0,i)*at.get(5,i)+xrhs);
	}
      }
      //    at.set_nlines(2*ngrid-1);
    }

    if (debug>0) {
      for(int i=0;i<2*ngrid-1;i+=9) {
	if (i==108) i-=8;
	cout.precision(4);
	cout.width(3);
	cout << i << " ";
	for(int j=0;j<6;j++) {
	  cout << at.get(j,i) << " ";
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
    ubvector sx(3), sy(3);
    int n_eq=3, n_b=2;

    monfact=lmonfact;

    //----------------------------------------------
    // Create object

    rel=new seminf_nr_relax(3,2,ngrid);
    rel->itmax=200;
    rel->rp=&rp;

    //----------------------------------------------
    // Construct a guess by shooting

    if (pf_index==1) {
      xstor[1]=0.0;
      ystor(1,1)=nn_left+np_left;
      ystor(2,1)=firstderiv;
      dx=initialstep/((double)ngrid);
    
      guessdone=false;
      for(int i=2;guessdone==false && i<=ngrid;i++) {
	xstor[i]=xstor[i-1]+dx;
      
	sx[1]=(1.0-protfrac)*ystor(1,i-1);
	sx[2]=protfrac*ystor(1,i-1);
	barn=ystor(1,i-1);
	mm_funct11 qqf=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	   (&seminf_nr::qnnqnpfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	nd.msolve(2,sx,qqf);
      
	ystor(1,i)=ystor(1,i-1)+dx*ystor(2,i-1);
	ystor(2,i)=ystor(2,i-1)+dx*(neutron.mu-mun)/(qnn);
      
	if (sx[2]<nndrip+npdrip+1.0e-6 || ystor(2,i)>0.0) {
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
	  rel->y(j,i)=yg1[j][i];
	}
      }
    }

    //----------------------------------------------
    // Solve

    rhsmode=false;
    int relret1=rel->solve(relaxconverge,1.0);
  
    //----------------------------------------------
    // Copy solution for next guess

    for(int i=1;i<=ngrid;i++) {
      xg1[i]=rel->x[i];
      for(int j=1;j<=n_eq;j++) {
	yg1[j][i]=rel->y(j,i);
      }
    }
  
    //----------------------------------------------
    // Rescale solution and set neutron and proton
    // densities:

    at.set_nlines(ngrid);
    for(int i=1;i<=rel->ngrid;i++) {

      sx[1]=(1.0-protfrac)*rel->y(1,i);
      sx[2]=protfrac*rel->y(1,i);

      barn=rel->y(1,i);
      mm_funct11 qqf=std::bind
	(std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	 (&seminf_nr::qnnqnpfun),this,std::placeholders::_1,
	 std::placeholders::_2,std::placeholders::_3);
      nd.msolve(2,sx,qqf);

      at.set(0,i-1,rel->x[i]*rel->y(3,i));
      at.set(1,i-1,sx[1]);
      at.set(2,i-1,sx[2]);
      at.set(3,i-1,rel->y(2,i));
      at.set(4,i-1,0.0);
      at.set(5,i-1,rel->y(3,i));
    }

    //----------------------------------------------
    // Store quantites at rhs for later use:
  
    xrhs=at[0][ngrid-1];
    if (npdrip==0.0) {
      nnrhs=at[1][ngrid-1];
      cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
      if (nnrhs<rhsmin) {
	nnrhs=rhsmin;
	cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
      }
    } else {
      nprhs=at[2][ngrid-1];
      cout << "RHS, nprhs=" << nprhs << endl;
      if (nprhs<rhsmin) {
	nprhs=rhsmin;
	cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
      }
    }

    //----------------------------------------------
    // Delete seminf_nr_relax

    delete rel;

    //---------------------------------------------
    // RHS
    
    if (fabs(protfrac-0.5)>1.0e-4) {
      rel=new seminf_nr_relax(3,1,ngrid);
      rel->rp=&rp;
    
      if (pf_index==1) {
	//----------------------------------------------
	// Construct a simple linear guess
      
	if (npdrip==0.0) {
	  for(int i=1;i<=ngrid;i++) {
	    rel->x[i]=((double)(i-1))/((double)(ngrid-1));
	    rel->y(1,i)=(nnrhs-nndrip)*(1.0-rel->x[i])+nndrip;
	    rel->y(2,i)=-0.08*(1.0-rel->x[i]);
	    rel->y(3,i)=0.001;
	  }
	} else {
	  for(int i=1;i<=ngrid;i++) {
	    rel->x[i]=((double)(i-1))/((double)(ngrid-1));
	    rel->y(1,i)=(nprhs-npdrip)*(1.0-rel->x[i])+npdrip;
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
	    rel->y(j,i)=yg2[j][i];
	  }
	}
      
      }
    
      //----------------------------------------------
      // Solve
    
      rhsmode=true;
      int relret2=rel->solve(relaxconverge,1.0);
    
      //----------------------------------------------
      // Copy solution for next guess
      
      for(int i=1;i<=ngrid;i++) {
	xg2[i]=rel->x[i];
	for(int j=1;j<=n_eq;j++) {
	  yg2[j][i]=rel->y(j,i);
	}
      }
    
      //----------------------------------------------
      // Store solution, rescale, and delete seminf_nr_relax
    
      at.set_nlines(2*ngrid-1);
      for(int i=1;i<=rel->ngrid;i++) {
	at.set(5,i-2+ngrid,rel->y(3,i));
	at.set(3,i-2+ngrid,rel->y(2,i));
	at.set(4,i-2+ngrid,0.0);
	at.set(2,i-2+ngrid,0.0);
	at.set(1,i-2+ngrid,rel->y(1,i));
	at.set(0,i-2+ngrid,rel->x[i]*rel->y(3,i)+xrhs);
      }
      delete rel;
    }
  
    return 0.0;

  }
#endif
  
  /** \brief Desc
   */
  int qnnqnpfun(size_t sn, const ubvector &sx, ubvector &sy) {

    neutron.n=sx[1];
    proton.n=sx[2];
    if (sx[1]<0.0 || sx[2]<0.0) return 1;
  
    eos->calc_e(neutron,proton,hb);

    sy[1]=neutron.n+proton.n-barn;
    sy[2]=neutron.mu-proton.mu-mun+mup;

    return 0;
  }

  /** \brief Desc
   */
  int ndripfun(size_t sn, const ubvector &sx, ubvector &sy) {
    double pleft, pright, munleft, munright;

    if (sx[1]<0.0 || sx[2]<0.0 || sx[3]<0.0) return 1;

    neutron.n=sx[1];
    proton.n=sx[2];
    sy[1]=proton.n-protfrac*(proton.n+neutron.n);
  
    eos->calc_e(neutron,proton,hb);
    pleft=hb.pr;
    munleft=neutron.mu;

    neutron.n=sx[3];
    proton.n=0.0;
  
    eos->calc_e(neutron,proton,hb);
    pright=hb.pr;
    munright=neutron.mu;

    sy[2]=pleft-pright;
    sy[3]=munleft-munright;

    return 0;
  }

  /** \brief Desc
   */
  int pdripfun(size_t sn, const ubvector &sx, ubvector &sy) {
    double pleft, pright, mupleft, mupright;

    if (sx[1]<0.0 || sx[2]<0.0 || sx[3]<0.0) return 1;

    neutron.n=sx[1];
    proton.n=sx[2];
    sy[1]=proton.n-protfrac*(proton.n+neutron.n);
  
    eos->calc_e(neutron,proton,hb);
    pleft=hb.pr;
    mupleft=proton.mu;

    neutron.n=0.0;
    proton.n=sx[3];
  
    eos->calc_e(neutron,proton,hb);
    pright=hb.pr;
    mupright=proton.mu;

    sy[2]=pleft-pright;
    sy[3]=mupleft-mupright;

    return 0;
  }

  /** \brief Future function for \ref o2scl::ode_it_solve
   */
  double difeq(size_t ieq, double x, ubmatrix_row &y) {
    if (rhsmode==false) {
      ubvector y2=y, dydx(5);
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
      ubvector y2(5), dydx(5);
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

  /** \brief Future function for \ref o2scl::ode_it_solve
   */
  double left(size_t ieq, double x, ubmatrix_row &y) {

    if (rhsmode==false) {
      
      double delta=nn_left*monfact/exp(big*expo);
      double epsi=delta*(-dmundn+qnn*expo*expo)/
	(dmundp-qnp*expo*expo);
      
      if (epsi*delta<0.0) {
	cout << "epsi and delta have different signs" << endl;
	exit(-1);
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

    // If rhsmode is true
    return y[0]-nnrhs;
  }

  /** \brief Future function for \ref o2scl::ode_it_solve
   */
  double right(size_t ieq, double x, ubmatrix_row &y) {
    if (rhsmode==false) {
      return y[1];
    }
    // If rhsmode is true
    if (ieq==0) {
      return y[0]-nndrip;
    }
    return y[1];
  }
  
public:
  
  seminf_nr() {
    rhslength=2.0;
  }

  //--------------------------------------------
  // Function prototypes
  
  int run(int argc, char *argv[]) {
    
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
    double surf, sbulk, sgrad, sssv_drop=0.0, *pflist;
    ubvector sx(5), sy(5);
    double sssv_jim=0.0, hns, delta, w0jl=0.0, wdjl=0.0, den, w0;
    // Locations of fixed relative densities
    double xn[3], xp[3];
    double thick;
    int npf;
    ofstream fout;
    double lon, lop, hin, hip, sqt;
    double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

    int verbose=1;
    
    // Not used ATM but probably useful later
    //ubvector rho, alpha, ebulk, egrad, esurf, thickint;
    //ubvector wdint, wd2int;
    //ubvector vqnn, vqnp, vqpp;
    
    nd.tol_abs=1.0e-11;
    nd.tol_rel=1.0e-8;
    nd.ntrial=100;
    nd.verbose=0;

    initialstep=7.0;

    //--------------------------------------------
    // Get model parameters

    model="skyrme";

    if (model==((string)"skyrme")) {

      eos=&eoss;
      skyrme_load(eoss,"NRAPR");
      
      // Determine values of Q_{nn}, Q_{pp}, and Q_{np}
      qnn=0.1875*(eoss.t1*(1.0-eoss.x1)-eoss.t2*(1.0+eoss.x2));
      qpp=qnn;
      qnp=0.125*(3.0*eoss.t1*(1.0+eoss.x1/2.0)-
		 eoss.t2*(1.0+eoss.x2/2.0));
      
    } else if (model==((string)"schematic")) {
      eos=&eosp;
    } else if (model==((string)"apr")) {
      eos=&eosa;
    } else if (model==((string)"gp")) {
      eos=&eosg;
    } else {
      cout << "Unknown EOS" << endl;
      exit(-1);
    }
  
    //--------------------------------------------
    // Test if Qnn=Qnp

    if (model!="gp" && model!="apr" && fabs(qnn-qnp)<1.0e-7) {
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
  
    flatden=false;

    //--------------------------------------------

    at.set_nlines(ngrid*2);
    at.line_of_names(((string)"x nn np nnp npp scale rho alpha ")+
		     "ebulk egrad esurf thickint wdint wd2int");

    xstor.resize(ngrid+1);
    ystor.resize(5+1,ngrid+1);

    minerr=1.0e-6;
    relaxconverge=1.0e-9;
    showrelax=true;
    relaxfile=false;
    firstderiv=-2.0e-4;
    rhsmin=1.0e-7;
    monfact=0.002;

    neutron.init(o2scl_settings.get_convert_units().convert
		 ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    proton.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_proton),2.0);
    neutron.non_interacting=false;
    proton.non_interacting=false;
    eos->set_n_and_p(neutron,proton);
    eos->set_thermo(hb);

    // Compute saturation density
    n0half=eos->fn0(0.0,dtemp);
    double r0=cbrt(0.75/pi/ n0half);

    pflist=new double[2];
    pflist[0]=0.45;
    pflist[1]=0.46;
    npf=2;
  
    for(pf_index=1;pf_index<=npf;pf_index++) {
      protfrac=pflist[pf_index-1];

      //------------------------------------------------------
      // Calculate properties of saturated nuclear matter with
      // appropriate value of proton fraction including the
      // possibility of drip particles if they exist
      
      nsat=eos->fn0(1.0-2.0*protfrac,dtemp);
      np_left=nsat*protfrac;
      nn_left=nsat-np_left;
      neutron.n=nn_left;
      proton.n=np_left;
      eos->calc_e(neutron,proton,hb);
      mun=neutron.mu;
      mup=proton.mu;
      nndrip=0.0;
      npdrip=0.0;

      if (mun> neutron.m) {
	cout << "Neutron drip" << endl;

	sx[0]=neutron.n;
	sx[1]=proton.n;
	sx[2]=0.01;
	mm_funct11 ndf=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	   (&seminf_nr::ndripfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	nd.msolve(3,sx,ndf);
	ndripfun(3,sx,sy);
	if (fabs(sy[0])>1.0e-4 || fabs(sy[1])>1.0e-4 || fabs(sy[2])>1.0e-4) {
	  cout << "Error in solution of ndripfun." << endl;
	  exit(-1);
	}

	nn_left=sx[0];
	np_left=sx[1];
	neutron.n=nn_left;
	proton.n=np_left;
	eos->calc_e(neutron,proton,hb);
	mun=neutron.mu;
	mup=proton.mu;
	nndrip=sx[2];
	npdrip=0.0;

      }

      if (mup> proton.m) {
	cout << "Proton drip" << endl;

	sx[0]=neutron.n;
	sx[1]=proton.n;
	sx[2]=0.01;
	mm_funct11 pdf=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	   (&seminf_nr::pdripfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	nd.msolve(3,sx,pdf);
	pdripfun(3,sx,sy);
	if (fabs(sy[0])>1.0e-4 || fabs(sy[1])>1.0e-4 || fabs(sy[2])>1.0e-4) {
	  cout << "Error in solution of pdripfun." << endl;
	  exit(-1);
	}

	nn_left=sx[0];
	np_left=sx[1];
	eos->calc_e(neutron,proton,hb);
	mun=neutron.mu;
	mup=proton.mu;
	nndrip=0.0;
	npdrip=sx[2];

      }

      nsat=nn_left+np_left;
      cout << "Saturation density at x=" <<  protfrac << ": " <<  nsat 
	   << endl;
      cout << "nn_left: " <<  nn_left << " np_left: " <<  np_left << endl;
      cout << "nndrip: " << nndrip << " npdrip: " << npdrip << endl;

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
      if (model=="apr") {
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
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
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
      if (model=="apr") {
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
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
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
      if (model=="apr") {
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
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
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
      if (model=="apr") {
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
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
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
      if (model=="apr") {
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
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
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
      if (model=="apr") {
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
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
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
	solve_qnn_neq_qnp(monfact,argc);
      }

      //--------------------------------------------
      // Calculate the value of x for nn*0.9, nn*0.5, 
      // nn*0.1, etc.

      verbose=2;
      if (verbose>1) cout << "Lookups. " << endl;
      cout << at.get_nlines() << endl;

      o2scl::interp<std::vector<double>,std::vector<double> > it(itp_linear);
      
      for(int i=0;i<3;i++) {
	xn[i]=it.eval(nndrip+(nn_left-nndrip)/10.0*((double)(i*4+1)),
		      at.get_nlines(),at.get_column("nn"),
		      at.get_column("x"));
	xp[i]=it.eval(npdrip+(np_left-npdrip)/10.0*((double)(i*4+1)),
		      at.get_nlines(),at.get_column("np"),
		      at.get_column("x"));
      }
      
      delta=1.0-2.0*protfrac;
      barn=nsat;
      hns=eos->fesym(barn);

      if (verbose>1) cout << "Integrands." << endl;
      for(int i=0;i<((int)at.get_nlines());i++) {
	// Fix negative values
	if (at[1][i]<0.0) at.set(1,i,0.0);
	if (at[2][i]<0.0) at.set(2,i,0.0);

	neutron.n=at[1][i];
	proton.n=at.get(2,i);

	eos->calc_e(neutron,proton,hb);

	if (model=="apr") {
	  eosa.gradient_qij2(neutron.n,proton.n,qnn,qnp,qpp,
			      dqnndnn,dqnndnp,dqnpdnn,
			      dqnpdnp,dqppdnn,dqppdnp);
	  at.set("qnn",i,qnn);
	  at.set("qnp",i,qnp);
	  at.set("qpp",i,qpp);
	}
	if (model=="gp") {
	  thermo thth;
	  eosg.gradient_qij(neutron,proton,thth,qnn,qnp,qpp,
			    dqnndnn,dqnndnp,dqnpdnn,
			    dqnpdnp,dqppdnn,dqppdnp);
	  at.set("qnn",i,qnn);
	  at.set("qnp",i,qnp);
	  at.set("qpp",i,qpp);
	} 

	at.set("rho",i,neutron.n+proton.n);
	at.set("alpha",i,neutron.n-proton.n);
	at.set("ebulk",i,hb.ed-mun*at.get(1,i)-mup*at.get(2,i));
	at.set("egrad",i,0.5*qnn*at.get(3,i)*at.get(3,i)+
	       0.5*qpp*at.get(4,i)*at.get(4,i)+
	       qnp*at.get(3,i)*at.get(4,i));
	at.set("esurf",i,at.get("ebulk",i)+at.get("egrad",i));
	at.set("thickint",i,at.get(1,i)/ nn_left-at.get(2,i)/np_left);
	if (fabs(protfrac-0.5)>1.0e-8) {
	  at.set("wdint",i,at.get("alpha",i)/delta-at.get("rho",i));
	  if (at.get("rho",i)>=1.0e-5) {
	    barn=at.get("rho",i);
	    at.set("wd2int",i,at.get("rho",i)*
		   (pow(at.get("alpha",i)/delta/at.get("rho",i),2.0)*
		    eos->fesym(barn)-hns));
	  }
	}
      }
    
      if (verbose>1) cout << "Integrals." << endl;

      interp<std::vector<double>,std::vector<double> > gi;
      surf=gi.integ(at.get(0,0),at[0][at.get_nlines()-1],at.get_nlines(),
		    at[0],(at.get_column("esurf")));
      sbulk=gi.integ(at.get(0,0),at[0][at.get_nlines()-1],at.get_nlines(),
		     at[0],(at.get_column("ebulk")));
      sgrad=gi.integ(at.get(0,0),at[0][at.get_nlines()-1],at.get_nlines(),
		     at[0],(at.get_column("egrad")));
      thick=gi.integ(at.get(0,0),at[0][at.get_nlines()-1],at.get_nlines(),
		     at[0],(at.get_column("thickint")));
      if (verbose>1) cout << "thick: " << thick << endl;
      if (fabs(protfrac-0.5)>1.0e-8) {
	wd=gi.integ(at.get(0,0),at[0][at.get_nlines()-1],at.get_nlines(),
		    at[0],(at.get_column("wdint")));
	wd*=hns;
      
	wd2=gi.integ(at.get(0,0),at[0][at.get_nlines()-1],at.get_nlines(),
		     at[0],(at.get_column("wd2int")));
      
	sssv_drop=4*pi*r0*r0*wd/hns;
      
	w0jl=0.0;
	wdjl=0.0;
	sssv_jim=0.0;
      
      }
    
      //--------------------------------------------
      // Add columns to check table:
    
      if (true && qnn_equals_qnp==false && 
	  (model=="skyrme" || model=="hcskyrme")) {
	double tx;
	ubvector ty(5), tdy(5);
	if (!at.is_column("nnpp")) at.new_column("nnpp");
	if (!at.is_column("nppp")) at.new_column("nppp");
	if (!at.is_column("rhsn")) at.new_column("rhsn");
	if (!at.is_column("rhsp")) at.new_column("rhsp");
	if (!at.is_column("lhs_n")) at.new_column("lhs_n");
	if (!at.is_column("lhs_a")) at.new_column("lhs_a");
	if (!at.is_column("rhs_n")) at.new_column("rhs_n");
	if (!at.is_column("rhs_a")) at.new_column("rhs_a");
	at.set("nnpp",0,0.0);
	at.set("nppp",0,0.0);
	at.set("rhsn",0,0.0);
	at.set("rhsp",0,0.0);
	at.set("lhs_n",0,0.0);
	at.set("lhs_a",0,0.0);
	at.set("rhs_n",0,0.0);
	at.set("rhs_a",0,0.0);

	rhsmode=false;

	for(int i=1;i<((int)at.get_nlines());i++) {
	  
	  at.set("nnpp",i,(at.get("nnp",i)-at.get("nnp",i-1))/
		 (at.get("x",i)-at.get("x",i-1)));
	  at.set("nppp",i,(at.get("npp",i)-at.get("npp",i-1))/
		 (at.get("x",i)-at.get("x",i-1)));
	
	  tx=(at.get("x",i)+at.get("x",i-1))/2.0;
	  ty[0]=(at.get("nn",i)+at.get("nn",i-1))/2.0;
	  ty[1]=(at.get("np",i)+at.get("np",i-1))/2.0;
	  ty[2]=(at.get("nnp",i)+at.get("nnp",i-1))/2.0;
	  ty[3]=(at.get("npp",i)+at.get("npp",i-1))/2.0;

	  if (ty[1]<=0.0) rhsmode=true;
	  derivs(tx,ty,tdy);
	  
	  at.set("rhsn",i,tdy[2]);
	  at.set("rhsp",i,tdy[3]);
	  
	  if (ty[1]>0.0) {
	    at.set("lhs_n",i,(neutron.mu-mun+proton.mu-mup)/2.0);
	    at.set("lhs_a",i,(neutron.mu-mun-proton.mu+mup)/2.0);
	    at.set("rhs_n",i,(at.get("nnpp",i)+at.get("nppp",i))/2.0*
		   (qnn+qnp));
	    at.set("rhs_a",i,(at.get("nnpp",i)-at.get("nppp",i))/2.0*
		   (qnn-qnp));
	  } else {
	    at.set("lhs_n",i,neutron.mu-mun);
	    at.set("rhs_n",i,at.get("nnpp",i)*(qnn));
	    at.set("lhs_a",i,0.0);
	    at.set("rhs_a",i,0.0);
	  }
	}
      }

      hdf_file hf;
      string tablename=((string)"nr")+std::to_string(pf_index);
      hf.open_or_create("nr.o2");
      hdf_output(hf,at,tablename);
      hf.close();
      cout << "Wrote solution to file 'nr.o2'" << endl;
      cout << endl;

      // Loop for next proton fraction
    }
  
    return 0;
  }
  
};


int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  seminf_nr sn;
  sn.run(argc,argv);
  return 0;
}
