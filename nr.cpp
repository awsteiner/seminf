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
/**
   
   \warning This code is very preliminary, and needs quite a bit of
   work.
   
   \warning Does not work with dripped protons yet.
   
   \note The first proton fraction cannot be 0.5, as the 
   algorithm cannot handle later values different from 0.5
   
   12/15/03 - Implemented a new solution method for Qnn!=Qnp for when
   nndrip becomes > 0. It seems that in this case, it is better not to
   automatically adjust the scale of the x-axis, but use a gradual
   broadening similar to what we have done using the relativistic
   code. This allows calculation down to really low proton fractions!
   We need to see if this is consistent with the analytic code for
   both Qnn<Qnp (hard) and Qnn>Qnp (probably works).
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

#include "relax.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

//--------------------------------------------
// Global variables

static const double big=3.0;
static const int ngrid=100;

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

//--------------------------------------------
// Structure for relax parameters

class relaxp {

public:

  /// Neutron chemical potential
  double mun;
  /// Proton chemical potential
  double mup;
  double nn0;
  double np0;
  /// Isoscalar gradient term for neutrons
  double qnn;
  /// Isoscalar gradient term for protons
  double qpp;
  /// Isovector gradient term
  double qnp;
  /// Saturation density of isospin-symmetric matter
  double n0half;
  bool showrelax;
  double dndnl;
  double dndnr;
  double dndpl;
  double dndpr;
  double dpdnl;
  double dpdnr;
  double dpdpl;
  double dpdpr;
  double lex1;
  double lex2;
  double rex1;
  double rex2;
  double barn;
  double dmundn;
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
  fermion n;
  /// Proton
  fermion p;
  /// Thermodynamic functions
  thermo hb;
};

class seminf_nr *snp;

/** \brief Relaxation class
 */
class seminf_nr_relax : public relax {

public:

  /// Desc
  relaxp *rp;
  
  /// Desc
  int ctr;

  /// Desc
  double errold;
  
  /// Desc
  seminf_nr_relax(int tne, int tnb, int tngrid);
  
  /// Desc
  int iter(int k, double err, double fac, ubvector_int &kmax,
	   ubvector &errmax);

  /// Desc
  int difeq(int k, int k1, int k2, int jsf, int is1, int isf);
  
};

/** \brief Desc
 */
class seminf_nr {
  
public:
  
  seminf_nr() {
    rhslength=2.0;
  }

  /// If true, we have good boundary conditions
  bool goodbound;
  /// Initial derivative for densities on LHS (typically negative)
  double firstderiv;
  
  ubvector xstor;
  ubmatrix ystor;
  double xg1[ngrid+1];
  double xg2[ngrid+1];
  double yg1[6][ngrid+1];
  double yg2[6][ngrid+1];
  double nnprhs;
  double npprhs;
  double relaxconverge;
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
  /// \name Interactions to use
  //@{
  eos_had_apr eosa;
  eos_had_potential eosg;
  eos_had_skyrme eoss;
  eos_had_schematic eosp;
  //@}
  /// To compute derivatives
  deriv_gsl<> df;
  seminf_nr_relax *rel;
  relaxp rp;
  
  //--------------------------------------------
  // Function prototypes
  
  int run(int argc, char *argv[]) {

    snp=this;

    double wd=0.0, wd2=0.0, dtemp;
    double rhsn, rhsp, det;
    double surf, sbulk, sgrad, sssv_drop=0.0, *pflist;
    ubvector sx(5), sy(5);
    double sssv_jim=0.0, hns, delta, w0jl=0.0, wdjl=0.0, den, w0;
    // Locations of fixed relative densities
    double xn[4], xp[4];
    double thick;
    int pf_index, npf;
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

      rp.eos=&eoss;
      skyrme_load(eoss,"NRAPR");
      
      // Determine values of Q_{nn}, Q_{pp}, and Q_{np}
      rp.qnn=0.1875*(eoss.t1*(1.0-eoss.x1)-eoss.t2*(1.0+eoss.x2));
      rp.qpp=rp.qnn;
      rp.qnp=0.125*(3.0*eoss.t1*(1.0+eoss.x1/2.0)-
		    eoss.t2*(1.0+eoss.x2/2.0));
      
    } else if (model==((string)"schematic")) {
      rp.eos=&eosp;
    } else if (model==((string)"apr")) {
      rp.eos=&eosa;
    } else if (model==((string)"gp")) {
      rp.eos=&eosg;
    } else {
      cout << "Unknown EOS" << endl;
      exit(-1);
    }
  
    //--------------------------------------------
    // Test if Qnn=Qnp

    if (model!="gp" && model!="apr" && fabs(rp.qnn-rp.qnp)<1.0e-7) {
      qnn_equals_qnp=true;
      cout << "Qnn=Qnp" << endl;
      if (true) {
	rp.qnn/=1.01;
	rp.qpp/=1.01;
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
    rp.showrelax=true;
    relaxfile=false;
    firstderiv=-2.0e-4;
    rhsmin=1.0e-7;
    monfact=0.002;

    rp.n.init(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    rp.p.init(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_proton),2.0);
    rp.n.non_interacting=false;
    rp.p.non_interacting=false;
    rp.eos->set_n_and_p(rp.n,rp.p);
    rp.eos->set_thermo(rp.hb);

    // Compute saturation density
    rp.n0half=rp.eos->fn0(0.0,dtemp);
    double r0=cbrt(0.75/pi/rp.n0half);

    pflist=new double[2];
    pflist[0]=0.45;
    pflist[1]=0.46;
    npf=2;
  
    for(pf_index=1;pf_index<=npf;pf_index++) {
      rp.pf_index=pf_index;
      rp.protfrac=pflist[pf_index-1];

      //------------------------------------------------------
      // Calculate properties of saturated nuclear matter with
      // appropriate value of proton fraction including the
      // possibility of drip particles if they exist
      
      rp.nsat=rp.eos->fn0(1.0-2.0*rp.protfrac,dtemp);
      rp.np0=rp.nsat*rp.protfrac;
      rp.nn0=rp.nsat-rp.np0;
      rp.n.n=rp.nn0;
      rp.p.n=rp.np0;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      rp.mun=rp.n.mu;
      rp.mup=rp.p.mu;
      nndrip=0.0;
      npdrip=0.0;

      if (rp.mun>rp.n.m) {
	cout << "Neutron drip" << endl;

	sx[0]=rp.n.n;
	sx[1]=rp.p.n;
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

	rp.nn0=sx[0];
	rp.np0=sx[1];
	rp.n.n=rp.nn0;
	rp.p.n=rp.np0;
	rp.eos->calc_e(rp.n,rp.p,rp.hb);
	rp.mun=rp.n.mu;
	rp.mup=rp.p.mu;
	nndrip=sx[2];
	npdrip=0.0;

      }

      if (rp.mup>rp.p.m) {
	cout << "Proton drip" << endl;

	sx[0]=rp.n.n;
	sx[1]=rp.p.n;
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

	rp.nn0=sx[0];
	rp.np0=sx[1];
	rp.eos->calc_e(rp.n,rp.p,rp.hb);
	rp.mun=rp.n.mu;
	rp.mup=rp.p.mu;
	nndrip=0.0;
	npdrip=sx[2];

      }

      rp.nsat=rp.nn0+rp.np0;
      cout << "Saturation density at x=" << rp.protfrac << ": " << rp.nsat 
	   << endl;
      cout << "nn0: " << rp.nn0 << " np0: " << rp.np0 << endl;
      cout << "nndrip: " << nndrip << " npdrip: " << npdrip << endl;

      //--------------------------------------------
      // Calculate derivatives useful for the 
      // matching to the analytical solutions

      // dmundn
   
      rp.n.n=rp.nn0;
      rp.p.n=rp.np0;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      hin=rp.n.mu;

      rp.n.n=rp.nn0-1.0e-4;
      rp.p.n=rp.np0;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      lon=rp.n.mu;
    
      rp.dmundn=1.0e4*(hin-lon);
  
      // dmundp
  
      rp.n.n=rp.nn0;
      rp.p.n=rp.np0;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      hin=rp.n.mu;

      rp.n.n=rp.nn0;
      rp.p.n=rp.np0-1.0e-4;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      lon=rp.n.mu;
    
      rp.dmundp=1.0e4*(hin-lon);
  
      //--------------------------------------------
      // Calculate derivatives of the RHS's numerically
  
      // At the left hand side
      rp.n.n=rp.nn0;
      rp.p.n=rp.np0;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      if (model=="apr") {
	eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
	rhsn=(rp.n.mu-rp.mun)/rp.qpp/2.0;
	rhsp=(rp.p.mu-rp.mup)/rp.qnn/2.0;
      } else {
	rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	rhsn/=det;
	rhsp/=det;
      }
      hin=rhsn;
      hip=rhsp;

      rp.n.n=rp.nn0-1.0e-4;
      rp.p.n=rp.np0;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      if (model=="apr") {
	eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
	rhsn=(rp.n.mu-rp.mun)/rp.qpp/2.0;
	rhsp=(rp.p.mu-rp.mup)/rp.qnn/2.0;
      } else {
	rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;

      rp.dndnl=1.0e4*(hin-lon);
      rp.dpdnl=1.0e4*(hip-lop);
  
      rp.n.n=rp.nn0;
      rp.p.n=rp.np0-1.0e-4;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      if (model=="apr") {
	eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
	rhsn=(rp.n.mu-rp.mun)/rp.qpp/2.0;
	rhsp=(rp.p.mu-rp.mup)/rp.qnn/2.0;
      } else {
	rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;
  
      rp.dndpl=1.0e4*(hin-lon);
      rp.dpdpl=1.0e4*(hip-lop);
  
      // At the right hand side
      rp.n.n=1.0e-4;
      rp.p.n=1.0e-4;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      if (model=="apr") {
	eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
	rhsn=(rp.n.mu-rp.mun)/rp.qpp/2.0;
	rhsp=(rp.p.mu-rp.mup)/rp.qnn/2.0;
      } else {
	rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	rhsn/=det;
	rhsp/=det;
      }
      hin=rhsn;
      hip=rhsp;
  
      rp.n.n=1.0e-8;
      rp.p.n=1.0e-4;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      if (model=="apr") {
	eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
	rhsn=(rp.n.mu-rp.mun)/rp.qpp/2.0;
	rhsp=(rp.p.mu-rp.mup)/rp.qnn/2.0;
      } else {
	rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;
  
      rp.dndnr=1.0e4*(hin-lon);
      rp.dpdnr=1.0e4*(hip-lop);
  
      rp.n.n=1.0e-4;
      rp.p.n=1.0e-8;
      rp.eos->calc_e(rp.n,rp.p,rp.hb);
      if (model=="apr") {
	eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
      } 
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
      } 
      det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
      if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
	rhsn=(rp.n.mu-rp.mun)/rp.qpp/2.0;
	rhsp=(rp.p.mu-rp.mup)/rp.qnn/2.0;
      } else {
	rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	rhsn/=det;
	rhsp/=det;
      }
      lon=rhsn;
      lop=rhsp;
  
      rp.dndpr=1.0e4*(hin-lon);
      rp.dpdpr=1.0e4*(hip-lop);
  
      goodbound=true;

      sqt=sqrt(rp.dndnl*rp.dndnl-2.0*rp.dpdpl*rp.dndnl+rp.dpdpl*rp.dpdpl+
	       4.0*rp.dndpl*rp.dpdnl);
      rp.lex1=sqrt(2.0*(rp.dndnl+rp.dpdpl)+2.0*sqt)/2.0;
      rp.lex2=sqrt(2.0*(rp.dndnl+rp.dpdpl)-2.0*sqt)/2.0;

      if (!std::isfinite(rp.lex1)) {
	goodbound=false;
	rp.lex1=0.0;
      }
      if (!std::isfinite(rp.lex2)) {
	expo=rp.lex1;
	rp.lex2=0.0;
	goodbound=false;
      } else if (rp.lex2<rp.lex1) {
	expo=rp.lex2;
      } else {
	expo=rp.lex1;
      }

      sqt=sqrt(rp.dndnr*rp.dndnr-2.0*rp.dpdpr*rp.dndnr+rp.dpdpr*rp.dpdpr+
	       4.0*rp.dndpr*rp.dpdnr);
      rp.rex1=-sqrt(2.0*(rp.dndnr+rp.dpdpr)+2.0*sqt)/2.0;
      rp.rex2=-sqrt(2.0*(rp.dndnr+rp.dpdpr)-2.0*sqt)/2.0;
  
      if (!std::isfinite(rp.rex1)) {
	goodbound=false;
	rp.rex1=0.0;
      }
      if (!std::isfinite(rp.rex2)) {
	goodbound=false;
	rp.rex2=0.0;
      }

      //--------------------------------------------
      // Solve
      
      if (qnn_equals_qnp) {
	solve_qnn_equal_qnp(monfact);
      } else {
	solve_qnn_neq_qnp(monfact,argc);
      }

      //--------------------------------------------
      // Calculate the value of x for nn*0.9, nn*0.5, 
      // nn*0.1, etc.

      if (verbose>1) cout << "Lookups. " << endl;
      for(int i=1;i<=3;i++) {
	xn[i]=lookup(at.get_nlines(),
		     nndrip+(rp.nn0-nndrip)/10.0*((double)(i*4-3)),
		     at[0],at[1]);
	xp[i]=lookup(at.get_nlines(),
		     npdrip+(rp.np0-npdrip)/10.0*((double)(i*4-3)),
		     at[0],at[2]);
      }

      delta=1.0-2.0*rp.protfrac;
      rp.barn=rp.nsat;
      hns=rp.eos->fesym(rp.barn);

      if (verbose>1) cout << "Integrands." << endl;
      for(int i=0;i<((int)at.get_nlines());i++) {
	// Fix negative values
	if (at[1][i]<0.0) at.set(1,i,0.0);
	if (at[2][i]<0.0) at.set(2,i,0.0);

	rp.n.n=at[1][i];
	rp.p.n=at.get(2,i);

	rp.eos->calc_e(rp.n,rp.p,rp.hb);

	if (model=="apr") {
	  eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			     dqnndnn,dqnndnp,dqnpdnn,
			     dqnpdnp,dqppdnn,dqppdnp);
	  at.set("qnn",i,rp.qnn);
	  at.set("qnp",i,rp.qnp);
	  at.set("qpp",i,rp.qpp);
	}
	if (model=="gp") {
	  thermo thth;
	  eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			    dqnndnn,dqnndnp,dqnpdnn,
			    dqnpdnp,dqppdnn,dqppdnp);
	  at.set("qnn",i,rp.qnn);
	  at.set("qnp",i,rp.qnp);
	  at.set("qpp",i,rp.qpp);
	} 

	at.set("rho",i,rp.n.n+rp.p.n);
	at.set("alpha",i,rp.n.n-rp.p.n);
	at.set("ebulk",i,rp.hb.ed-rp.mun*at.get(1,i)-rp.mup*at.get(2,i));
	at.set("egrad",i,0.5*rp.qnn*at.get(3,i)*at.get(3,i)+
	       0.5*rp.qpp*at.get(4,i)*at.get(4,i)+
	       rp.qnp*at.get(3,i)*at.get(4,i));
	at.set("esurf",i,at.get("ebulk",i)+at.get("egrad",i));
	at.set("thickint",i,at.get(1,i)/rp.nn0-at.get(2,i)/rp.np0);
	if (fabs(rp.protfrac-0.5)>1.0e-8) {
	  at.set("wdint",i,at.get("alpha",i)/delta-at.get("rho",i));
	  if (at.get("rho",i)>=1.0e-5) {
	    rp.barn=at.get("rho",i);
	    at.set("wd2int",i,at.get("rho",i)*
		   (pow(at.get("alpha",i)/delta/at.get("rho",i),2.0)*
		    rp.eos->fesym(rp.barn)-hns));
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
      if (fabs(rp.protfrac-0.5)>1.0e-8) {
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
	  ty[1]=(at.get("nn",i)+at.get("nn",i-1))/2.0;
	  ty[2]=(at.get("np",i)+at.get("np",i-1))/2.0;
	  ty[3]=(at.get("nnp",i)+at.get("nnp",i-1))/2.0;
	  ty[4]=(at.get("npp",i)+at.get("npp",i-1))/2.0;

	  if (ty[2]<=0.0) rhsmode=true;
	  derivs(tx,ty,tdy);

	  at.set("rhsn",i,tdy[3]);
	  at.set("rhsp",i,tdy[4]);

	  if (ty[2]>0.0) {
	    at.set("lhs_n",i,(rp.n.mu-rp.mun+rp.p.mu-rp.mup)/2.0);
	    at.set("lhs_a",i,(rp.n.mu-rp.mun-rp.p.mu+rp.mup)/2.0);
	    at.set("rhs_n",i,(at.get("nnpp",i)+at.get("nppp",i))/2.0*
		   (rp.qnn+rp.qnp));
	    at.set("rhs_a",i,(at.get("nnpp",i)-at.get("nppp",i))/2.0*
		   (rp.qnn-rp.qnp));
	  } else {
	    at.set("lhs_n",i,rp.n.mu-rp.mun);
	    at.set("rhs_n",i,at.get("nnpp",i)*(rp.qnn));
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

  /** \brief Desc
   */
  double lookup(int n, double yy0, const std::vector<double> &x,
		const std::vector<double> &y) {

    for(int i=2;i<=n;i++) {
      if ((y[i]>=yy0 && y[i-1]<yy0) || (y[i]<yy0 && y[i-1]>=yy0)) {
	double x0=x[i-1]+(x[i]-x[i-1])*(yy0-y[i-1])/(y[i]-y[i-1]);
	return x0;
      }
    }

    return 0.0;
  }

  /** \brief Desc
   */
  int derivs(double sx, const ubvector &sy, ubvector &dydx) {
    double rhsn, rhsp=0.0, det;
    double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

    rp.n.n=sy[1];
    rp.p.n=sy[2];
    if (model=="apr") {
      eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
    }
    if (model=="gp") {
      thermo thth;
      if (false) {
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
      } else {
	rp.qnn=100.0/hc_mev_fm;
	rp.qpp=100.0/hc_mev_fm;
	rp.qnp=90.0/hc_mev_fm;
	dqnndnn=0.0;
	dqnndnp=0.0;
	dqnpdnn=0.0;
	dqnpdnp=0.0;
	dqppdnn=0.0;
	dqppdnp=0.0;
      }
    }

    rp.eos->calc_e(rp.n,rp.p,rp.hb);
#ifdef DEBUG
    cout << "A: " << rp.n.n << " " << rp.p.n << " " << rp.qnn << endl;
    cout << "B: " << rp.n.m << " " << rp.p.m << endl;
    cout << "C: " << rp.n.mu << " " << rp.p.mu << endl;
#endif
    /*
      double mun=rp.n.mu, mup=rp.p.mu, msn=rp.n.ms, msp=rp.p.ms;
      double hed=rp.hb->ed, hpr=rp.hb->pr, part;
      if (rp.pf_index==1) part=1000.0;
      else if (rp.pf_index==2) part=100.0;
      else if (rp.pf_index==3) part=40.0;
      else if (rp.pf_index==4) part=20.0;
      else if (rp.pf_index==5) part=5.0;
      else if (rp.pf_index==6) part=3.0;
      else if (rp.pf_index==7) part=1.0;
      else part=0.0;
      eosa.calc_e(rp.n,rp.p,rp.hb);
    
      rp.n.mu=(rp.n.mu*part+mun)/(part+1.0);
      rp.p.mu=(rp.p.mu*part+mup)/(part+1.0);
      rp.n.ms=(rp.n.ms*part+msn)/(part+1.0);
      rp.p.ms=(rp.p.ms*part+msp)/(part+1.0);
      rp.hb->ed=(rp.hb->ed*part+hed)/(part+1.0);
      rp.hb->pr=(rp.hb->pr*part+hpr)/(part+1.0);
    */

    if (rp.n.n<=0.0) {
      if (rp.p.n<=0.0) {
	dydx[1]=0.0;
	dydx[2]=0.0;
	dydx[3]=0.0;
	dydx[4]=0.0;
      } else {
	dydx[1]=0.0;
	dydx[2]=sy[4];
	dydx[3]=0.0;
	if (model!="apr" && model!="gp") {
	  dydx[4]=(rp.p.mu-rp.mup)/rp.qpp;
	} else {
	  dydx[4]=(rp.p.mu-rp.mup+0.5*dqppdnp*sy[4]*sy[4])/rp.qpp;
	}
      }
    } else if (rp.p.n<=0.0) {
      dydx[1]=sy[3];
      dydx[2]=0.0;
      if (model!="apr" && model!="gp") {
	dydx[3]=(rp.n.mu-rp.mun)/rp.qnn;
      } else {
	dydx[3]=(rp.n.mu-rp.mun+0.5*dqnndnn*sy[3]*sy[3])/rp.qnn;
      }
      dydx[4]=0.0;
    } else {
      if (rhsmode==false) {
	if (model!="apr" && model!="gp") {
	  det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
	  rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	  rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	  rhsn/=det;
	  rhsp/=det;
#ifdef DEBUG
	  cout << "K: " << rhsn << " " << rhsp << " " << rp.n.mu << " "
	       << rp.qpp << endl;
	  cout << "L: " << rp.qnp << " " << rp.qnn << " "
	       << rp.p.mu << " " << rp.mun << endl;
	  #endif
	} else {
	  det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
	  rhsn=(rp.n.mu-rp.mun-0.5*dqnndnn*sy[3]*sy[3]+dqnndnp*sy[3]*sy[4]-
		dqnpdnp*sy[4]*sy[4]+0.5*dqppdnn*sy[4]*sy[4])*rp.qpp-
	    (rp.p.mu-rp.mup+0.5*dqnndnp*sy[3]*sy[3]-dqnpdnn*sy[3]*sy[3]-
	     dqppdnn*sy[3]*sy[4]-0.5*dqppdnp*sy[4]*sy[4])*rp.qnp;
	  rhsp=(rp.p.mu-rp.mup+0.5*dqnndnp*sy[3]*sy[3]-dqnpdnn*sy[3]*sy[3]-
		dqppdnn*sy[3]*sy[4]-0.5*dqppdnp*sy[4]*sy[4])*rp.qnn-
	    (rp.n.mu-rp.mun-0.5*dqnndnn*sy[3]*sy[3]+dqnndnp*sy[3]*sy[4]-
	     dqnpdnp*sy[4]*sy[4]+0.5*dqppdnn*sy[4]*sy[4])*rp.qnp;
	  rhsn/=det;
	  rhsp/=det;
	}
      } else {
	if (model!="apr" && model!="gp") {
	  rhsn=(rp.n.mu-rp.mun)/rp.qnn;
	  rhsp=0.0;
	} else {
	  rhsn=(rp.n.mu-rp.mun+0.5*dqnndnn*sy[3]*sy[3])/rp.qnn;
	  rhsp=0.0;
	}
      }
    
      dydx[1]=sy[3];
      dydx[2]=sy[4];
      dydx[3]=rhsn;
      dydx[4]=rhsp;
    }
  
    //--------------------------------------------
    // Return sensible results 

    if (!std::isfinite(dydx[3])) {
      cout << "3 not finite." << endl;
      dydx[3]=0.0;
    }
    if (!std::isfinite(dydx[4])) {
      cout << "4 not finite." << endl;
      dydx[4]=0.0;
    }
  
    return 0;
  }

  /** \brief Desc
   */
  int derivs2(double sx, const ubvector &sy, ubvector &dydx) {
    double rhsn, rhsp=0.0, det;
    double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

    rp.n.n=sy[0];
    rp.p.n=sy[1];
    if (model=="apr") {
      eosa.gradient_qij2(rp.n.n,rp.p.n,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
    }
    if (model=="gp") {
      thermo thth;
      if (false) {
	eosg.gradient_qij(rp.n,rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
      } else {
	rp.qnn=100.0/hc_mev_fm;
	rp.qpp=100.0/hc_mev_fm;
	rp.qnp=90.0/hc_mev_fm;
	dqnndnn=0.0;
	dqnndnp=0.0;
	dqnpdnn=0.0;
	dqnpdnp=0.0;
	dqppdnn=0.0;
	dqppdnp=0.0;
      }
    }

    rp.eos->calc_e(rp.n,rp.p,rp.hb);
#ifdef DEBUG
    cout << "A2: " << rp.n.n << " " << rp.p.n << " " << rp.qnn << endl;
    cout << "B2: " << rp.n.m << " " << rp.p.m << endl;
    cout << "C2: " << rp.n.mu << " " << rp.p.mu << endl;
    #endif
    /*
      double mun=rp.n.mu, mup=rp.p.mu, msn=rp.n.ms, msp=rp.p.ms;
      double hed=rp.hb->ed, hpr=rp.hb->pr, part;
      if (rp.pf_index==1) part=1000.0;
      else if (rp.pf_index==2) part=100.0;
      else if (rp.pf_index==3) part=40.0;
      else if (rp.pf_index==4) part=20.0;
      else if (rp.pf_index==5) part=5.0;
      else if (rp.pf_index==6) part=3.0;
      else if (rp.pf_index==7) part=1.0;
      else part=0.0;
      eosa.calc_e(rp.n,rp.p,rp.hb);
    
      rp.n.mu=(rp.n.mu*part+mun)/(part+1.0);
      rp.p.mu=(rp.p.mu*part+mup)/(part+1.0);
      rp.n.ms=(rp.n.ms*part+msn)/(part+1.0);
      rp.p.ms=(rp.p.ms*part+msp)/(part+1.0);
      rp.hb->ed=(rp.hb->ed*part+hed)/(part+1.0);
      rp.hb->pr=(rp.hb->pr*part+hpr)/(part+1.0);
    */

    if (rp.n.n<=0.0) {
      if (rp.p.n<=0.0) {
	dydx[0]=0.0;
	dydx[1]=0.0;
	dydx[2]=0.0;
	dydx[3]=0.0;
      } else {
	dydx[0]=0.0;
	dydx[1]=sy[3];
	dydx[2]=0.0;
	if (model!="apr" && model!="gp") {
	  dydx[3]=(rp.p.mu-rp.mup)/rp.qpp;
	} else {
	  dydx[3]=(rp.p.mu-rp.mup+0.5*dqppdnp*sy[3]*sy[3])/rp.qpp;
	}
      }
    } else if (rp.p.n<=0.0) {
      dydx[0]=sy[2];
      dydx[1]=0.0;
      if (model!="apr" && model!="gp") {
	dydx[2]=(rp.n.mu-rp.mun)/rp.qnn;
      } else {
	dydx[2]=(rp.n.mu-rp.mun+0.5*dqnndnn*sy[2]*sy[2])/rp.qnn;
      }
      dydx[3]=0.0;
    } else {
      if (rhsmode==false) {
	if (model!="apr" && model!="gp") {
	  det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
	  rhsn=(rp.n.mu-rp.mun)*rp.qpp-(rp.p.mu-rp.mup)*rp.qnp;
	  rhsp=(rp.p.mu-rp.mup)*rp.qnn-(rp.n.mu-rp.mun)*rp.qnp;
	  rhsn/=det;
	  rhsp/=det;
#ifdef DEBUG
	  cout << "K2: " << rhsn << " " << rhsp << " " << rp.n.mu << " "
	       << rp.qpp << endl;
	  cout << "L2: " << rp.qnp << " " << rp.qnn << " "
	       << rp.p.mu << " " << rp.mun << endl;
	  #endif
	} else {
	  det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
	  rhsn=(rp.n.mu-rp.mun-0.5*dqnndnn*sy[2]*sy[2]+dqnndnp*sy[2]*sy[3]-
		dqnpdnp*sy[3]*sy[3]+0.5*dqppdnn*sy[3]*sy[3])*rp.qpp-
	    (rp.p.mu-rp.mup+0.5*dqnndnp*sy[2]*sy[2]-dqnpdnn*sy[2]*sy[2]-
	     dqppdnn*sy[2]*sy[3]-0.5*dqppdnp*sy[3]*sy[3])*rp.qnp;
	  rhsp=(rp.p.mu-rp.mup+0.5*dqnndnp*sy[2]*sy[2]-dqnpdnn*sy[2]*sy[2]-
		dqppdnn*sy[2]*sy[3]-0.5*dqppdnp*sy[3]*sy[3])*rp.qnn-
	    (rp.n.mu-rp.mun-0.5*dqnndnn*sy[2]*sy[2]+dqnndnp*sy[2]*sy[3]-
	     dqnpdnp*sy[3]*sy[3]+0.5*dqppdnn*sy[3]*sy[3])*rp.qnp;
	  rhsn/=det;
	  rhsp/=det;
	}
      } else {
	if (model!="apr" && model!="gp") {
	  rhsn=(rp.n.mu-rp.mun)/rp.qnn;
	  rhsp=0.0;
	} else {
	  rhsn=(rp.n.mu-rp.mun+0.5*dqnndnn*sy[2]*sy[2])/rp.qnn;
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

  /** \brief Desc
   */
  double solve_qnn_neq_qnp(double lmonfact, int argc) {
    bool guessdone;
    int debug=0;
    double dx, xrhs, delta, epsi;
    ubvector y(5), dydx(5);
    int ilast=0, interpi;
    ofstream itout;
    char ch;

    monfact=lmonfact;

    //----------------------------------------------
    // Create object

    if (flatden) rel=new seminf_nr_relax(5,3,ngrid);
    else rel=new seminf_nr_relax(5,4,ngrid);
    rel->itmax=100;
    rel->rp=&rp;

    if (rp.pf_index==1) {

      //----------------------------------------------
      // Construct a guess by shooting
    
      xstor[1]=0.0;
      ystor(1,1)=rp.nn0;
      ystor(2,1)=rp.np0;
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
      
	for(int j=1;j<=4;j++) {
	  y[j]=ystor(j,i-1);
	}
	derivs(xstor[i-1],y,dydx);
	for(int j=1;j<=4;j++) ystor(j,i)=ystor(j,i-1)+dx*dydx[j];
      
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
	for(int j=1;j<=4;j++) {
	  rel->y(j,i)=ystor(j,interpi)+
	    (ystor(j,interpi+1)-ystor(j,interpi))*
	    (rel->x[i]-xstor[interpi])/(xstor[interpi+1]-xstor[interpi]);
	}
	rel->y(5,i)=(rel->x[ngrid]);
	rel->x[i]=((double)(i-1))/((double)(ngrid-1));
      }

      if (debug>1) {
	cout.precision(4);
	for(int iz=1;iz<=ngrid;iz++) {
	  cout.width(3);
	  cout << iz << " " << rel->x[iz] << " " 
	       << rel->y(1,iz) << " " 
	       << rel->y(2,iz) << " " 
	       << rel->y(3,iz) << " " 
	       << rel->y(4,iz) << " " 
	       << rel->y(5,iz) << endl;
	}
	cout.precision(6);
      }
    
    } else {

      //----------------------------------------------
      // Use last solution for guess
    
      for(int i=1;i<=ngrid;i++) {
	rel->x[i]=xg1[i];
	for(int j=1;j<=rel->ne;j++) {
	  rel->y(j,i)=yg1[j][i];
	}
      }

    }

    //----------------------------------------------
    // Try ode_it_solve

    if (argc>=2) {
      ubvector ox(ngrid);
      ubmatrix oy(ngrid,rel->ne);
      for(int i=1;i<=ngrid;i++) {
	ox[i-1]=rel->x[i];
	for(int j=1;j<=rel->ne;j++) {
	  oy(i-1,j-1)=rel->y(j,i);
	}
      }
      ode_it_funct11 f_derivs=std::bind
	(std::mem_fn<double(size_t,double,ubmatrix_row &)>
	 (&seminf_nr::difeq2),this,std::placeholders::_1,std::placeholders::_2,
	 std::placeholders::_3);       
      ode_it_funct11 f_left=std::bind
	(std::mem_fn<double(size_t,double,ubmatrix_row &)>
	 (&seminf_nr::left2),this,std::placeholders::_1,std::placeholders::_2,
	 std::placeholders::_3);       
      ode_it_funct11 f_right=std::bind
	(std::mem_fn<double(size_t,double,ubmatrix_row &)>
	 (&seminf_nr::right2),this,std::placeholders::_1,std::placeholders::_2,
	 std::placeholders::_3);   
      ode_it_solve<> oit;
      ubmatrix A(ngrid*rel->ne,ngrid*rel->ne);
      ubvector rhs(ngrid*rel->ne), dy(ngrid*rel->ne);
      oit.verbose=2;
      oit.solve(ngrid,rel->ne,rel->nb,ox,oy,f_derivs,f_left,f_right,
		A,rhs,dy);
      exit(-1);
    }

    //----------------------------------------------
    // Solve

    rhsmode=false;
    int relret1=rel->solve(relaxconverge,1.0);
  
    //----------------------------------------------
    // Copy solution for next guess

    for(int i=1;i<=ngrid;i++) {
      xg1[i]=rel->x[i];
      for(int j=1;j<=rel->ne;j++) {
	yg1[j][i]=rel->y(j,i);
      }
    }
  
    //----------------------------------------------
    // Rescale solution

    for(int i=1;i<=rel->ngrid;i++) {
      rel->x[i]=rel->x[i]*rel->y(5,i);
    }

    //----------------------------------------------
    // Store quantites at rhs for later use

    xrhs=rel->x[ngrid];
    if (npdrip==0.0) {
      nnrhs=rel->y(1,ngrid);
      cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
      if (nnrhs<rhsmin) {
	nnrhs=rhsmin;
	cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
      }
    } else {
      nprhs=rel->y(2,ngrid);
      cout << "RHS, nprhs=" << nprhs << endl;
      if (nprhs<rhsmin) {
	nprhs=rhsmin;
	cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
      }
    }

    //----------------------------------------------
    // Store solution and delete seminf_nr_relax

    for(int i=1;i<=ngrid;i++) {
      at.set(0,i-1,rel->x[i]);
      for(int j=1;j<=5;j++) {
	at.set(j,i-1,rel->y(j,i));
      }
    }
    delete rel;
    at.set_nlines(ngrid);

    //----------------------------------------------
    // RHS
    
    if (fabs(rp.protfrac-0.5)>1.0e-4) {

      at.set_nlines(ngrid*2-1);

      rel=new seminf_nr_relax(3,1,ngrid);
      rel->itmax=100;
      rel->rp=&rp;

      if (rp.pf_index==1) {
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
	  for(int j=1;j<=3;j++) {
	    rel->y(j,i)=yg2[j][i];
	  }
	}

      }

      //----------------------------------------------
      // Solve
  
      rhsmode=true;
      cout << "second solve" << endl;

      if (nndrip>0.0) {
	cout << "Rhs length: " << rhslength << endl;
	int relret2=rel->solve(relaxconverge,1.0);

	while (-rel->y(2,ngrid)>1.0e-4) {
	  rhslength*=1.2;

	  cout << "Rhs length: " << rhslength << endl;
	  relret2=rel->solve(relaxconverge,1.0);
	}
      } else {
	int relret2=rel->solve(relaxconverge,1.0);
      }
    
      //----------------------------------------------
      // Copy solution for next guess

      for(int i=1;i<=ngrid;i++) {
	xg2[i]=rel->x[i];
	for(int j=1;j<=3;j++) {
	  yg2[j][i]=rel->y(j,i);
	}
      }
  
      //----------------------------------------------
      // Store solution and delete seminf_nr_relax

      for(int i=1;i<=ngrid;i++) {
	at.set(0,i+ngrid-2,rel->x[i]);
	for(int j=1;j<=3;j++) {
	  at.set(j,i+ngrid-2,rel->y(j,i));
	}
      }
      delete rel;

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
	for(int i=1;i<=rel->ngrid;i++) {
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

    return 0.0;

  }

  /** \brief Desc
   */
  double solve_qnn_equal_qnp(double lmonfact) {
    bool guessdone, debug=false;
    double dx, y[5], dydx[5], xrhs;
    int ilast=0, interpi;
    ofstream itout;
    char ch;
    ubvector sx(3), sy(3);

    monfact=lmonfact;

    //----------------------------------------------
    // Create object

    rel=new seminf_nr_relax(3,2,ngrid);
    rel->itmax=200;
    rel->rp=&rp;

    //----------------------------------------------
    // Construct a guess by shooting

    if (rp.pf_index==1) {
      xstor[1]=0.0;
      ystor(1,1)=rp.nn0+rp.np0;
      ystor(2,1)=firstderiv;
      dx=initialstep/((double)ngrid);
    
      guessdone=false;
      for(int i=2;guessdone==false && i<=ngrid;i++) {
	xstor[i]=xstor[i-1]+dx;
      
	sx[1]=(1.0-rp.protfrac)*ystor(1,i-1);
	sx[2]=rp.protfrac*ystor(1,i-1);
	rp.barn=ystor(1,i-1);
	mm_funct11 qqf=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	   (&seminf_nr::qnnqnpfun),snp,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	nd.msolve(2,sx,qqf);
      
	ystor(1,i)=ystor(1,i-1)+dx*ystor(2,i-1);
	ystor(2,i)=ystor(2,i-1)+dx*(rp.n.mu-rp.mun)/(rp.qnn);
      
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
	for(int j=1;j<=rel->ne;j++) {
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
      for(int j=1;j<=rel->ne;j++) {
	yg1[j][i]=rel->y(j,i);
      }
    }
  
    //----------------------------------------------
    // Rescale solution and set neutron and proton
    // densities:

    at.set_nlines(ngrid);
    for(int i=1;i<=rel->ngrid;i++) {

      sx[1]=(1.0-rp.protfrac)*rel->y(1,i);
      sx[2]=rp.protfrac*rel->y(1,i);

      rp.barn=rel->y(1,i);
      mm_funct11 qqf=std::bind
	(std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	 (&seminf_nr::qnnqnpfun),snp,std::placeholders::_1,
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
    
    if (fabs(rp.protfrac-0.5)>1.0e-4) {
      rel=new seminf_nr_relax(3,1,ngrid);
      rel->rp=&rp;
    
      if (rp.pf_index==1) {
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
	  for(int j=1;j<=rel->ne;j++) {
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
	for(int j=1;j<=rel->ne;j++) {
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

  /** \brief Desc
   */
  int qnnqnpfun(size_t sn, const ubvector &sx, ubvector &sy) {

    rp.n.n=sx[1];
    rp.p.n=sx[2];
    if (sx[1]<0.0 || sx[2]<0.0) return 1;
  
    rp.eos->calc_e(rp.n,rp.p,rp.hb);

    sy[1]=rp.n.n+rp.p.n-rp.barn;
    sy[2]=rp.n.mu-rp.p.mu-rp.mun+rp.mup;

    return 0;
  }

  /** \brief Desc
   */
  int ndripfun(size_t sn, const ubvector &sx, ubvector &sy) {
    double pleft, pright, munleft, munright;

    if (sx[1]<0.0 || sx[2]<0.0 || sx[3]<0.0) return 1;

    rp.n.n=sx[1];
    rp.p.n=sx[2];
    sy[1]=rp.p.n-rp.protfrac*(rp.p.n+rp.n.n);
  
    rp.eos->calc_e(rp.n,rp.p,rp.hb);
    pleft=rp.hb.pr;
    munleft=rp.n.mu;

    rp.n.n=sx[3];
    rp.p.n=0.0;
  
    rp.eos->calc_e(rp.n,rp.p,rp.hb);
    pright=rp.hb.pr;
    munright=rp.n.mu;

    sy[2]=pleft-pright;
    sy[3]=munleft-munright;

    return 0;
  }

  /** \brief Desc
   */
  int pdripfun(size_t sn, const ubvector &sx, ubvector &sy) {
    double pleft, pright, mupleft, mupright;

    if (sx[1]<0.0 || sx[2]<0.0 || sx[3]<0.0) return 1;

    rp.n.n=sx[1];
    rp.p.n=sx[2];
    sy[1]=rp.p.n-rp.protfrac*(rp.p.n+rp.n.n);
  
    rp.eos->calc_e(rp.n,rp.p,rp.hb);
    pleft=rp.hb.pr;
    mupleft=rp.p.mu;

    rp.n.n=0.0;
    rp.p.n=sx[3];
  
    rp.eos->calc_e(rp.n,rp.p,rp.hb);
    pright=rp.hb.pr;
    mupright=rp.p.mu;

    sy[2]=pleft-pright;
    sy[3]=mupleft-mupright;

    return 0;
  }

  /** \brief Future function for \ref o2scl::ode_it_solve
   */
  double difeq2(size_t ieq, double x, ubmatrix_row &y) {
    ubvector y2=y, dydx(5);
#ifdef DEBUG
    cout << "H2: " << y2(0) << " " << y2(1) << " " << y2(2) << " "
	 << y2(3) << " " << y2(4) << endl;
    #endif
    derivs2(x,y2,dydx);
#ifdef DEBUG
    cout << "J2: " << dydx(0) << " " << dydx(1) << " " << dydx(2) << " "
	 << dydx(3) << " " << dydx(4) << endl;
    char ch;
    cin >> ch;
    #endif
    return dydx[ieq]*y2[4];
  }

  /** \brief Future function for \ref o2scl::ode_it_solve
   */
  double left2(size_t ieq, double x, ubmatrix_row &y) {
    
    double delta=rp.nn0*monfact/exp(big*expo);
    double epsi=delta*(-rp.dmundn+rp.qnn*expo*expo)/
      (rp.dmundp-rp.qnp*expo*expo);

    if (epsi*delta<0.0) {
      cout << "epsi and delta have different signs" << endl;
      exit(-1);
    }

    /*
      cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
      cout << rp.nn0 << " " << rp.np0 << " "
      << delta << " " << epsi << " " << big << " "
      << expo << endl;
    */
    
    double ret;
    if (ieq==0) ret=y[0]-(rp.nn0-delta*exp(big*expo));
    else if (ieq==1) ret=y[1]-(rp.np0-epsi*exp(big*expo));
    else if (ieq==2) ret=y[2]+delta*expo*exp(big*expo);
    else ret=y[3]+epsi*expo*exp(big*expo);
    //cout << "I " << ieq << " " << ret << endl;
    return ret;
  }

  /** \brief Future function for \ref o2scl::ode_it_solve
   */
  double right2(size_t ieq, double x, ubmatrix_row &y) {
    return y[1];
  }
  
};

// Constructor
seminf_nr_relax::seminf_nr_relax(int tne, int tnb, int tngrid) :
  relax(tne,tnb,tngrid) {
}

// Code to be executed every iteration
int seminf_nr_relax::iter(int k, double err, double fac, ubvector_int &kmax,
			  ubvector &ermax) {
  
  ofstream itout;
  char s1[80];

  //--------------------------------------------
  // Output relaxation steps to file
  
  if (relaxfile) {
    string soutt="rel"+itos(k)+".out"+
      itos(rp->pf_index);
    itout.open(soutt.c_str());
    itout.setf(ios::scientific);
    for(int i=1;i<=ngrid;i++) {
      itout << x[i] << " ";
      for(int j=1;j<=ne;j++) itout << (y)(j,i) << " ";
      itout << endl;
    }
    itout.close();
  }

  //--------------------------------------------
  // Track relaxation on cout
  
  if (rp->showrelax) {
    if (k==1) {
      cout << "\tIter.\tError\t\tDer           Center       rhs" << endl;
    }
    if (k>=1 && k<=itmax) {
      cout << "\t" << k << "  \t" << err << "  \t" << y(3,1) << " ";
      ubmatrix_row amc(y,1);
      std::vector<double> amc2;
      std::vector<double> x2;
      vector_copy(amc,amc2);
      vector_copy(x,x2);
      cout << snp->lookup(ngrid,y(1,1)/2.0,x2,amc2) << " " 
	   << x[ngrid] << endl;
    }
  }
  
  //--------------------------------------------
  // Correct if density is negative
  
  if (!rhsmode && !qnn_equals_qnp) {
    for(int i=1;i<=ngrid;i++) {
      if ((this->y)(1,i)<0.0) (this->y)(3,i)=0.0;
      if ((this->y)(2,i)<0.0) (this->y)(4,i)=0.0;
    }
  }

  //--------------------------------------------
  // Notice if we're not making any progress

  if (k==1) {
    ctr=0;
    errold=0.0;
  } else if (fabs(err-errold)/fabs(err)<1.0e-4) {
    ctr++;
  }
  
  //--------------------------------------------
  // Check too see if we have failed

  if (k>=itmax) {
    cout << "Too many iterations." << endl;
    return 10;
  } else if (err>1.0e4) {
    cout << "Error exceeded maximum." << endl;
    return 11;
  } else if (ctr>=5) {
    if (err<minerr) {
      cout << "Relaxation making no progress. \"err\" is small." << endl;
      return 12;
    } else {
      cout << "Relaxation making no progress. \"err\" is large." << endl;
      return 13;
    }
  } else {
    errold=err;
  }

  return 0;
}

/** \brief Desc
 */
int seminf_nr_relax::difeq(int k, int k1, int k2, int jsf, int is1,
			   int isf) {

  double dx, exx;
  ubvector exy(8), exdy(8);
  ubvector sx(4);
  double delta, epsi, tmp;

  //--------------------------------------------
  // LHS for Qnn==Qnp

  if (qnn_equals_qnp==true && rhsmode==false) {

    delta=(rp->nn0+rp->np0)*monfact/exp(big*expo);
    mm_funct11 qqf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&seminf_nr::qnnqnpfun),snp,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);

    if (k==k1) {
      s(2,jsf)=y(1,1)-(rp->nn0+rp->np0-delta*exp(big*expo));
      s(3,jsf)=y(2,1)+delta*expo*exp(big*expo);
    } else if (k>k2 && npdrip==0.0) {
      sx[1]=(1.0-rp->protfrac)*y(1,ngrid);
      sx[2]=rp->protfrac*y(1,ngrid);
      rp->barn=y(1,ngrid);
      nd.msolve(2,sx,qqf);
      s(1,jsf)=sx[2];
    } else if (k>k2 && nndrip==0.0) {
      sx[1]=(1.0-rp->protfrac)*y(1,ngrid);
      sx[2]=rp->protfrac*y(1,ngrid);
      rp->barn=y(1,ngrid);
      nd.msolve(2,sx,qqf);
      s(1,jsf)=sx[1];
    } else {
      dx=x[k]-x[k-1];
      for(int i=1;i<=3;i++) exy[i]=(y(i,k)+y(i,k-1))/2.0;

      sx[1]=(1.0-rp->protfrac)*exy[1];
      sx[2]=rp->protfrac*exy[1];
      rp->barn=exy[1];
      nd.msolve(2,sx,qqf);
      
      if (!std::isfinite(rp->n.mu)) {
	rp->n.mu=rp->n.m;
      }
      if (!std::isfinite(rp->p.mu)) {
	rp->p.mu=rp->p.m;
      }

      s(1,jsf)=y(1,k)-y(1,k-1)-dx*exy[3]*exy[2];
      s(2,jsf)=y(2,k)-y(2,k-1)-
	dx*exy[3]*(rp->n.mu-rp->mun)/rp->qnn;
      s(3,jsf)=y(3,k)-y(3,k-1);
    }

    //--------------------------------------------
    // LHS for Qnn!=Qnp for both flatden==true and
    // flatden==false and for APR

  } else if (rhsmode==false) {

    if (k==k1) {
      delta=rp->nn0*monfact/exp(big*expo);
      epsi=delta*(-rp->dmundn+rp->qnn*expo*expo)/
	(rp->dmundp-rp->qnp*expo*expo);
      
      if (epsi*delta<0.0) {
	cout << "epsi and delta have different signs" << endl;
	exit(-1);
      }

      if (flatden) {
	s(3,jsf)=y(1,1)-(rp->nn0-delta*exp(big*expo));
	s(4,jsf)=y(2,1)-(rp->np0-epsi*exp(big*expo));
	s(5,jsf)=y(3,1)+delta*expo*exp(big*expo);
      } else {
	s(2,jsf)=y(1,1)-(rp->nn0-delta*exp(big*expo));
	s(3,jsf)=y(2,1)-(rp->np0-epsi*exp(big*expo));
	s(4,jsf)=y(3,1)+delta*expo*exp(big*expo);
	s(5,jsf)=y(4,1)+epsi*expo*exp(big*expo);
	/*
	cout << y(1,1) << " " << y(2,1) << " "
	<< y(3,1) << " " << y(4,1) << endl;
	cout << rp->nn0 << " " << rp->np0 << " "
	<< delta << " " << epsi << " " << big << " "
	<< expo << endl;
	cout << "I " << 0 << " " << s(2,jsf) << endl;
	cout << "I " << 0 << " " << s(3,jsf) << endl;
	cout << "I " << 0 << " " << s(4,jsf) << endl;
	cout << "I " << 0 << " " << s(5,jsf) << endl;
	exit(-1);
	*/
      }
      
    } else if (k>k2) {
      if (npdrip==0.0) {
	if (flatden) {
	  s(1,jsf)=y(2,ngrid);
	  s(2,jsf)=y(4,ngrid);
	} else {
	  s(1,jsf)=y(2,ngrid);
	}
      } else if (nndrip==0.0) {
	if (flatden) {
	  s(1,jsf)=y(1,ngrid);
	  s(2,jsf)=y(3,ngrid);
	} else {
	  s(1,jsf)=y(1,ngrid);
	}
      }
    } else {
      dx=x[k]-x[k-1];
      
      for(int i=1;i<=5;i++) exy[i]=(y(i,k)+y(i,k-1))/2.0;
      exx=(x[k]+x[k-1])/2.0*exy[5];

#ifdef DEBUG
      cout << "H: ";
      for(int i=1;i<=5;i++) {
	cout << exy[i] << " ";
      }
      cout << endl;
      #endif
      
      snp->derivs(exx,exy,exdy);
      exdy[5]=0.0;

#ifdef DEBUG
      cout << "J: ";
      for(int i=1;i<=5;i++) {
	cout << exdy[i] << " ";
      }
      cout << endl;

      exit(-1);
      #endif

      for(int i=1;i<=5;i++) {
	s(i,jsf)=y(i,k)-y(i,k-1)-dx*exy[5]*exdy[i];
      }

    }

    //--------------------------------------------
    // RHS for both Qnn==Qnp and Qnn!=Qnp and APR

  } else {

    if (k==k1) {
      if (npdrip==0.0) s(3,jsf)=y(1,1)-nnrhs;
      else s(3,jsf)=y(1,1)-nprhs;
    } else if (k>k2) {
      if (npdrip==0.0) {
	s(1,jsf)=y(1,ngrid)-nndrip;
	if (nndrip>0.0) {
	  s(2,jsf)=y(3,ngrid)-rhslength;
	} else {
	  s(2,jsf)=y(2,ngrid);
	}
      } else {
	s(1,jsf)=y(1,ngrid)-npdrip;
	s(2,jsf)=y(2,ngrid);
      }
    } else {
      dx=x[k]-x[k-1];
    
      exy[1]=(y(1,k)+y(1,k-1))/2.0;
      exy[2]=0.0;
      exy[3]=(y(2,k)+y(2,k-1))/2.0;
      exy[4]=0.0;
      exy[5]=(y(3,k)+y(3,k-1))/2.0;
      exx=(x[k]+x[k-1])/2.0*exy[5];
      snp->derivs(exx,exy,exdy);

      // 10/6/04
      // I think this statement was missing in previous versions
      // and caused nasty things to happen. This should be
      // double-checked. FIXME!
      exdy[5]=0.0;

      s(1,jsf)=y(1,k)-y(1,k-1)-dx*exy[5]*exdy[1];
      s(2,jsf)=y(2,k)-y(2,k-1)-dx*exy[5]*exdy[3];
      s(3,jsf)=y(3,k)-y(3,k-1)-dx*exy[5]*exdy[5];
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  seminf_nr sn;
  sn.run(argc,argv);
  return 0;
}
