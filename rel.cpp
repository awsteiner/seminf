/*
  -------------------------------------------------------------------
  
  Copyright (C) 2002-2017, Andrew W. Steiner
  
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

#include <o2scl/eos_had_rmf.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

/** \brief Semi-infinite nuclear matter for relativistic mean-field
    models in the Thomas-Fermi approximation

    In the semi-infinite nuclear matter approximation, the 
    meson field equations are
    \f{eqnarray*}
    \frac{d^2 \sigma}{d z^2} &=& 
    m_{\sigma}^2 \sigma - g_{\sigma} \left( n_{s n} + n_{s p} \right)
    + b M g_{\sigma}^3 \sigma^2 + c g_{\sigma}^4 \sigma^3 -
    g_{\rho}^2 \rho^2 \frac{\partial f}{\partial \sigma} \nonumber \\
    \frac{d^2 \omega}{d z^2} &=&
    m_{\omega}^2 \omega - g_{\omega} \left(n_n+n_p\right)
    + \frac{\zeta}{6} g_{\omega}^4 \omega^3 + g_{\rho}^2 \rho^2 
    \frac{\partial f}{\partial \omega} \nonumber \\
    \frac{d^2 \rho}{d z^2} &=&
    m_{\rho}^2 \rho + \frac{1}{2} g_{\rho} \left(n_n-n_p\right)
    + 2 g_{\rho}^2 \rho f + \frac{\xi}{6} g_{\rho}^4 \rho^3
    \f}
    in the same notation as \ref o2scl::eos_had_rmf .

    The algorithm works by solving first over a small interval in
    coordinate space and slowly spreading the solution out until the
    derivatives of the meson fields at the left boundary drop below a
    specified tolerance.
    
    \note wd2int() is broken since fesym() tends to fail at low
    densities. For this reason, it wasn't used in \ref Steiner05ia .
*/
class seminf_rel {

protected:
  
  /** \brief The grid size
   */
  int ngrid;
  /** \brief The saturation density at the current proton
      fraction
  */
  double nsat;
  /** \brief The current proton fraction
   */
  double protfrac;
  /** \brief The current neutron chemical potential
   */
  double mun;
  /** \brief The current proton chemical potential
   */
  double mup;
  /// \name The meson fields at the LHS
  //@{
  double sigma_left;
  double omega_left;
  double rho_left;
  //@}
  /** \brief The hadronic EOS
   */
  eos_had_rmf rmf_eos;
  /** \brief Thermodynamic variables
   */
  thermo hb;
  /** \brief Thermodynamic variables
   */
  thermo tht;
  /** \brief Neutron
   */
  fermion neutron;
  /** \brief Proton
   */
  fermion proton;
  
  /** \brief Equations for the boundary conditions in 
      the case of a neutron drip
   */
  int ndripfun(size_t sn, const si_vector_t &sx, si_vector_t &sy) {
    double pleft, pright, munleft, munright;

    neutron.n=sx[0];
    proton.n=sx[1];
    rmf_eos.calc_e(neutron,proton,hb);

    // Ensure the proton fraction on the LHS matches
    // the value in protfrac
    sy[0]=proton.n-protfrac*(proton.n+neutron.n);
    pleft=hb.pr;
    munleft=neutron.mu;

    neutron.n=sx[2];
    proton.n=0.001;
    rmf_eos.calc_e(neutron,proton,hb);

    pright=hb.pr;
    munright=neutron.mu;

    // Ensure the pressures and neutron chemical potentials
    // on the LHS and RHS match
    sy[1]=pleft-pright;
    sy[2]=munleft-munright;

    return 0;
  }

  /** \brief Differential equations to solve
   */
  double difeq(size_t ieq, double x, si_matrix_row_t &y) {
    
    neutron.mu=mun;
    proton.mu=mup;

    double fn, fn2, fn3;
    rmf_eos.calc_eq_p(neutron,proton,y[0],y[1],y[2],fn,fn2,fn3,hb);

    if (ieq==0) return y[3];
    else if (ieq==1) return y[4];
    else if (ieq==2) return y[5];
    else if (ieq==3) return fn;
    else if (ieq==4) return fn2;
    return fn3;
    
  }

  /** \brief LHS boundary conditions
   */
  double left(size_t ieq, double x, si_matrix_row_t &y) {
    if (ieq==0) return y[0]-sigma_left;
    if (ieq==1) return y[1]-omega_left;
    return y[2]-rho_left;
  }

  /** \brief RHS boundary conditions
   */
  double right(size_t ieq, double x, si_matrix_row_t &y) {
    if (ieq==0) return y[0];
    if (ieq==1) return y[1];
    return y[2];
  }
  
public:

  seminf_rel() {
    ngrid=100;
  }
  
  /** \brief Main
   */
  int run(int argc, char *argv[]) {
    
    bool convergeflag;
    bool flattendone;
    bool summaryout=true;
    bool iterfile=false;
    bool outputiter=true;
    bool debug=true;
    double conve=1.0e-10;
    double protfrac=0.0;
    double fact=1.04;
    double kfn;
    double kfp;
    double finalconverge=1.0e-12;
    double xcent;
    double derivlimit=0.04;
    int flattenit=70;
    static const int ne=6, nb=3;
    int lastit;
    int npoints=64;
    si_vector_t xstor(ngrid);
    si_matrix_t ystor(ne,ngrid);

    // Output quantities
    double wd=0.0;
    double wd2=0.0;
    double w0=0.0;
    double wdjl=0.0;
    double w0jl=0.0;
    double w02=0.0;
    double sssv1=0.0;
    double sssv2=0.0;
    double thick=0.0;
    double surf;
    double surf2;
    double sbulk;
    double sgrad;
    double xn[3];
    double xp[3];
    
    double temp;
    double nint;
    double dx=7.0/((double)ngrid);
    double xinterp;
    int ilast=0;
    int interp;

    //--------------------------------------------
    // Euqation solver

    mroot_hybrids<mm_funct,si_vector_t,si_matrix_t,jac_funct> nd;
    nd.tol_abs=1.0e-15;
    nd.tol_rel=1.0e-12;
    nd.ntrial=100;

    //--------------------------------------------
    // Equation of state and particle initializations

    rmf_load(rmf_eos,"RAPR");
    neutron.init(o2scl_settings.get_convert_units().convert
		 ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    proton.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_proton),2.0);
    neutron.non_interacting=false;
    proton.non_interacting=false;
    
    //--------------------------------------------
    // Initializations for ODE solver
    
    ode_it_solve<ode_it_funct,si_vector_t,si_matrix_t,si_matrix_row_t,
		 si_vector_t,si_matrix_t> oit;
    ode_it_funct f_derivs=std::bind
      (std::mem_fn<double(size_t,double,si_matrix_row_t &)>
       (&seminf_rel::difeq),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);       
    ode_it_funct f_left=std::bind
      (std::mem_fn<double(size_t,double,si_matrix_row_t &)>
       (&seminf_rel::left),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);       
    ode_it_funct f_right=std::bind
      (std::mem_fn<double(size_t,double,si_matrix_row_t &)>
       (&seminf_rel::right),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);   

#ifdef USE_EIGEN
    o2scl_linalg::linear_solver_eigen_colQR
      <Eigen::VectorXd,Eigen::MatrixXd> sol;
    oit.set_solver(sol);
#endif

    si_vector_t ox(ngrid);
    si_matrix_t oy(ngrid,ne);
    si_matrix_t A(ngrid*ne,ngrid*ne);
    si_vector_t rhs(ngrid*ne), dy(ngrid*ne);
    if (outputiter) oit.verbose=1;
    else oit.verbose=0;
    
    //--------------------------------------------
    // Create output table
    
    table<> at(ngrid);
    at.line_of_names("x sigma omega rho sigmap omegap rhop ");
    at.line_of_names(((string)"nn np n alpha nprime esurf esurf2 ebulk ")+
		     "egrad thickint wdint wd2int qpq ");
  
    //--------------------------------------------
    // Main loop over requested proton fractions

    for(int pf_index=1;pf_index<=2;pf_index++) {
      
      if (pf_index==1) {
	protfrac=0.5;
      } else {
	protfrac=0.49;
      }

      //-----------------------------------------------------------
      // Now calculate for saturation nuclear matter with proton
      // fraction given in parameter file
      
      double eoatemp;
      rmf_eos.set_n_and_p(neutron,proton);
      rmf_eos.set_thermo(tht);
      
      // Compute the saturation density for this proton fraction
      nsat=rmf_eos.fn0(1.0-2.0*protfrac,eoatemp);
      cout << "Saturation density at x=" << protfrac << ": "
	   << nsat << endl;

      // Compute the meson fields and chemical potentials at this
      // density
      double nn_left=neutron.n;
      double np_left=proton.n;
      neutron.n=nsat*(1.0-protfrac);
      proton.n=nsat*protfrac;
      rmf_eos.calc_e(neutron,proton,tht);
      rmf_eos.get_fields(sigma_left,omega_left,rho_left);
      mun=neutron.mu;
      mup=proton.mu;
      
      // The nuclear radius parameter at this saturation density
      double rad0=cbrt(0.75/pi/nsat);

      // Neutron drip case. Not working.
      if (mun>neutron.m) {

	protfrac=0.3;

	si_vector_t px(3);
	px[0]=0.07;
	px[1]=0.07;
	px[2]=0.04;
	
	mm_funct mff=std::bind
	  (std::mem_fn<int(size_t,const si_vector_t &,si_vector_t &)>
	   (&seminf_rel::ndripfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	cout << nd.msolve(3,px,mff) << endl;

	cout << px[0] << " " << px[1] << " " << px[2] << endl;
	cout << "Neutron drip case not working." << endl;
	exit(-1);
      }

      if (debug) {
	cout << "sigma_left: " << sigma_left
	     << "  omega_left: " << omega_left
	     << " rho_left: " << rho_left << endl;
	cout << " mun: " << mun << " mup: " << mup << endl;
	cout << endl;
      }

      if (true) {
      
	//----------------------------------------------
	// Construct LHS for initial guess
      
	ox[0]=0.0;
	oy(0,0)=sigma_left;
	oy(0,1)=omega_left;
	oy(0,2)=rho_left;
	oy(0,3)=-1.0e-3;
	oy(0,4)=-1.0e-3;
	if (fabs(protfrac-0.5)<0.0001) oy(0,5)=0.0;
	else oy(0,5)=1.0e-4;
	
	//----------------------------------------------
	// Integrate to determine initial guess

	neutron.mu=mun;
	proton.mu=mup;
	neutron.n=0.08;
	proton.n=0.08;

	bool guessdone=false;
	for(int i=1;i<ngrid;i++) {
	  ox[i]=ox[i-1]+dx;
	
	  if (guessdone==false) {

	    double f1, f2, f3;
	    rmf_eos.calc_eq_p(neutron,proton,oy(i-1,0),oy(i-1,1),
			      oy(i-1,2),f1,f2,f3,hb);
      
	    oy(i,0)=oy(i-1,0)+dx*oy(i-1,3);
	    oy(i,1)=oy(i-1,1)+dx*oy(i-1,4);
	    oy(i,2)=oy(i-1,2)+dx*oy(i-1,5);
	    oy(i,3)=oy(i-1,3)+dx*f1;
	    oy(i,4)=oy(i-1,4)+dx*f2;
	    oy(i,5)=oy(i-1,5)+dx*f3;

	    if (oy(i,0)<0.0 || oy(i,1)<0.0) {
	      guessdone=true;
	      ilast=i;
	      i=ngrid+10;
	    }
	  } else {
	    for(int j=0;j<6;j++) oy(i,j)=0.0;
	  }
	}

	//--------------------------------------------------------------
	// Stretch the solution over the entire grid. Start at the RHS,
	// so that we can overwrite the original data. This works until
	// we get to the left hand side, where we adjust accordingly
	
	ilast=18;
	for(int i=ngrid-1;i>=0;i--) {
	  interp=(int)(((double)(i+1))/((double)ngrid)*((double)(ilast+1)))-1;
	  xinterp=dx*((double)(i+1))/((double)ngrid)*((double)(ilast+1));
	  ox[i]=xinterp;
	  if (i>9) {
	    for(int j=0;j<6;j++) {
	      oy(i,j)=oy(interp,j)+(xinterp-ox[interp])/dx*
		(oy(interp+1,j)-oy(interp,j));
	    }
	  } else {
	    oy(i,0)=sigma_left;
	    oy(i,1)=omega_left;
	    oy(i,2)=rho_left;
	    for(int j=3;j<6;j++) oy(i,j)=0.0;
	  }
	}

	//----------------------------------------------
	// Now we center the x-axis on zero, which makes it
	// easier to expand the grid later

	xcent=ox[ngrid/2-1];
	for(int i=0;i<ngrid;i++) {
	  ox[i]-=xcent;
	}

	// End of calculation of initial guess
      } else {

	// Load initial guess from file
	hdf_file hf;
	hf.open("rel.o2");
	string tname;
	if (pf_index==1) tname="rel1";
	else tname="rel2";
	hdf_input(hf,at,tname);
	hf.close();
	for(int i=0;i<ngrid;i++) {
	  ox[i]=at.get("x",i);
	  oy(i,0)=at.get("sigma",i);
	  oy(i,1)=at.get("omega",i);
	  oy(i,2)=at.get("rho",i);
	  oy(i,3)=at.get("sigmap",i);
	  oy(i,4)=at.get("omegap",i);
	  oy(i,5)=at.get("rhop",i);
	}

      }

      convergeflag=true;
      flattendone=false;

      int j;
      for(j=1;j<=flattenit && flattendone==false && 
	    convergeflag==true;j++) {
      
	//--------------------------------------------
	// Store most recent result
	for(int i=0;i<ngrid;i++) {
	  xstor[i]=ox[i];
	  for(int k=0;k<ne;k++) {
	    ystor(k,i)=oy(i,k);
	  }
	}
	
	//----------------------------------------------
	// Try ode_it_solve
	
	oit.tol_rel=conve;
	oit.solve(ngrid,ne,nb,ox,oy,f_derivs,f_left,f_right,
		  A,rhs,dy);
	  
	if (outputiter) {
	  cout << "j: ";
	  cout.width(2);
	  cout << j << "  x(left): " << ox(ngrid-1)
	       << " sigma'(left): " << oy(0,3) << endl;
	  cout << " omega'(left): " << oy(0,4);
	  cout.setf(ios::showpos);
	  cout << "   rho'(left): " << oy(0,5);
	  cout.unsetf(ios::showpos);
	  cout << " conflag: " << convergeflag << endl;
	}
      
	if (fabs(oy(0,3))<derivlimit && fabs(oy(0,4))<derivlimit &&
	    fabs(oy(0,5))<derivlimit) {
	  flattendone=true;
	  if (outputiter) cout << endl;
	} else {
	  // Why does this work? The alternative of simply extending
	  // the LHS doesn't seem to work.
	  if (convergeflag==true) {
	    for(int i=0;i<ngrid;i++) {
	      ox[i]*=fact;
	    }
	  }
	}
	
      }
      lastit=j-1;

      if (flattendone==false) {
	cout << "Could not get vanishing derivative at boundary." << endl;
      } else if (convergeflag==false) {
	cout << "Relaxation did not converge." << endl;

	//--------------------------------------------
	// If convergence failed, use previous result
	// for calculations

	for(int i=0;i<ngrid;i++) {
	  ox[i]=xstor[i];
	  for(int k=0;k<ne;k++) {
	    oy(i,k)=ystor(k,i);
	  }
	}
    
      } else {
	cout << "Going to final solution." << endl;
	oit.tol_rel=finalconverge;
	oit.solve(ngrid,ne,nb,ox,oy,f_derivs,f_left,f_right,
		  A,rhs,dy);
      }

      //--------------------------------------------
      // Store final solution to table

      for(int i=0;i<ngrid;i++) {
	at.set("x",i,ox[i]);
	at.set("sigma",i,oy(i,0));
	at.set("omega",i,oy(i,1));
	at.set("rho",i,oy(i,2));
	at.set("sigmap",i,oy(i,3));
	at.set("omegap",i,oy(i,4));
	at.set("rhop",i,oy(i,5));
      }
      
      //--------------------------------------------
      // Calculate rhon, rhop, energy, and ebulk at
      // every point in profile
      
      double delta=1.0-2.0*protfrac;
      double hns=rmf_eos.fesym(nsat,delta);

      neutron.mu=mun;
      proton.mu=mup;

      for(int i=0;i<ngrid;i++) {

	double f1, f2, f3;
	rmf_eos.calc_eq_p(neutron,proton,oy(i,0),oy(i,1),
			  oy(i,2),f1,f2,f3,hb);

	at.set("nn",i,neutron.n);
	at.set("np",i,proton.n);
	at.set("n",i,at.get("nn",i)+at.get("np",i));
	at.set("alpha",i,at.get("nn",i)-at.get("np",i));
	
	at.set("ebulk",i,hb.ed-mup*at.get("np",i)-mun*at.get("nn",i));
	at.set("egrad",i,0.5*(oy(i,3)*oy(i,3)-oy(i,4)*oy(i,4)-
			      oy(i,5)*oy(i,5)));
	
	at.set("esurf",i,at.get("ebulk",i)+at.get("egrad",i));
	if (at.get("n",i)>0.0) {
	  at.set("esurf2",i,at.get("esurf",i));
	} else {
	  at.set("esurf2",i,0.0);
	}
      }

      at.delete_column("nprime");
      at.deriv("x","n","nprime");

      for(int i=0;i<ngrid;i++) {
	
	if (at.get("nprime",i)!=0.0) {

	  // Here qpq=Q_nn+Q_np is calculated from the bulk energy.
	  // One could calculate this from the gradient part of the
	  // surface energy as well. There is not much difference.
	  
	  at.set("qpq",i,(hb.ed-mup*at.get("np",i)-
			  mun*at.get("nn",i))*4.0/at.get("nprime",i)/
		 at.get("nprime",i));
	  
	  //qpq[i]=(oy(i,3)*oy(i,3)-
	  // oy(i,4)*oy(i,4)-
	  // oy(i,5)*oy(i,5))*2.0/at.get("nprime",i)/
	  // at.get("nprime",i);
	  
	} else {
	  at.set("qpq",i,0.0);
	}
	
	at.set("thickint",i,(at.get("nn",i)/nn_left-at.get("np",i)/np_left));
	
	if (pf_index>=2) {
	  at.set("wdint",i,at.get("alpha",i)/delta-at.get("n",i));
	  if (at.get("n",i)>0.0) {
	    at.set("wd2int",i,at.get("n",i)*
		   (pow(at.get("alpha",i)/delta/at.get("n",i),2.0)*
		    rmf_eos.fesym(at.get("n",i))-hns));
	  } else {
	    at.set("wd2int",i,0.0);
	  }
	}
      }
    
      //----------------------------------------------------------
      // Calculate the value of x for nn*0.9, nn*0.5, nn*0.1, etc.
      
      o2scl::interp<std::vector<double>,si_vector_t> it(itp_linear);
      
      for(int i=0;i<3;i++) {
	xn[i]=it.eval(nn_left/10.0*((double)(i*4+1)),ngrid,
		      at.get_column("nn"),ox);
	xp[i]=it.eval(np_left/10.0*((double)(i*4+1)),ngrid,
		      at.get_column("np"),ox);
      }

      //----------------------------------------------------------
      // Integrals
      
      o2scl::interp<std::vector<double>,std::vector<double> > it2;

      const std::vector<double> &x_vec=at.get_column("x");
      double x_left=x_vec[0];
      double x_right=x_vec[at.get_nlines()-1];
      size_t istt=at.get_nlines()-1;

      surf=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		     at.get_column("esurf"));
      surf2=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		      at.get_column("esurf2"));
      sbulk=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		      at.get_column("ebulk"));
      sgrad=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		      at.get_column("egrad"));
      thick=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		      at.get_column("thickint"));
      
      if (pf_index>=2) {
	
	wd=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		     at.get_column("wdint"));
	wd*=hns;
	
	wd2=it2.integ(x_left,x_right,at.get_nlines(),x_vec,
		      at.get_column("wd2int"));
	
	sssv1=4.0*pi*rad0*rad0*wd/hns;
	
	cout << "wd, wd2, sssv1: " << wd << " " << wd2 << " " << sssv1 << endl;
      } else {
	w0=surf;
	w02=surf2;
	cout << "w0, w02: " << w0 << " " << w02 << endl;
      }

#ifdef NEVER_DEFINED
      
      //--------------------------------------------
      // Jim's formulas for surface tension, etc.

      if (pf_index==1) {
	w0jl=0.0;
	wdjl=0.0;

	inte_qag_gsl gl;
	double xint=0.5*nsat, tweight;
      
	for(int i=0;i<npoints*2;i++) {
	  if (i>npoints) {
	    nint=xint+xint*gl->get_abscissa(i-npoints);
	    tweight=xint*gl->get_weight(i-npoints);
	  } else {
	    nint=xint-xint*gl->get_abscissa(i);
	    tweight=xint*gl->get_weight(i);
	  }

	  si_matrix_row_t ar1(*rel->y,1);
	  si_matrix_row_t ar2(*rel->y,2);
	  si_matrix_row_t ar3(*rel->y,3);
	
	  neutron.n=lookup(ngrid,nint,rhon,rho);
	  proton.n=lookup(ngrid,nint,rhop,rho);
	  double sig=lookup(ngrid,nint,ar1,rho);
	  double ome=lookup(ngrid,nint,ar2,rho);
	  double rhof=lookup(ngrid,nint,ar3,rho);
	  double qq=lookup(ngrid,nint,qpq,rho);

	  neutron.kffromden();
	  proton.kffromden();
	  neutron.ms=neutron.m-gs*sig;
	  proton.ms=neutron.ms;
	  neutron.nu=sqrt(neutron.kf*neutron.kf+neutron.ms*neutron.ms);
	  proton.nu=sqrt(proton.kf*proton.kf+proton.ms*proton.ms);
	
	  double estmp;
	
	  rmf_eos.calc_p(neutron,proton,sig,ome,rhof,f1,f2,f3,hb);
	  //hb.ed=-hb.pr+neutron.n*neutron.mu+proton.n*proton.mu;
	
	  estmp=rmf_eos.fesym(nint);
	  // fesym() automatically computes the bulk energy density
	  // and puts the result into tht.ed:
	  hb.ed=tht.ed;

	  if (hb.ed-mun*neutron.n-mup*proton.n>0.0 && qq>0.0) {
	    wdjl+=tweight*sqrt(qq)*nint*(hns/estmp-1.0)/
	      sqrt(hb.ed-mun*neutron.n-mup*proton.n);
	    w0jl+=tweight*sqrt(fabs(qq*(hb.ed-mun*neutron.n-
					mup*proton.n)));
	    if (!finite(wdjl)) {
	      cout << "wdjl not finite." << endl;
	      cout << hb.ed-mun*neutron.n-mup*proton.n << endl;
	      cout << hns << " " << qq << endl;
	      exit(-1);
	    }
	  }
	}
      
	//delete gl;
	wdjl*=hns/2.0;
	sssv2=4*pi*rad0*rad0*wdjl/hns;
      }
      
#endif
      
      if (true) {
	hdf_file hf;
	string tablename=((string)"rel")+std::to_string(pf_index);
	hf.open_or_create("rel.o2");
	hdf_output(hf,at,tablename);
	hf.close();
	cout << "Wrote solution to file 'rel.o2'" << endl;
	cout << endl;
      }
      
    }

    return 0;
  }

};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  seminf_rel sn;
  sn.run(argc,argv);

  if (argc>=2 && ((std::string)argv[1])=="check") {
    test_mgr t;
    t.set_output_level(2);
    hdf_file hf;
    string name;
    table_units<> tab, tab_expected;

    name="rel1";

    hf.open("rel.o2");
    hdf_input(hf,tab,name);
    hf.close();

    hf.open("rel_save.o2");
    hdf_input(hf,tab_expected,name);
    hf.close();

    t.test_rel_nonzero_table(tab,tab_expected,1.0e-12,1.0e-8,"table 1");

    name="rel2";

    hf.open("rel.o2");
    hdf_input(hf,tab,name);
    hf.close();

    hf.open("rel_save.o2");
    hdf_input(hf,tab_expected,name);
    hf.close();

    t.test_rel_nonzero_table(tab,tab_expected,1.0e-12,1.0e-8,"table 2");
    
    if (!t.report()) {
      exit(-1);
    }
  }
  
  return 0;
}
