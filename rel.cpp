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
/*
  \note Initial guess reading might be broken
  
  \note wd2int() is broken since fesym() tends to 
  fail at low densities. This is ok, however, since wd(drop2) is
  not used in any of the figures.
*/

#include <cmath>
#include <iostream>
#include <functional>

#include <o2scl/table.h>
#include <o2scl/constants.h>
#include <o2scl/part.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/lib_settings.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mm_funct.h>
#include <o2scl/interp.h>
#include <o2scl/ode_it_solve.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

/** \brief Desc
 */
class ode_it_solve2 :
  public ode_it_solve<ode_it_funct11,ubvector,ubmatrix,ubmatrix_row,
		      ubvector,ubmatrix> {

protected:
  
  /** \brief Desc
   */
  double eps;
  
public:
  
  /** \brief Desc
   */
  ode_it_solve2() : ode_it_solve<ode_it_funct11,ubvector,ubmatrix,
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

/** \brief Desc
 */
class seminf_rel {

protected:
  
  /** \brief The grid size
   */
  int ngrid;
  double nsat;
  double protfrac;
  /** \brief The current neutron chemical potential
   */
  double mun;
  /** \brief The current proton chemical potential
   */
  double mup;
  double phi0;
  double v0;
  double r0;
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
  
  /** \brief Desc
   */
  int nucmat(size_t nv, const ubvector &ex, ubvector &ey) {
    double f1, f2, f3, sig, ome, rho;

    neutron.nu=ex[0];
    proton.nu=ex[1];
    sig=ex[2];
    ome=ex[3];
    rho=ex[4];
    
    rmf_eos.calc_eq_p(neutron,proton,sig,ome,rho,f1,f2,f3,hb);

    ey[0]=proton.n+neutron.n-nsat;
    ey[1]=proton.n-nsat*protfrac;
    ey[2]=f1;
    ey[3]=f2;
    ey[4]=f3;

    return 0;
  }

  /** \brief Desc
   */
  int ndripfun(size_t sn, const ubvector &sx, ubvector &sy) {
    double pleft, pright, munleft, munright;
  
    neutron.n=sx[0];
    proton.n=sx[1];
    rmf_eos.calc_e(neutron,proton,hb);
  
    sy[0]=proton.n-protfrac*(proton.n+neutron.n);
    pleft=hb.pr;
    munleft=neutron.mu;

    neutron.n=sx[2];
    proton.n=0.001;
    rmf_eos.calc_e(neutron,proton,hb);

    pright=hb.pr;
    munright=neutron.mu;

    sy[1]=pleft-pright;
    sy[2]=munleft-munright;

    return 0;
  }

  /** \brief Differential equations to solve
   */
  double difeq2(size_t ieq, double x, ubmatrix_row &y) {

    double gw=rmf_eos.cw*rmf_eos.mw;
    double gr=rmf_eos.cr*rmf_eos.mr;

    double phi=y[0];
    double v=y[1];
    double r=y[2];
    double phip=y[3];
    double vp=y[4];
    double rprime=y[5];
    
    neutron.nu=mun-gw*v+0.5*gr*r;
    proton.nu=mup-gw*v-0.5*gr*r;

    double fn, fn2, fn3;
    rmf_eos.calc_eq_p(neutron,proton,phi,v,r,fn,fn2,fn3,hb);

    if (ieq==0) {
      return phip;
    } else if (ieq==1) {
      return vp;
    } else if (ieq==2) {
      return rprime;
    } else if (ieq==3) {
      return fn;
    } else if (ieq==4) {
      return fn2;
    }
    return fn3;
    
  }

  /** \brief LHS boundary conditions
   */
  double left2(size_t ieq, double x, ubmatrix_row &y) {
    if (ieq==0) {
      return y[0]-phi0;
    }
    if (ieq==1) {
      return y[1]-v0;
    }
    return y[2]-r0;
  }

  /** \brief RHS boundary conditions
   */
  double right2(size_t ieq, double x, ubmatrix_row &y) {
    if (ieq==0) return y[0];
    if (ieq==1) return y[1];
    return y[2];
  }
  
  /** \brief Desc
   */
  template<class vec_t, class vec2_t> 
  double lookup(int n, double yy0, const vec_t &x,
		const vec2_t &y) {

    for(int i=2;i<=n;i++) {
      if ((y[i]>=yy0 && y[i-1]<yy0) || (y[i]<yy0 && y[i-1]>=yy0)) {
	double x0=x[i-1]+(x[i]-x[i-1])*(yy0-y[i-1])/(y[i]-y[i-1]);
	return x0;
      }
    }

    return 0.0;
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
    double xn[4];
    double xp[4];
    double fact=1.04;
    double gs;
    double gw;
    double gr;
    double mstar;
    double kfn;
    double kfp;
    double finalconverge=1.0e-12;
    double xcent;
    double derivlimit=0.04;
    int flattenit=70;
    static const int ne=6, nb=3;
    int lastit;
    int npoints=64;
    ubvector xstor(ngrid);
    ubmatrix ystor(ne,ngrid);

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
    
    double temp;
    double nint;
    double f1;
    double f2;
    double f3;
    double dx=7.0/((double)ngrid);
    double xinterp;
    int ilast=0;
    int interp;

    //--------------------------------------------
    // Euqation solver

    mroot_hybrids<> nd;
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
    // Initializations for ODE solver class
    
    ode_it_solve2 oit;
    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&seminf_rel::difeq2),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&seminf_rel::left2),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&seminf_rel::right2),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3);   
    ubvector ox(ngrid);
    ubmatrix oy(ngrid,ne);
    ubmatrix A(ngrid*ne,ngrid*ne);
    ubvector rhs(ngrid*ne), dy(ngrid*ne);
    if (outputiter) oit.verbose=1;
    else oit.verbose=0;
    
    //--------------------------------------------
    // Create output table

    table<> at(1000);
    at.line_of_names("x sigma omega rho sigmap omegap rhop ");
    at.line_of_names(((string)"nn np n alpha nprime esurf esurf2 ebulk ")+
		     "egrad thickint wdint wd2int qpq ");
    at.set_nlines(ngrid);
  
    //--------------------------------------------
    // Main loop over requested proton fractions

    for(int pf_index=1;pf_index<=2;pf_index++) {
      
      if (pf_index==1) {
	protfrac=0.5;
      } else {
	protfrac=0.49;
      }

      //--------------------------------------------
      // Now calculate for saturation nuclear matter with proton
      // fraction given in parameter file 
  
      protfrac=protfrac;
      double eoatemp;
      rmf_eos.set_n_and_p(neutron,proton);
      rmf_eos.set_thermo(tht);
    
      nsat=rmf_eos.fn0(1.0-2.0*protfrac,eoatemp);
      rmf_eos.get_fields(phi0,v0,r0);
    
      cout << "Saturation density at x=" << protfrac << ": "
	   << nsat << endl;

      double rad0=0.0;
      if (pf_index==1) rad0=cbrt(0.75/pi/nsat);
    
      mun=neutron.mu;
      mup=proton.mu;
      double nn_left=neutron.n;
      double np_left=proton.n;

      // Neutron drip case. Not working.
      if (mun>neutron.m) {

	protfrac=0.3;

	ubvector px(3);
	px[0]=0.07;
	px[1]=0.07;
	px[2]=0.04;
	
	mm_funct11 mff=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	   (&seminf_rel::ndripfun),this,std::placeholders::_1,
	   std::placeholders::_2,std::placeholders::_3);
	cout << nd.msolve(3,px,mff) << endl;

	cout << px[0] << " " << px[1] << " " << px[2] << endl;
	cout << "Neutron drip case not working." << endl;
	exit(-1);
      }

      if (debug) {
	cout << "phi0: " << phi0 << "  v0: " << v0 << " r0: "
	     << r0 << endl;
	cout << " mun: " << mun << " mup: " << mup << endl;
	cout << endl;
      }

      gs=rmf_eos.cs*rmf_eos.ms;
      gw=rmf_eos.cw*rmf_eos.mw;
      gr=rmf_eos.cr*rmf_eos.mr;
  
      //----------------------------------------------
      // Construct LHS for initial guess
      
      ox[0]=0.0;
      oy(0,0)=phi0;
      oy(0,1)=v0;
      oy(0,2)=r0;
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
	  neutron.nu=mun-gw*oy(i-1,1)+0.5*gr*oy(i-1,2);
	  proton.nu=mup-gw*oy(i-1,1)-0.5*gr*oy(i-1,2);
	    
	  neutron.mu=mun;
	  proton.mu=mup;
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
	  oy(i,0)=phi0;
	  oy(i,1)=v0;
	  oy(i,2)=r0;
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
      // Output characteristics of final solution

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
    
      double hns=rmf_eos.fesym(nsat,protfrac);
      double delta=1.0-2.0*protfrac;
    
      for(int i=0;i<ngrid;i++) {
	mstar=rmf_eos.mnuc-gs*oy(i,0);
    
	if (mun-gw*oy(i,1)+gr*oy(i,2)/2.0<mstar) {
	  kfn=0.0;
	} else {
	  kfn=sqrt(pow(mun-gw*oy(i,1)+gr*oy(i,2)/2.0,2.0)-
		   mstar*mstar);
	}
	if (mup-gw*oy(i,1)-gr*oy(i,2)/2.0<mstar) {
	  kfp=0.0;
	} else {
	  kfp=sqrt(pow(mup-gw*oy(i,1)-gr*oy(i,2)/2.0,2.0)-
		   mstar*mstar);
	}

	at.set("nn",i,pow(kfn,3.0)/3.0/pi2);
	at.set("np",i,pow(kfp,3.0)/3.0/pi2);
	at.set("n",i,at.get("nn",i)+at.get("np",i));
	at.set("alpha",i,at.get("nn",i)-at.get("np",i));

	neutron.ms=mstar;
	proton.ms=mstar;
	neutron.nu=sqrt(kfn*kfn+neutron.ms*neutron.ms);
	proton.nu=sqrt(kfp*kfp+proton.ms*proton.ms);
	rmf_eos.calc_eq_p(neutron,proton,oy(i,0),oy(i,1),
			  oy(i,2),f1,f2,f3,hb);
	hb.ed=-hb.pr+neutron.n*neutron.mu+proton.n*proton.mu;

	at.set("ebulk",i,hb.ed-mup*at.get("np",i)-
	       mun*at.get("nn",i));
	at.set("egrad",i,0.5*(oy(i,3)*oy(i,3)-oy(i,4)*oy(i,4)-
			      oy(i,5)*oy(i,5)));

	at.set("esurf",i,at.get("ebulk",i)+at.get("egrad",i));
	if (at.get("n",i)>0.0) {
	  at.set("esurf2",i,at.get("esurf",i));
	} else {
	  at.set("esurf2",i,0.0);
	}
      
	// There doesn't seem to be much difference between this approach
	// and the one below
	at.set("nprime",i,(kfn*(neutron.nu*(-gw*oy(i,4)+gr/2.0*oy(i,5))+
				mstar*gs*oy(i,3))+
			   kfp*(proton.nu*(-gw*oy(i,4)-gr/2.0*oy(i,5))+
				mstar*gs*oy(i,3)))/pi2);
	//if (i==0) rhoprime[i]=0;
	//else rhoprime[i]=(rho[i]-rho[i-1])/(rel->x[i]-rel->x[i-1]);
	
	if (at.get("nprime",i)!=0.0) {
	  // Here qpq=Q_nn+Q_np is calculated from the bulk energy. One
	  // Could calculate this from the gradient part of the surface
	  // energy as well. There is not much difference.
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
	  if (at.get("n",i)<0.0) {
	    at.set("wd2int",i,at.get("n",i)*
		   (pow(at.get("alpha",i)/delta/at.get("n",i),2.0)*
		    rmf_eos.fesym(at.get("n",i))-hns));
	  } else {
	    at.set("wd2int",i,0.0);
	  }
	}
      }
    
      //--------------------------------------------
      // Calculate the value of x for nn*0.9, nn*0.5, 
      // nn*0.1, etc.

      mstar=rmf_eos.mnuc-gs*phi0;
      kfn=sqrt(pow(mun-gw*v0+gr*r0/2.0,2.0)-mstar*mstar);
      kfp=sqrt(pow(mup-gw*v0-gr*r0/2.0,2.0)-mstar*mstar);

      for(int i=1;i<=3;i++) {
	xn[i]=lookup(ngrid,pow(kfn,3.0)/3.0/pi2/10.0*((double)(i*4-3)),
		     ox,at.get_column("nn"));
	xp[i]=lookup(ngrid,pow(kfp,3.0)/3.0/pi2/10.0*((double)(i*4-3)),
		     ox,at.get_column("np"));
      }

      cout << "Integrals." << endl;
      o2scl::interp<std::vector<double>,std::vector<double> > gi;

      const std::vector<double> &xav=at[0];
      size_t istt=at.get_nlines()-1;

      surf=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
		    at.get_column("esurf"));
      surf2=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
		     at.get_column("esurf2"));
      sbulk=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
		     at.get_column("ebulk"));
      sgrad=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
		     at.get_column("egrad"));
      thick=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
		     at.get_column("thickint"));

      if (pf_index>=2) {

	wd=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
		    at.get_column("wdint"));
	wd*=hns;

	wd2=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column("x"),
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

	  ubmatrix_row ar1(*rel->y,1);
	  ubmatrix_row ar2(*rel->y,2);
	  ubmatrix_row ar3(*rel->y,3);
	
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
      
      hdf_file hf;
      string tablename=((string)"rel")+std::to_string(pf_index);
      hf.open_or_create("rel.o2");
      hdf_output(hf,at,tablename);
      hf.close();
      cout << "Wrote solution to file 'rel.o2'" << endl;
      cout << endl;
      
    }

    return 0;
  }

};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  seminf_rel sn;
  sn.run(argc,argv);
  return 0;
}
