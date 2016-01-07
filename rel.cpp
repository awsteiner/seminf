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
  1/03 - Added option of using proton fraction from lead nucleus

  11/03 - Initial guess reading broken

  2/20/04 - wd2int() is broken since fesym() tends to 
  fail at low densities. This is ok, however, since wd(drop2) is
  not used in any of the figures.

  5/5/05 - Converted to new o2scl vectors, etc.
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

#include "relax.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

static const int ngrid=700;

//--------------------------------------------
// Class definition for relaxation routine

class s3relax : public relax {
public:
  s3relax(int tne, int tnb, int tngrid);
  int iter(int k, double err, double fac, ubvector_int &kmax,
	   ubvector &errmax);
  int difeq(int k, int k1, int k2, int jsf, int is1, int isf);
};

s3relax::s3relax(int tne, int tnb, int tngrid) : 
  relax(tne,tnb,tngrid) {}

//--------------------------------------------
// Structure for relax parameters

typedef struct relaxp_s {
  double mun;
  double mup;
  double maxdev;
  double phi0;
  double v0;
  double r0;
  bool showdiff;
  bool showrelax;
  bool convergeflag;
  bool debug;
  eos_had_rmf rmf_eos;
  fermion n;
  fermion p;
  thermo hb;
} relaxp;

class seminf_rel {
  
public:
  
  //--------------------------------------------
  // Structure for newt parameters

  typedef struct newtp_s {
    // barn and protfrac are always the baryon number
    // and the proton fraction at -infinity
    double nsat;
    double protfrac;
    eos_had_rmf rmf_eos;
    fermion n;
    fermion p;
    thermo hb;
  } newtp;

  thermo tht;

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

  //--------------------------------------------
  // Main

  int run(int argv, char *argc[]) {

    cout.setf(ios::scientific);

    bool flattendone;
    bool summaryout;
    bool iterfile;
    bool altsym;
    bool outputiter;
    bool debug=true;
    char ch;
    double conve;
    double slowc;
    double jlint;
    double senuc;
    double surf;
    double protfrac=0.0;
    double lamv;
    double lam4;
    double xn[4];
    double xp[4];
    double fact;
    double gs;
    double gw;
    double gr;
    double mstar;
    double kfn;
    double kfp;
    double finalconverge;
    double xcent;
    double derivlimit;
    ubvector px(7);
    double wd=0.0;
    double wd2=0.0;
    double w0=0.0;
    ifstream fin; 
    int flattenit;
    int sz;
    int sz2;
    static const int ne=6, nb=3;
    int lastit, npoints=64;
    ofstream fout;
    double xstor[ngrid+1], ystor[ne+1][ngrid+1];
    string dirname;

    double sssv1=0.0;
    double sssv2=0.0;
    double temp;
    double temppf;
    double qq;
    double sig;
    double ome;
    double rhof;
    double hns;
    double wdjl=0.0;
    double nint;
    double w0jl=0.0;
    double w02=0.0;
    double nn1;
    double np1;
    double nn2;
    double np2;
    double hbarn;
    double *ebulklam;
    double nn0;
    double np0;
    double sbulk;
    double sgrad;
    double *rhoprime;
    double surf2;
    double f1;
    double f2;
    double f3;
    double thick=0.0;
    double n0half;
    double sprotfrac;
    double delta;
    double *qpq;
    bool guessdone=false;
    double dx=0.01;
    double fn;
    double fn2;
    double fn3;
    double xinterp;
    int ilast=0;
    int interp;
    bool usehc;
    string hcmodel;
    table<> at(1000);

    //--------------------------------------------
    // Equation of state initializations

    mroot_hybrids<> nd;
    nd.tol_abs=1.0e-15;
    nd.tol_rel=1.0e-12;
    nd.ntrial=100;
    mroot_hybrids<> nd2;
    nd2.tol_abs=1.0e-15;
    nd2.tol_rel=1.0e-12;
    nd2.ntrial=100;

    eos_had_rmf rmf_eos;
    newtp np;
    np.n.init(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    np.p.init(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_proton),2.0);
    np.n.non_interacting=false;
    np.p.non_interacting=false;

    //--------------------------------------------
    // Initializations for relaxation class

    s3relax *rel;

    rel=new s3relax(ne,nb,ngrid);
  
    at.line_of_names("x sigma omega rho sigmap omegap rhop ");
  
    slowc=1.0;
  
    relaxp rp;
    rp.debug=debug;
    rp.showdiff=false;
    rp.n=np.n;
    rp.p=np.p;
    rp.hb=np.hb;

    double r0=cbrt(0.75/pi/n0half);

    for(int pf_index=1;pf_index<=2;pf_index++) {

      if (pf_index==1) {
	protfrac=0.5;
      } else {
	protfrac=sprotfrac;
      }

      //--------------------------------------------
      // Now calculate for saturation nuclear matter with proton
      // fraction given in parameter file 
  
      np.protfrac=protfrac;
      double eoatemp;
      rmf_eos.set_n_and_p(np.n,np.p);
      rmf_eos.set_thermo(tht);
    
      np.nsat=rmf_eos.fn0(protfrac,eoatemp);
      rmf_eos.get_fields(rp.phi0,rp.v0,rp.r0);
    
      cout << "Saturation density at x=" << protfrac << ": " << np.nsat << endl;

      rp.mun=np.n.mu;
      rp.mup=np.p.mu;
      nn0=np.n.n;
      np0=np.p.n;

      if (rp.mun>np.n.m) {

	np.protfrac=0.3;

	px[0]=0.07;
	px[1]=0.07;
	px[2]=0.04;

	/*
	  mm_funct11 mff=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	  (&seminf_rel::ndripfun),this,std::placeholders::_1,
	  std::placeholders::_2,std::placeholders::_3);
	*/
	mm_funct11 mff;
	cout << nd.msolve(3,px,mff) << endl;

	cout << px[0] << " " << px[1] << " " << px[2] << endl;
	cout << "Exiting." << endl;
	exit(-1);
      }

      if (debug) {
	cout << "phi0 v0 r0 mun mup" << endl 
	     << rp.phi0 << " " << rp.v0 << " " << rp.r0 << " " << rp.mun
	     << " " << rp.mup << endl;
      }

      gs=rmf_eos.cs*rmf_eos.ms;
      gw=rmf_eos.cw*rmf_eos.mw;
      gr=rmf_eos.cr*rmf_eos.mr;
  
      //----------------------------------------------
      // Construct guess

      if (argv<3) {
	rel->x[1]=0.0;
	rel->y(1,1)=rp.phi0;
	rel->y(2,1)=rp.v0;
	rel->y(3,1)=rp.r0;
	rel->y(4,1)=-1.0e-3;
	rel->y(5,1)=-1.0e-3;
	if (fabs(protfrac-0.5)<0.0001) rel->y(6,1)=0.0;
	else rel->y(6,1)=1.0e-4;

	//----------------------------------------------
	// Use simple integration

	for(int i=2;i<=ngrid;i++) {
	  rel->x[i]=rel->x[i-1]+dx;

	  if (guessdone==false) {
	    rp.n.nu=rp.mun-gw*rel->y(2,i-1)+0.5*gr*rel->y(3,i-1);
	    rp.p.nu=rp.mup-gw*rel->y(2,i-1)-0.5*gr*rel->y(3,i-1);
	    
	    rmf_eos.calc_eq_p(rp.n,rp.p,rel->y(1,i-1),rel->y(2,i-1),
			      rel->y(3,i-1),fn,fn2,fn3,rp.hb);
      
	    rel->y(1,i)=rel->y(1,i-1)+dx*rel->y(4,i-1);
	    rel->y(2,i)=rel->y(2,i-1)+dx*rel->y(5,i-1);
	    rel->y(3,i)=rel->y(3,i-1)+dx*rel->y(6,i-1);
	    rel->y(4,i)=rel->y(4,i-1)+dx*fn;
	    rel->y(5,i)=rel->y(5,i-1)+dx*fn2;
	    rel->y(6,i)=rel->y(6,i-1)+dx*fn3;

	    if (rel->y(1,i)<0.0 || rel->y(2,i)<0.0) {
	      guessdone=true;
	      ilast=i-1;
	      i=ngrid+10;
	    }
	  } else {
	    for(int j=1;j<=6;j++) rel->y(j,i)=0.0;
	  }
	}

	//----------------------------------------------
	// Stretch the solution over the entire grid
	// Start at the RHS, so that we can overwrite
	// the original data. This works until we get
	// to the left hand side, where we adjust accordingly

	for(int i=ngrid;i>=1;i--) {
	  interp=(int)(((double)i)/((double)ngrid)*((double)ilast));
	  xinterp=dx*((double)i)/((double)ngrid)*((double)ilast);
	  rel->x[i]=xinterp;
	  if (i>10) {
	    for(int j=1;j<=6;j++) 
	      rel->y(j,i)=rel->y(j,interp)+
		(xinterp-rel->x[interp])/dx*
		(rel->y(j,interp+1)-rel->y(j,interp));
	  } else {
	    rel->y(1,i)=rp.phi0;
	    rel->y(2,i)=rp.v0;
	    rel->y(3,i)=rp.r0;
	    for(int j=4;j<=6;j++) rel->y(j,i)=0.0;
	  }
	}
      }
    
      //----------------------------------------------
      // Now we center the x-axis on zero, which makes it
      // easier to expand the grid later

      xcent=rel->x[ngrid/2];
      for(int i=1;i<=ngrid;i++) {
	rel->x[i]-=xcent;
      }

      rp.convergeflag=true;
      flattendone=false;

      int j;
      for(j=1;j<=flattenit && flattendone==false && 
	    rp.convergeflag==true;j++) {
      
	//--------------------------------------------
	// Store most recent result
	for(int i=1;i<=ngrid;i++) {
	  xstor[i]=rel->x[i];
	  for(int k=1;k<=ne;k++) {
	    ystor[k][i]=rel->y(k,i);
	  }
	}
  
	//--------------------------------------------
	// Solve equations
      
	rel->solve(conve,slowc);
  
	if (outputiter) {
	  cout << j << " " << rel->x[ngrid]-rel->x[1] << " " 
	       << rel->y(4,1) << " " << rel->y(5,1) << " " 
	       << rel->y(6,1) << " " 
	       << rp.convergeflag << endl;
	}
      
	if (fabs(rel->y(4,1))<derivlimit && 
	    fabs(rel->y(5,1))<derivlimit && 
	    fabs(rel->y(6,1))<derivlimit) {
	  flattendone=true;
	  if (outputiter) cout << endl;
	} else {
	  // Why does this work?
	  // The alternative of simply extending the LHS doesn't 
	  // seem to work.
	  if (rp.convergeflag==true) {
	    for(int i=1;i<=ngrid;i++) {
	      rel->x[i]*=fact;
	    }
	  }
	}

	if (iterfile) {
	  string soutt=((string)argc[1])+"/it"+ itos(j)+".out";
	  fout.open(soutt.c_str());
	  fout.setf(ios::scientific);
	  fout << "x sigma omega rho sigmap omegap rhop" << endl;
	  for(int i=1;i<=ngrid;i++) {
	    fout << rel->x[i] << " ";
	    for(int k=1;k<=6;k++) fout << rel->y(k,i) << " ";
	    fout << endl;
	  }
	  fout.close();
	}
      }
      lastit=j-1;

      if (flattendone==false) {
	cout << "Could not get vanishing derivative at boundary." << endl;
      } else if (rp.convergeflag==false) {
	cout << "Relaxation did not converge." << endl;

	//--------------------------------------------
	// If convergence failed, use previous result
	// for calculations

	for(int i=1;i<=ngrid;i++) {
	  rel->x[i]=xstor[i];
	  for(int k=1;k<=ne;k++) {
	    rel->y(k,i)=ystor[k][i];
	  }
	}
    
      } else {
	cout << "Going to final solution." << endl;
	conve=finalconverge;
	rel->solve(conve,slowc);
      }

      //--------------------------------------------
      // Output characteristics of final solution
    
      at.line_of_names(((string)"nn np n alpha nprime esurf esurf2 ebulk ")+
		       "egrad thickint wdint wd2int qpq ");
      at.set_nlines(ngrid);

      //--------------------------------------------
      // Calculate rhon, rhop, energy, and ebulk at
      // every point in profile
    
      hns=rmf_eos.fesym(np.nsat,protfrac);
      delta=1.0-2.0*protfrac;
    
      for(int i=1;i<=ngrid;i++) {
	mstar=rmf_eos.mnuc-gs*rel->y(1,i);
    
	if (rp.mun-gw*rel->y(2,i)+gr*rel->y(3,i)/2.0<mstar) kfn=0.0;
	else kfn=sqrt(pow(rp.mun-gw*rel->y(2,i)+gr*rel->y(3,i)/2.0,2.0)-
		      mstar*mstar);
	if (rp.mup-gw*rel->y(2,i)-gr*rel->y(3,i)/2.0<mstar) kfp=0.0;
	else kfp=sqrt(pow(rp.mup-gw*rel->y(2,i)-gr*rel->y(3,i)/2.0,2.0)-
		      mstar*mstar);

	at.set("nn",i,pow(kfn,3.0)/3.0/pi2);
	at.set("np",i,pow(kfp,3.0)/3.0/pi2);
	at.set("n",i,at.get("nn",i)+at.get("np",i));
	at.set("alpha",i,at.get("nn",i)-at.get("np",i));

	np.n.ms=mstar;
	np.p.ms=mstar;
	np.n.nu=sqrt(kfn*kfn+np.n.ms*np.n.ms);
	np.p.nu=sqrt(kfp*kfp+np.p.ms*np.p.ms);
	rmf_eos.calc_eq_p(np.n,np.p,rel->y(1,i),rel->y(2,i),
			  rel->y(3,i),f1,f2,f3,np.hb);
	np.hb.ed=-np.hb.pr+np.n.n*np.n.mu+np.p.n*np.p.mu;

	at.set("ebulk",i,np.hb.ed-rp.mup*at.get("np",i)-
	       rp.mun*at.get("nn",i));
	at.set("egrad",i,0.5*(rel->y(4,i)*rel->y(4,i)-
			      rel->y(5,i)*rel->y(5,i)-
			      rel->y(6,i)*rel->y(6,i)));

	at.set("esurf",i,at.get("ebulk",i)+at.get("egrad",i));
	if (at.get("n",i)>0.0) {
	  at.set("esurf2",i,at.get("esurf",i));
	} else {
	  at.set("esurf2",i,0.0);
	}
      
	// There doesn't seem to be much difference between this approach
	// and the one below
	rhoprime[i]=(kfn*(np.n.nu*(-gw*rel->y(5,i)+gr/2.0*rel->y(6,i))+
			  mstar*gs*rel->y(4,i))+
		     kfp*(np.p.nu*(-gw*rel->y(5,i)-gr/2.0*rel->y(6,i))+
			  mstar*gs*rel->y(4,i)))/pi2;
	//if (i==0) rhoprime[i]=0;
	//else rhoprime[i]=(rho[i]-rho[i-1])/(rel->x[i]-rel->x[i-1]);
    
	if (rhoprime[i]!=0.0) {
	  // Here qpq=Q_nn+Q_np is calculated from the bulk energy. One
	  // Could calculate this from the gradient part of the surface
	  // energy as well. There is not much difference.
	  qpq[i]=(np.hb.ed-rp.mup*at.get("np",i)-rp.mun*at.get("nn",i))*4.0/rhoprime[i]/
	    rhoprime[i];
	  //qpq[i]=(rel->y(4,i)*rel->y(4,i)-
	  // rel->y(5,i)*rel->y(5,i)-
	  // rel->y(6,i)*rel->y(6,i))*2.0/rhoprime[i]/rhoprime[i];
	} else {
	  qpq[i]=0.0;
	}

	at.set("thickint",i,(at.get("nn",i)/nn0-at.get("np",i)/np0));
      
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

      mstar=rmf_eos.mnuc-gs*rp.phi0;
      kfn=sqrt(pow(rp.mun-gw*rp.v0+gr*rp.r0/2.0,2.0)-mstar*mstar);
      kfp=sqrt(pow(rp.mup-gw*rp.v0-gr*rp.r0/2.0,2.0)-mstar*mstar);

      for(int i=1;i<=3;i++) {
	xn[i]=lookup(ngrid,pow(kfn,3.0)/3.0/pi2/10.0*((double)(i*4-3)),
		     rel->x,at.get_column("nn"));
	xp[i]=lookup(ngrid,pow(kfp,3.0)/3.0/pi2/10.0*((double)(i*4-3)),
		     rel->x,at.get_column("np"));
      }

      cout << "Integrals." << endl;
      o2scl::interp<std::vector<double>,std::vector<double> > gi;

      const std::vector<double> &xav=at[0];
      size_t istt=at.get_nlines()-1;
      
      surf=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		    at.get_column("esurf"));
      surf2=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		    at.get_column("esurf2"));
      sbulk=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		    at.get_column("ebulk"));
      sgrad=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		    at.get_column("egrad"));
      thick=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		    at.get_column("thickint"));

      if (pf_index>=2) {

	wd=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		       at.get_column("wdint"));
	wd*=hns;

	wd2=gi.integ(xav[0],xav[istt],at.get_nlines(),at.get_column(0),
		     at.get_column("wd2int"));
	
	sssv1=4.0*pi*r0*r0*wd/hns;
      } else {
	w0=surf;
	w02=surf2;
      }

#ifdef NEVER_DEFINED
      
      //--------------------------------------------
      // Jim's formulas for surface tension, etc.

      if (pf_index==1) {
	w0jl=0.0;
	wdjl=0.0;

	inte_qag_gsl gl;
	double xint=0.5*np.nsat, tweight;
      
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
	
	  np.n.n=lookup(ngrid,nint,rhon,rho);
	  np.p.n=lookup(ngrid,nint,rhop,rho);
	  sig=lookup(ngrid,nint,ar1,rho);
	  ome=lookup(ngrid,nint,ar2,rho);
	  rhof=lookup(ngrid,nint,ar3,rho);
	  qq=lookup(ngrid,nint,qpq,rho);

	  np.n.kffromden();
	  np.p.kffromden();
	  np.n.ms=np.n.m-gs*sig;
	  np.p.ms=np.n.ms;
	  np.n.nu=sqrt(np.n.kf*np.n.kf+np.n.ms*np.n.ms);
	  np.p.nu=sqrt(np.p.kf*np.p.kf+np.p.ms*np.p.ms);
	
	  double estmp;
	
	  np.rmf_eos.calc_p(np.n,np.p,sig,ome,rhof,f1,f2,f3,np.hb);
	  //np.hb.ed=-np.hb.pr+np.n.n*np.n.mu+np.p.n*np.p.mu;
	
	  estmp=rmf_eos.fesym(nint);
	  // fesym() automatically computes the bulk energy density
	  // and puts the result into tht.ed:
	  np.hb.ed=tht.ed;

	  if (np.hb.ed-rp.mun*np.n.n-rp.mup*np.p.n>0.0 && qq>0.0) {
	    wdjl+=tweight*sqrt(qq)*nint*(hns/estmp-1.0)/
	      sqrt(np.hb.ed-rp.mun*np.n.n-rp.mup*np.p.n);
	    w0jl+=tweight*sqrt(fabs(qq*(np.hb.ed-rp.mun*np.n.n-
					rp.mup*np.p.n)));
	    if (!finite(wdjl)) {
	      cout << "wdjl not finite." << endl;
	      cout << np.hb.ed-rp.mun*np.n.n-rp.mup*np.p.n << endl;
	      cout << hns << " " << qq << endl;
	      exit(-1);
	    }
	  }
	}
      
	//delete gl;
	wdjl*=hns/2.0;
	sssv2=4*pi*r0*r0*wdjl/hns;
      }
      
#endif
      
    }

    delete rel;

    return 0;
  }

  int nucmat(size_t nv, const ubvector &ex, ubvector &ey) {
    double f1,f2,f3,sig,ome,rho;
    char ch;
    newtp *sp;//=(newtp *)pa;
    fermion &n=sp->n, &p=sp->p;

    n.nu=ex[0];
    p.nu=ex[1];
    sig=ex[2];
    ome=ex[3];
    rho=ex[4];

    sp->rmf_eos.calc_eq_p(n,p,sig,ome,rho,f1,f2,f3,sp->hb);

    ey[0]=p.n+n.n-sp->nsat;
    ey[1]=p.n-sp->nsat*sp->protfrac;
    ey[2]=f1;
    ey[3]=f2;
    ey[4]=f3;

    return 0;
  }

  int ndripfun(size_t sn, const ubvector &sx, ubvector &sy) {
    double pleft, pright, munleft, munright;
  
    newtp *sp;//=(newtp *)pa;
    fermion &n=sp->n, &p=sp->p;
    thermo &hb=sp->hb;

    n.n=sx[0];
    p.n=sx[1];
    sp->rmf_eos.calc_e(n,p,hb);
  
    sy[0]=p.n-sp->protfrac*(p.n+n.n);
    pleft=hb.pr;
    munleft=n.mu;

    n.n=sx[2];
    p.n=0.001;
    sp->rmf_eos.calc_e(n,p,hb);

    pright=hb.pr;
    munright=n.mu;

    sy[1]=pleft-pright;
    sy[2]=munleft-munright;

    return 0;
  }

};

int s3relax::iter(int k, double err, double fac, ubvector_int &kmax,
		  ubvector &ermax) {
  ofstream itout;
  relaxp *rp;//=(relaxp *)pa;
  char ch;

  if (rp->showrelax) {
    if (k==1) {
      cout << "\tIter.\tError\t\tFac.         Center       rhs" << endl;
    }
    if (k>=1 && k<=itmax) {
      cout << "\t" << k << "  \t" << err << "  \t" << fac << " ";
      //      cout << lookup(ngrid,y(1,1)/2.0,x,(*y)[1]) << " " 
      cout << x[ngrid] << endl;
    }
  }

  if (k>=itmax) {
    cout << "Didn't converge. k: " << k << " itmax: " << itmax << endl;
    rp->convergeflag=false;
  } else {
    rp->convergeflag=true;
  }

  if (err>rp->maxdev) {
    cout << "Reached maxdev." << endl;
    exit(-1);
    rp->convergeflag=false;

    for(int i=1;i<=ngrid;i++) {
      x[i]=6.0*((double)(i-ngrid/2))/ngrid;
      //      y[i]=rp->phi0;
    }

  }

  return 0;
}

int s3relax::difeq(int k, int k1, int k2, int jsf, int is1, int isf) {
  double phi, v, phip, vp, fn, fn2, dx, mstar, rp, r;
  double fn3, kfn, kfp;
  double gs, gw, gr, hbarn;
  relaxp *rpa;//=(relaxp *)pa;
  char ch;

  gs=rpa->rmf_eos.cs*rpa->rmf_eos.ms;
  gw=rpa->rmf_eos.cw*rpa->rmf_eos.mw;
  gr=rpa->rmf_eos.cr*rpa->rmf_eos.mr;
  
  // 10/9/03 - Seems like this is repetitive. 
  // s[4][ne+indexv[1]] should be calculated already in calcderiv??
  if (k==k1) {
    if (nb==3) {
      s(4,ne+indexv[1])=1.0;
      s(4,ne+indexv[2])=0.0;
      s(4,ne+indexv[3])=0.0;
      s(4,ne+indexv[4])=0.0;
      s(4,ne+indexv[5])=0.0;
      s(4,ne+indexv[6])=0.0;
      s(4,jsf)=y(1,1)-rpa->phi0;
      s(5,ne+indexv[1])=0.0;
      s(5,ne+indexv[2])=1.0;
      s(5,ne+indexv[3])=0.0;
      s(5,ne+indexv[4])=0.0;
      s(5,ne+indexv[5])=0.0;
      s(5,ne+indexv[6])=0.0;
      s(5,jsf)=y(2,1)-rpa->v0;
      s(6,ne+indexv[1])=0.0;
      s(6,ne+indexv[2])=0.0;
      s(6,ne+indexv[3])=1.0;
      s(6,ne+indexv[4])=0.0;
      s(6,ne+indexv[5])=0.0;
      s(6,ne+indexv[6])=0.0;
      s(6,jsf)=y(3,1)-rpa->r0;
    } else {
      cout << "What?" << endl;
      exit(-1);
    }

  } else if (k>k2) {
    if (nb==3) {
      s(1,ne+indexv[1])=1.0;
      s(1,ne+indexv[2])=0.0;
      s(1,ne+indexv[3])=0.0;
      s(1,ne+indexv[4])=0.0;
      s(1,ne+indexv[5])=0.0;
      s(1,ne+indexv[6])=0.0;
      s(1,jsf)=y(1,ngrid);
      s(2,ne+indexv[1])=0.0;
      s(2,ne+indexv[2])=1.0;
      s(2,ne+indexv[3])=0.0;
      s(2,ne+indexv[4])=0.0;
      s(2,ne+indexv[5])=0.0;
      s(2,ne+indexv[6])=0.0;
      s(2,jsf)=y(2,ngrid);
      s(3,ne+indexv[1])=0.0;
      s(3,ne+indexv[2])=0.0;
      s(3,ne+indexv[3])=1.0;
      s(3,ne+indexv[4])=0.0;
      s(3,ne+indexv[5])=0.0;
      s(3,ne+indexv[6])=0.0;
      s(3,jsf)=y(3,ngrid);
    }

  } else {
    dx=x[k]-x[k-1];
    phi=(y(1,k)+y(1,k-1))/2.0;
    v=(y(2,k)+y(2,k-1))/2.0;
    r=(y(3,k)+y(3,k-1))/2.0;
    phip=(y(4,k)+y(4,k-1))/2.0;
    vp=(y(5,k)+y(5,k-1))/2.0;
    rp=(y(6,k)+y(6,k-1))/2.0;
    
    rpa->n.nu=rpa->mun-gw*v+0.5*gr*r;
    rpa->p.nu=rpa->mup-gw*v-0.5*gr*r;
    
    rpa->rmf_eos.calc_eq_p((rpa->n),(rpa->p),phi,v,r,fn,fn2,fn3,
			   (rpa->hb));
    s(1,jsf)=y(1,k)-y(1,k-1)-dx*phip;
    s(2,jsf)=y(2,k)-y(2,k-1)-dx*vp;
    s(3,jsf)=y(3,k)-y(3,k-1)-dx*rp;
    s(4,jsf)=y(4,k)-y(4,k-1)-dx*fn;
    s(5,jsf)=y(5,k)-y(5,k-1)-dx*fn2;
    s(6,jsf)=y(6,k)-y(6,k-1)-dx*fn3;

  }

  return 0;
}


int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  seminf_rel sn;
  sn.run(argc,argv);
  return 0;
}