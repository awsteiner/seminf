/*
  Does not work with dripped protons yet.

  Only zero temperature for now (but some code for finite 
  temperature has been added).

  11/6/03 - changed dx=14.0/((double)ngrid); to
  dx=7.0/((double)ngrid); in monfun and monfun2 because it helped
  convergence.

  Could rework integrations to utilize the facilities of table
  integration.

  Need to double check diff eq's for APR.

  12/15/03 - Implemented a new solution method for Qnn!=Qnp for when
  nndrip becomes > 0. It seems that in this case, it is better not to
  automatically adjust the scale of the x-axis, but use a gradual
  broadening similar to what we have done using the relativistic
  code. This allows calculation down to really low proton fractions!
  We need to see if this is consistent with the analytic code for both
  Qnn<Qnp (hard) and Qnn>Qnp (probably works).

  NOTE: The first proton fraction cannot be 0.5, since then 'smn'
  gets confused on further values different from 0.5

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

#include "relax.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
typedef int (*mm_funct)(size_t, const ubvector &, ubvector &);

//--------------------------------------------
// Class definition for relaxation routine

class smnrelax : public relax {
public:
  int ctr;
  double errold;
  smnrelax(int tne, int tnb, int tngrid);
  int iter(int k, double err, double fac, int *kmax, ubvector &errmax);
  int difeq(int k, int k1, int k2, int jsf, int is1, int isf);
  
};

// Constructor
smnrelax::smnrelax(int tne, int tnb, int tngrid) : relax(tne,tnb,tngrid) {
}

//--------------------------------------------
// Random global variables

const double big=3.0;
const int ngrid=100;
bool relaxfile, goodbound, qnn_equals_qnp, rhsmode, flatden;
string dirname;
double firstderiv, expo, rhslength=2.0;
double temper, nndrip, npdrip, minerr, pbpf;
ubvector xstor;
ubmatrix ystor;
double xg1[ngrid+1], xg2[ngrid+1], yg1[6][ngrid+1], yg2[6][ngrid+1];
double nnrhs, nnprhs, npprhs, nprhs, monfact, relaxconverge, rhsmin;
double initialstep;
int iend, relret1, relret2;
string model, pbmodel, hcmodel;
smnrelax *rel;
inte_qagiu_gsl<> gl2;
table<ubvector> *at;
eos_had_apr eosa;
eos_had_potential eosg;
deriv_gsl<> df;

//--------------------------------------------
// Structure for relax parameters

typedef struct relaxp_s {
  double mun, mup, nn0, np0, qnn, qpp, qnp, n0half;
  bool showrelax;
  double dndnl, dndnr, dndpl, dndpr, dpdnl, dpdnr, dpdpl, dpdpr;
  double lex1, lex2, rex1, rex2, dmundn, barn, dmundp, nsat, protfrac;
  int iii;
  eos_had_base *eos;
  fermion *n, *p;
  thermo *hb;
  mroot_hybrids<mm_funct> *nd;
} relaxp;

//--------------------------------------------
// Function prototypes

double lookup(int n, double y0, const ubvector &x, const ubvector &y);
double monfun(double sx, int ngrid, void *pa);
double monfun2(double sx, int ngrid, void *pa);
int derivs(double sx, const ubvector &sy, ubvector &dydx, void *pa);
int qnnqnpfun(size_t sn, const ubvector &sx, ubvector &sy);
int ndripfun(size_t sn, const ubvector &sx, ubvector &sy);
int pdripfun(size_t sn, const ubvector &sx, ubvector &sy);

int main(int argc, char *argv[]) {
  
  double wd=0.0, wd2=0.0, dtemp;
  double rhsn, rhsp, det;
  double surf, sbulk, sgrad, sssv_drop=0.0, *pflist;
  ubvector sx(5), sy(5);
  double sssv_jim=0.0, hns, delta, w0jl=0.0, wdjl=0.0, den, w0, xn[4], xp[4];
  double thick;
  int i, ngrid, iii, choice, ngl=32, j, npf;
  ifstream parin, guess; 
  ofstream fout;
  double lon, lop, hin, hip, sqt;
  double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;
  const ubvector *rho, *alpha, *ebulk, *egrad, *esurf, *thickint;
  const ubvector *wdint, *wd2int;
  const ubvector *vqnn, *vqnp, *vqpp;
  
  inte_qag_gsl<> gl;
  eos_had_skyrme eoss;
  eos_had_schematic eosp;
  mroot_hybrids<mm_funct> nd;
  relaxp rp;

  if (argc<2) {
    cout << "Need directory name." << endl;
    exit(-1);
  } else if (argc==3) {
    cout << "Reading inital guess from file: " << argv[2] << endl;
  }
  dirname=argv[1];

  cout.setf(ios::scientific);

  nd.tol_abs=1.0e-11;
  nd.tol_rel=1.0e-8;
  nd.ntrial=100;
  nd.verbose=0;

  rp.nd=&nd;

  initialstep=7.0;

  //--------------------------------------------
  // Get model parameters

  //if (word_io.get_def(pa,"model",model)!=0) model="skyrme";
  model="skyrme";

  eos_had_apr eosa;

  skyrme_load(eoss,"SLy4");

  /*
    if (model=="skyrme" || model=="hcskyrme") {
    eoss=new eos_had_skyrme;
    rp.eos=eoss;
    if (model=="skyrme") {
    cout << "Skyrme EOS" << endl;
    double_io.get_def(pa,"t0",eoss->t0);
    double_io.get_def(pa,"t1",eoss->t1);
    double_io.get_def(pa,"t2",eoss->t2);
    double_io.get_def(pa,"t3",eoss->t3);
    eoss->t0/=hc_mev_fm;
    eoss->t1/=hc_mev_fm;
    eoss->t2/=hc_mev_fm;
    eoss->t3/=hc_mev_fm;
    double_io.get_def(pa,"x0",eoss->x0);
    double_io.get_def(pa,"x1",eoss->x1);
    double_io.get_def(pa,"x2",eoss->x2);
    double_io.get_def(pa,"x3",eoss->x3);
    double_io.get_def(pa,"a",eoss->a);
    double_io.get_def(pa,"b",eoss->b);
    double_io.get_def(pa,"alpha",eoss->alpha);
    } else {
    cout << "Hard-coded Skyrme EOS" << endl;
    if (word_io.get_def(pa,"hcmodel",hcmodel)!=0) {
    cout << "No hard-coded model specified." << endl;
    exit(-1);
    } else {
    eoss->load(hcmodel);
    }
    }
    
    text_out_file tout(&cout);
    pa.fout(&tout);

    rp.qnn=0.1875*(eoss->t1*(1.0-eoss->x1)-eoss->t2*(1.0+eoss->x2));
    rp.qpp=rp.qnn;
    rp.qnp=0.125*(3.0*eoss->t1*(1.0+eoss->x1/2.0)-eoss->t2*(1.0+eoss->x2/2.0));

    if (bool_io.get_def(pa,"xpb",xpb)!=0) xpb=false;

  */
  
  //} else if (model=="schematic") {

  /*
    cout << "Schematic EOS" << endl;
    eosp=new eos_had_schematic;
    rp.eos=eosp;
    double_io.get_def(pa,"eoa",eosp->eoa);
    if (double_io.get_def(pa,"kprime",eosp->kprime)!=0) eosp->kprime=0.0;
    double_io.get_def(pa,"comp",eosp->comp);
    double_io.get_def(pa,"esym",eosp->esym);
    eosp->eoa/=hc_mev_fm;
    eosp->kprime/=hc_mev_fm;
    eosp->comp/=hc_mev_fm;
    eosp->esym/=hc_mev_fm;
    eosp->b=eosp->esym;
    eosp->a=0.0;
  */

  /*
    double_io.get_def(pa,"n0",eosp->n0);
    double_io.get_def(pa,"qnn",rp.qnn);
    double_io.get_def(pa,"qnp",rp.qnp);
    double_io.get_def(pa,"qpp",rp.qpp);
    rp.qnn/=hc_mev_fm;
    rp.qnp/=hc_mev_fm;
    rp.qpp/=hc_mev_fm;
  */
  
  /*
    } else if (model=="apr") {
    cout << "APR EOS" << endl;
    rp.eos=eosa;
    if (int_io.get_def(pa,"choice",choice)!=0) choice=1;
    eosa->select(choice);

    } else if (model=="gp") {
    cout << "General Potential EOS" << endl;
    eosg=new eos_had_potential;
    rp.eos=eosg;
    double_io.get_def(pa,"x",eosg.x,0.0);
    double_io.get_def(pa,"Au",eosg.Au,-95.98/hc_mev_fm);
    double_io.get_def(pa,"Al",eosg.Al,-120.57/hc_mev_fm);
    double_io.get_def(pa,"rho0",eosg.rho0,0.16);
    double_io.get_def(pa,"B",eosg.B,106.35/hc_mev_fm);
    double_io.get_def(pa,"sigma",eosg.sigma,4.0/3.0);
    double_io.get_def(pa,"Cl",eosg.Cl,-11.70/hc_mev_fm);
    double_io.get_def(pa,"Cu",eosg.Cu,-103.40/hc_mev_fm);
    double_io.get_def(pa,"Lambda",eosg.Lambda,cbrt(1.5*pi2*0.16));
    eosg.form=eosg.mdi_form;
    } else {
    cout << "Unknown EOS" << endl;
    exit(-1);
    }
  */
  
  //--------------------------------------------
  //

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
  
  flatden=false;

  //--------------------------------------------
  // Rearrange memory for relaxation

  at=new table<ubvector>(ngrid*2);
  at->line_of_names(((string)"x nn np nnp npp scale rho alpha ")+
		    "ebulk egrad esurf thickint wdint wd2int");

  rho=&at->get_column("rho");
  alpha=&at->get_column("alpha");
  ebulk=&at->get_column("ebulk");
  egrad=&at->get_column("egrad");
  esurf=&at->get_column("esurf");
  thickint=&at->get_column("thickint");
  wdint=&at->get_column("wdint");
  wd2int=&at->get_column("wd2int");

  if (model=="apr" || model=="gp") {
    at->new_column("qnn");
    at->new_column("qnp");
    at->new_column("qpp");
    vqnn=&at->get_column("qnn");
    vqnp=&at->get_column("qnp");
    vqpp=&at->get_column("qpp");
  }
  
  xstor.resize(ngrid+1);
  ystor.resize(5+1,ngrid+1);

  minerr=1.0e-6;
  relaxconverge=1.0e-9;
  rp.showrelax=true;
  relaxfile=false;
  firstderiv=-2.0e-4;
  temper=0.0;
  rhsmin=1.0e-7;
  monfact=0.002;

  rp.n=new fermion(939.0/hc_mev_fm,2.0);
  rp.p=new fermion(939.0/hc_mev_fm,2.0);
  rp.n->non_interacting=false;
  rp.p->non_interacting=false;
  rp.eos->set_n_and_p(*rp.n,*rp.p);
  rp.hb=new thermo;
  rp.eos->set_thermo(*rp.hb);
  rp.n0half=rp.eos->fn0(0.5,dtemp);
  double r0=cbrt(0.75/pi/rp.n0half);

  /*
    if (xpb==true) {
    pflist=new double[2];
    pflist[0]=0.48;
    pflist[1]=pbpf;
    npf=2;
    } else if (double_io.get(pa,"protfraclist",pflist,npf)!=0) {
    double_io.get_def(pa,"protfrac",pflist[1]);
    npf=2;
    }
  */

  pflist=new double[2];
  pflist[0]=0.5;
  npf=2;
  
  for(iii=1;iii<=npf;iii++) {
    rp.iii=iii;
    rp.protfrac=pflist[iii-1];

    //--------------------------------------------
    // Calculate properties of saturated nuclear 
    // matter with appropriate value of proton fraction
    // including the possibility of drip particles
    // if they exist
    
    rp.nsat=rp.eos->fn0(rp.protfrac,dtemp);
    rp.np0=rp.nsat*rp.protfrac;
    rp.nn0=rp.nsat-rp.np0;
    rp.n->n=rp.nn0;
    rp.p->n=rp.np0;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    rp.mun=rp.n->mu;
    rp.mup=rp.p->mu;
    nndrip=0.0;
    npdrip=0.0;

    if (rp.mun>rp.n->m) {
      cout << "Neutron drip" << endl;

      sx[0]=rp.n->n;
      sx[1]=rp.p->n;
      sx[2]=0.01;
      mm_funct ndf(ndripfun);
      nd.msolve(3,sx,ndf);
      ndripfun(3,sx,sy);
      if (fabs(sy[0])>1.0e-4 || fabs(sy[1])>1.0e-4 || fabs(sy[2])>1.0e-4) {
	cout << "Error in solution of ndripfun." << endl;
	exit(-1);
      }

      rp.nn0=sx[0];
      rp.np0=sx[1];
      rp.n->n=rp.nn0;
      rp.p->n=rp.np0;
      rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
      rp.mun=rp.n->mu;
      rp.mup=rp.p->mu;
      nndrip=sx[2];
      npdrip=0.0;

    }

    if (rp.mup>rp.p->m) {
      cout << "Proton drip" << endl;

      sx[0]=rp.n->n;
      sx[1]=rp.p->n;
      sx[2]=0.01;
      mm_funct pdf(pdripfun);
      nd.msolve(3,sx,pdf);
      pdripfun(3,sx,sy);
      if (fabs(sy[0])>1.0e-4 || fabs(sy[1])>1.0e-4 || fabs(sy[2])>1.0e-4) {
	cout << "Error in solution of pdripfun." << endl;
	exit(-1);
      }

      rp.nn0=sx[0];
      rp.np0=sx[1];
      rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
      rp.mun=rp.n->mu;
      rp.mup=rp.p->mu;
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
   
    rp.n->n=rp.nn0;
    rp.p->n=rp.np0;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    hin=rp.n->mu;

    rp.n->n=rp.nn0-1.0e-4;
    rp.p->n=rp.np0;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    lon=rp.n->mu;
    
    rp.dmundn=1.0e4*(hin-lon);
  
    // dmundp
  
    rp.n->n=rp.nn0;
    rp.p->n=rp.np0;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    hin=rp.n->mu;

    rp.n->n=rp.nn0;
    rp.p->n=rp.np0-1.0e-4;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    lon=rp.n->mu;
    
    rp.dmundp=1.0e4*(hin-lon);
  
    //--------------------------------------------
    // Calculate derivatives of the RHS's numerically
  
    // At the left hand side
    rp.n->n=rp.nn0;
    rp.p->n=rp.np0;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    if (model=="apr") {
      eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
    } 
    if (model=="gp") {
      thermo thth;
      eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,
			 dqnpdnp,dqppdnn,dqppdnp);
    } 
    det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
    if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
      rhsn=(rp.n->mu-rp.mun)/rp.qpp/2.0;
      rhsp=(rp.p->mu-rp.mup)/rp.qnn/2.0;
    } else {
      rhsn=(rp.n->mu-rp.mun)*rp.qpp-(rp.p->mu-rp.mup)*rp.qnp;
      rhsp=(rp.p->mu-rp.mup)*rp.qnn-(rp.n->mu-rp.mun)*rp.qnp;
      rhsn/=det;
      rhsp/=det;
    }
    hin=rhsn;
    hip=rhsp;

    rp.n->n=rp.nn0-1.0e-4;
    rp.p->n=rp.np0;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    if (model=="apr") {
      eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
    } 
    if (model=="gp") {
      thermo thth;
      eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,
			 dqnpdnp,dqppdnn,dqppdnp);
    } 
    det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
    if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
      rhsn=(rp.n->mu-rp.mun)/rp.qpp/2.0;
      rhsp=(rp.p->mu-rp.mup)/rp.qnn/2.0;
    } else {
      rhsn=(rp.n->mu-rp.mun)*rp.qpp-(rp.p->mu-rp.mup)*rp.qnp;
      rhsp=(rp.p->mu-rp.mup)*rp.qnn-(rp.n->mu-rp.mun)*rp.qnp;
      rhsn/=det;
      rhsp/=det;
    }
    lon=rhsn;
    lop=rhsp;

    rp.dndnl=1.0e4*(hin-lon);
    rp.dpdnl=1.0e4*(hip-lop);
  
    rp.n->n=rp.nn0;
    rp.p->n=rp.np0-1.0e-4;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    if (model=="apr") {
      eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
    } 
    if (model=="gp") {
      thermo thth;
      eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,
			 dqnpdnp,dqppdnn,dqppdnp);
    } 
    det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
    if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
      rhsn=(rp.n->mu-rp.mun)/rp.qpp/2.0;
      rhsp=(rp.p->mu-rp.mup)/rp.qnn/2.0;
    } else {
      rhsn=(rp.n->mu-rp.mun)*rp.qpp-(rp.p->mu-rp.mup)*rp.qnp;
      rhsp=(rp.p->mu-rp.mup)*rp.qnn-(rp.n->mu-rp.mun)*rp.qnp;
      rhsn/=det;
      rhsp/=det;
    }
    lon=rhsn;
    lop=rhsp;
  
    rp.dndpl=1.0e4*(hin-lon);
    rp.dpdpl=1.0e4*(hip-lop);
  
    // At the right hand side
    rp.n->n=1.0e-4;
    rp.p->n=1.0e-4;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    if (model=="apr") {
      eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
    } 
    if (model=="gp") {
      thermo thth;
      eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,
			 dqnpdnp,dqppdnn,dqppdnp);
    } 
    det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
    if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
      rhsn=(rp.n->mu-rp.mun)/rp.qpp/2.0;
      rhsp=(rp.p->mu-rp.mup)/rp.qnn/2.0;
    } else {
      rhsn=(rp.n->mu-rp.mun)*rp.qpp-(rp.p->mu-rp.mup)*rp.qnp;
      rhsp=(rp.p->mu-rp.mup)*rp.qnn-(rp.n->mu-rp.mun)*rp.qnp;
      rhsn/=det;
      rhsp/=det;
    }
    hin=rhsn;
    hip=rhsp;
  
    rp.n->n=1.0e-8;
    rp.p->n=1.0e-4;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    if (model=="apr") {
      eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
    } 
    if (model=="gp") {
      thermo thth;
      eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,
			 dqnpdnp,dqppdnn,dqppdnp);
    } 
    det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
    if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
      rhsn=(rp.n->mu-rp.mun)/rp.qpp/2.0;
      rhsp=(rp.p->mu-rp.mup)/rp.qnn/2.0;
    } else {
      rhsn=(rp.n->mu-rp.mun)*rp.qpp-(rp.p->mu-rp.mup)*rp.qnp;
      rhsp=(rp.p->mu-rp.mup)*rp.qnn-(rp.n->mu-rp.mun)*rp.qnp;
      rhsn/=det;
      rhsp/=det;
    }
    lon=rhsn;
    lop=rhsp;
  
    rp.dndnr=1.0e4*(hin-lon);
    rp.dpdnr=1.0e4*(hip-lop);
  
    rp.n->n=1.0e-4;
    rp.p->n=1.0e-8;
    rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);
    if (model=="apr") {
      eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			  dqnndnn,dqnndnp,dqnpdnn,
			  dqnpdnp,dqppdnn,dqppdnp);
    } 
    if (model=="gp") {
      thermo thth;
      eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			 dqnndnn,dqnndnp,dqnpdnn,
			 dqnpdnp,dqppdnn,dqppdnp);
    } 
    det=rp.qnn*rp.qpp-rp.qnp*rp.qnp;
    if (qnn_equals_qnp==true || (model=="apr" && fabs(det)<1.0e-5)) {
      rhsn=(rp.n->mu-rp.mun)/rp.qpp/2.0;
      rhsp=(rp.p->mu-rp.mup)/rp.qnn/2.0;
    } else {
      rhsn=(rp.n->mu-rp.mun)*rp.qpp-(rp.p->mu-rp.mup)*rp.qnp;
      rhsp=(rp.p->mu-rp.mup)*rp.qnn-(rp.n->mu-rp.mun)*rp.qnp;
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

    if (!finite(rp.lex1)) {
      goodbound=false;
      rp.lex1=0.0;
    }
    if (!finite(rp.lex2)) {
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
  
    if (!finite(rp.rex1)) {
      goodbound=false;
      rp.rex1=0.0;
    }
    if (!finite(rp.rex2)) {
      goodbound=false;
      rp.rex2=0.0;
    }

    //--------------------------------------------
    // Solve

    if (qnn_equals_qnp) {
      monfun2(monfact,ngrid,&rp);
    } else {
      monfun(monfact,ngrid,&rp);
    }

    //--------------------------------------------
    // Calculate the value of x for nn*0.9, nn*0.5, 
    // nn*0.1, etc.

    cout << "Lookups. " << endl;
    for(i=1;i<=3;i++) {
      xn[i]=lookup(at->get_nlines(),
		   nndrip+(rp.nn0-nndrip)/10.0*((double)(i*4-3)),
		   (*at)[0],(*at)[1]);
      xp[i]=lookup(at->get_nlines(),
		   npdrip+(rp.np0-npdrip)/10.0*((double)(i*4-3)),
		   (*at)[0],(*at)[2]);
    }

    delta=1.0-2.0*rp.protfrac;
    rp.barn=rp.nsat;
    hns=rp.eos->fesym(rp.barn);

    cout << "Integrands." << endl;
    for(i=0;i<((int)at->get_nlines());i++) {
      // Fix negative values
      if ((*at)[1][i]<0.0) at->set(1,i,0.0);
      if ((*at)[2][i]<0.0) at->set(2,i,0.0);

      rp.n->n=(*at)[1][i];
      rp.p->n=at->get(2,i);

      rp.eos->calc_e(*rp.n,*rp.p,*rp.hb);

      if (model=="apr") {
	eosa.gradient_qij2(rp.n->n,rp.p->n,rp.qnn,rp.qnp,rp.qpp,
			    dqnndnn,dqnndnp,dqnpdnn,
			    dqnpdnp,dqppdnn,dqppdnp);
	at->set("qnn",i,rp.qnn);
	at->set("qnp",i,rp.qnp);
	at->set("qpp",i,rp.qpp);
      }
      if (model=="gp") {
	thermo thth;
	eosg.gradient_qij(*rp.n,*rp.p,thth,rp.qnn,rp.qnp,rp.qpp,
			   dqnndnn,dqnndnp,dqnpdnn,
			   dqnpdnp,dqppdnn,dqppdnp);
	at->set("qnn",i,rp.qnn);
	at->set("qnp",i,rp.qnp);
	at->set("qpp",i,rp.qpp);
      } 

      at->set("rho",i,rp.n->n+rp.p->n);
      at->set("alpha",i,rp.n->n-rp.p->n);
      at->set("ebulk",i,rp.hb->ed-rp.mun*at->get(1,i)-rp.mup*at->get(2,i));
      at->set("egrad",i,0.5*rp.qnn*at->get(3,i)*at->get(3,i)+
	      0.5*rp.qpp*at->get(4,i)*at->get(4,i)+
	      rp.qnp*at->get(3,i)*at->get(4,i));
      at->set("esurf",i,at->get("ebulk",i)+at->get("egrad",i));
      at->set("thickint",i,at->get(1,i)/rp.nn0-at->get(2,i)/rp.np0);
      if (fabs(rp.protfrac-0.5)>1.0e-8) {
	at->set("wdint",i,(*alpha)[i]/delta-(*rho)[i]);
	if ((*rho)[i]>=1.0e-5) {
	  rp.barn=(*rho)[i];
	  at->set("wd2int",i,(*rho)[i]*(pow((*alpha)[i]/delta/(*rho)[i],2.0)*
					rp.eos->fesym(rp.barn)-hns));
	}
      }
    }
    
    cout << "Integrals." << endl;

    interp<> gi;
    surf=gi.integ(at->get(0,0),(*at)[0][at->get_nlines()-1],at->get_nlines(),
		  (*at)[0],(at->get_column("esurf")));
    sbulk=gi.integ(at->get(0,0),(*at)[0][at->get_nlines()-1],at->get_nlines(),
		   (*at)[0],(at->get_column("ebulk")));
    sgrad=gi.integ(at->get(0,0),(*at)[0][at->get_nlines()-1],at->get_nlines(),
		   (*at)[0],(at->get_column("egrad")));
    thick=gi.integ(at->get(0,0),(*at)[0][at->get_nlines()-1],at->get_nlines(),
		   (*at)[0],(at->get_column("thickint")));
    cout << "thick: " << thick << endl;
    if (fabs(rp.protfrac-0.5)>1.0e-8) {
      wd=gi.integ(at->get(0,0),(*at)[0][at->get_nlines()-1],at->get_nlines(),
		  (*at)[0],(at->get_column("wdint")));
      wd*=hns;
      
      wd2=gi.integ(at->get(0,0),(*at)[0][at->get_nlines()-1],at->get_nlines(),
		   (*at)[0],(at->get_column("wd2int")));
      
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
      at->new_column("nnpp");
      at->new_column("nppp");
      at->new_column("rhsn");
      at->new_column("rhsp");
      at->new_column("lhs_n");
      at->new_column("lhs_a");
      at->new_column("rhs_n");
      at->new_column("rhs_a");
      at->set("nnpp",0,0.0);
      at->set("nppp",0,0.0);
      at->set("rhsn",0,0.0);
      at->set("rhsp",0,0.0);
      at->set("lhs_n",0,0.0);
      at->set("lhs_a",0,0.0);
      at->set("rhs_n",0,0.0);
      at->set("rhs_a",0,0.0);

      rhsmode=false;

      for(i=1;i<((int)at->get_nlines());i++) {

	at->set("nnpp",i,(at->get("nnp",i)-at->get("nnp",i-1))/
		(at->get("x",i)-at->get("x",i-1)));
	at->set("nppp",i,(at->get("npp",i)-at->get("npp",i-1))/
		(at->get("x",i)-at->get("x",i-1)));
	
	tx=(at->get("x",i)+at->get("x",i-1))/2.0;
	ty[1]=(at->get("nn",i)+at->get("nn",i-1))/2.0;
	ty[2]=(at->get("np",i)+at->get("np",i-1))/2.0;
	ty[3]=(at->get("nnp",i)+at->get("nnp",i-1))/2.0;
	ty[4]=(at->get("npp",i)+at->get("npp",i-1))/2.0;

	if (ty[2]<=0.0) rhsmode=true;
	derivs(tx,ty,tdy,&rp);

	at->set("rhsn",i,tdy[3]);
	at->set("rhsp",i,tdy[4]);

	if (ty[2]>0.0) {
	  at->set("lhs_n",i,(rp.n->mu-rp.mun+rp.p->mu-rp.mup)/2.0);
	  at->set("lhs_a",i,(rp.n->mu-rp.mun-rp.p->mu+rp.mup)/2.0);
	  at->set("rhs_n",i,(at->get("nnpp",i)+at->get("nppp",i))/2.0*
		  (rp.qnn+rp.qnp));
	  at->set("rhs_a",i,(at->get("nnpp",i)-at->get("nppp",i))/2.0*
		  (rp.qnn-rp.qnp));
	} else {
	  at->set("lhs_n",i,rp.n->mu-rp.mun);
	  at->set("rhs_n",i,at->get("nnpp",i)*(rp.qnn));
	  at->set("lhs_a",i,0.0);
	  at->set("rhs_a",i,0.0);
	}
      }
    }

#ifdef REMOVING_OLD_OUTPUT
    
    //--------------------------------------------
    // Make solution output file
    
    table_io.out_one(((string)argc[1])+"/final.out"+itos(iii),
		     "at",at);
    
    //--------------------------------------------
    // Make summary output file

    cout << "Summary output." << endl;
    string soutt=((string)argc[1])+"/sisum.out"+itos(iii);
    fout.open(soutt.c_str());
    fout.setf(ios::scientific);
    fout << "int relret1 " << relret1 << endl;
    fout << "int relret2 " << relret2 << endl;
    fout << "double nsat " << rp.nsat << endl;
    fout << "double nn0 " << rp.nn0 << endl;
    fout << "double np0 " << rp.np0 << endl;
    fout << "double mun " << rp.mun << endl;
    fout << "double mup " << rp.mup << endl;
    if (model!="apr" && model!="gp") {
      fout << "double qnn " << rp.qnn << endl;
      fout << "double qpp " << rp.qpp << endl;
      fout << "double qnp " << rp.qnp << endl;
    }
    fout << "double nnp " << (*at)[3][1] << endl;
    fout << "double npp " << (*at)[4][1] << endl;
    fout << "double Tn " << xn[1]-xn[3] << endl;
    fout << "double Tp " << xp[1]-xp[3] << endl;
    fout << "double XnXp " << xn[2]-xp[2] << endl;
    fout << "double XnXp2 " << thick << endl;
    fout << "double surf " << surf << endl;
    fout << "double sbulk " << sbulk << endl;
    fout << "double sgrad " << sgrad << endl;
    fout << "double lex1 " << rp.lex1 << endl;
    fout << "double lex2 " << rp.lex2 << endl;
    fout << "double expo " << expo << endl;
    fout << "double rex1 " << rp.rex1 << endl;
    fout << "double rex2 " << rp.rex2 << endl;
    fout << "double monfact " << monfact << endl;
    fout << "double nnrhs " << nnrhs << endl;
    fout << "double nprhs " << nprhs << endl;
    fout << "double rhsmin " << rhsmin << endl;
    fout << "double nndrip " << nndrip << endl;
    fout << "double npdrip " << npdrip << endl;

    rp.barn=rp.n0half;
    fout << "double esym_n0 " << rp.eos->fesym(rp.n0half) << endl;
    fout << "double esym_ns " << hns << endl;

    fout << "bool xpb " << xpb << endl;
    if (xpb==true && iii==2) {
      fout << "double pbpf " << pbpf << endl;
    }
    if (fabs(rp.protfrac-0.5)>1.0e-8) {
      fout << "double wd_droplet " << wd << endl;
      fout << "double wd_droplet2 " << wd2 << endl;
      fout << "double wdjl " << wdjl << endl;
      fout << "double w0jl " << w0jl << endl;
      fout << "double sssv_drop " << sssv_drop << endl;
      fout << "double sssv_jim " << sssv_jim << endl;
    }
    fout << "double protfrac " << rp.protfrac << endl;
    fout << "bool flatden " << flatden << endl;

    fout.close();

#endif
    
  }
  
  delete at;

  delete rp.hb;
  delete rp.n;
  delete rp.p;
  
  return 0;
}
  
//--------------------------------------------
// Code to be executed every iteration

int smnrelax::iter(int k, double err, double fac, int *kmax, ubvector &ermax) {
  ofstream itout;
  relaxp *rp;//=(relaxp *)pa;
  char s1[80];

  int i,j;

  //--------------------------------------------
  // Output relaxation steps to file
  
  if (true || relaxfile) {
    string soutt=dirname+"/rel"+itos(k)+".out"+
      itos(rp->iii);
    itout.open(soutt.c_str());
    itout.setf(ios::scientific);
    for(i=1;i<=ngrid;i++) {
      itout << (*x)[i] << " ";
      for(j=1;j<=ne;j++) itout << (*y)(j,i) << " ";
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
      cout << "\t" << k << "  \t" << err << "  \t" << (*y)(3,1) << " ";
      ubmatrix_row amc(*y,1);
      cout << lookup(ngrid,(*y)(1,1)/2.0,*x,amc) << " " 
	   << (*x)[ngrid] << endl;
    }
  }
  
  //--------------------------------------------
  // Correct if density is negative
  
  if (!rhsmode && !qnn_equals_qnp) {
    for(i=1;i<=ngrid;i++) {
      if ((*rel->y)(1,i)<0.0) (*rel->y)(3,i)=0.0;
      if ((*rel->y)(2,i)<0.0) (*rel->y)(4,i)=0.0;
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

double lookup(int n, double yy0, const ubvector &x, const ubvector &y) {
  double x0;
  int i;

  for(i=2;i<=n;i++) {
    if ((y[i]>=yy0 && y[i-1]<yy0) || (y[i]<yy0 && y[i-1]>=yy0)) {
      x0=x[i-1]+(x[i]-x[i-1])*(yy0-y[i-1])/(y[i]-y[i-1]);
      return x0;
    }
  }

  return 0.0;
}

int derivs(double sx, const ubvector &sy, ubvector &dydx, void *pa) {
  relaxp *rp=(relaxp *)pa;
  double rhsn, rhsp=0.0, det;
  double dqnndnn, dqnndnp, dqnpdnn, dqnpdnp, dqppdnn, dqppdnp;

  rp->n->n=sy[1];
  rp->p->n=sy[2];
  if (model=="apr") {
    eosa.gradient_qij2(rp->n->n,rp->p->n,rp->qnn,rp->qnp,rp->qpp,
			dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
  }
  if (model=="gp") {
    thermo thth;
    if (false) {
      eosg.gradient_qij(*rp->n,*rp->p,thth,rp->qnn,rp->qnp,rp->qpp,
			 dqnndnn,dqnndnp,dqnpdnn,dqnpdnp,dqppdnn,dqppdnp);
    } else {
      rp->qnn=100.0/hc_mev_fm;
      rp->qpp=100.0/hc_mev_fm;
      rp->qnp=90.0/hc_mev_fm;
      dqnndnn=0.0;
      dqnndnp=0.0;
      dqnpdnn=0.0;
      dqnpdnp=0.0;
      dqppdnn=0.0;
      dqppdnp=0.0;
    }
  }
  
  rp->eos->calc_e(*rp->n,*rp->p,*rp->hb);
  /*
    double mun=rp->n->mu, mup=rp->p->mu, msn=rp->n->ms, msp=rp->p->ms;
    double hed=rp->hb->ed, hpr=rp->hb->pr, part;
    if (rp->iii==1) part=1000.0;
    else if (rp->iii==2) part=100.0;
    else if (rp->iii==3) part=40.0;
    else if (rp->iii==4) part=20.0;
    else if (rp->iii==5) part=5.0;
    else if (rp->iii==6) part=3.0;
    else if (rp->iii==7) part=1.0;
    else part=0.0;
    eosa.calc_e(*rp->n,*rp->p,*rp->hb);
    
    rp->n->mu=(rp->n->mu*part+mun)/(part+1.0);
    rp->p->mu=(rp->p->mu*part+mup)/(part+1.0);
    rp->n->ms=(rp->n->ms*part+msn)/(part+1.0);
    rp->p->ms=(rp->p->ms*part+msp)/(part+1.0);
    rp->hb->ed=(rp->hb->ed*part+hed)/(part+1.0);
    rp->hb->pr=(rp->hb->pr*part+hpr)/(part+1.0);
  */

  if (rp->n->n<=0.0) {
    if (rp->p->n<=0.0) {
      dydx[1]=0.0;
      dydx[2]=0.0;
      dydx[3]=0.0;
      dydx[4]=0.0;
    } else {
      dydx[1]=0.0;
      dydx[2]=sy[4];
      dydx[3]=0.0;
      if (model!="apr" && model!="gp") {
	dydx[4]=(rp->p->mu-rp->mup)/rp->qpp;
      } else {
	dydx[4]=(rp->p->mu-rp->mup+0.5*dqppdnp*sy[4]*sy[4])/rp->qpp;
      }
    }
  } else if (rp->p->n<=0.0) {
    dydx[1]=sy[3];
    dydx[2]=0.0;
    if (model!="apr" && model!="gp") {
      dydx[3]=(rp->n->mu-rp->mun)/rp->qnn;
    } else {
      dydx[3]=(rp->n->mu-rp->mun+0.5*dqnndnn*sy[3]*sy[3])/rp->qnn;
    }
    dydx[4]=0.0;
  } else {
    if (rhsmode==false) {
      if (model!="apr" && model!="gp") {
	det=rp->qnn*rp->qpp-rp->qnp*rp->qnp;
	rhsn=(rp->n->mu-rp->mun)*rp->qpp-(rp->p->mu-rp->mup)*rp->qnp;
	rhsp=(rp->p->mu-rp->mup)*rp->qnn-(rp->n->mu-rp->mun)*rp->qnp;
	rhsn/=det;
	rhsp/=det;
      } else {
	det=rp->qnn*rp->qpp-rp->qnp*rp->qnp;
	rhsn=(rp->n->mu-rp->mun-0.5*dqnndnn*sy[3]*sy[3]+dqnndnp*sy[3]*sy[4]-
	      dqnpdnp*sy[4]*sy[4]+0.5*dqppdnn*sy[4]*sy[4])*rp->qpp-
	  (rp->p->mu-rp->mup+0.5*dqnndnp*sy[3]*sy[3]-dqnpdnn*sy[3]*sy[3]-
	   dqppdnn*sy[3]*sy[4]-0.5*dqppdnp*sy[4]*sy[4])*rp->qnp;
	rhsp=(rp->p->mu-rp->mup+0.5*dqnndnp*sy[3]*sy[3]-dqnpdnn*sy[3]*sy[3]-
	      dqppdnn*sy[3]*sy[4]-0.5*dqppdnp*sy[4]*sy[4])*rp->qnn-
	  (rp->n->mu-rp->mun-0.5*dqnndnn*sy[3]*sy[3]+dqnndnp*sy[3]*sy[4]-
	   dqnpdnp*sy[4]*sy[4]+0.5*dqppdnn*sy[4]*sy[4])*rp->qnp;
	rhsn/=det;
	rhsp/=det;
      }
    } else {
      if (model!="apr" && model!="gp") {
	rhsn=(rp->n->mu-rp->mun)/rp->qnn;
	rhsp=0.0;
      } else {
	rhsn=(rp->n->mu-rp->mun+0.5*dqnndnn*sy[3]*sy[3])/rp->qnn;
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

  if (!finite(dydx[3])) {
    cout << "3 not finite." << endl;
    dydx[3]=0.0;
  }
  if (!finite(dydx[4])) {
    cout << "4 not finite." << endl;
    dydx[4]=0.0;
  }
  
  return 0;
}

int smnrelax::difeq(int k, int k1, int k2, int jsf, int is1, int isf) {

  int i;
  double dx, exx;
  ubvector exy(8), exdy(8);
  ubvector sx(4);
  relaxp *rp;//=(relaxp *)pa;
  double delta, epsi, tmp;

  //--------------------------------------------
  // LHS for Qnn==Qnp

  if (qnn_equals_qnp==true && rhsmode==false) {

    delta=(rp->nn0+rp->np0)*monfact/exp(big*expo);
    mm_funct qqf(qnnqnpfun);

    if (k==k1) {
      (*s)(2,jsf)=(*y)(1,1)-(rp->nn0+rp->np0-delta*exp(big*expo));
      (*s)(3,jsf)=(*y)(2,1)+delta*expo*exp(big*expo);
    } else if (k>k2 && npdrip==0.0) {
      sx[1]=(1.0-rp->protfrac)*(*y)(1,ngrid);
      sx[2]=rp->protfrac*(*y)(1,ngrid);
      rp->barn=(*y)(1,ngrid);
      rp->nd->msolve(2,sx,qqf);
      (*s)(1,jsf)=sx[2];
    } else if (k>k2 && nndrip==0.0) {
      sx[1]=(1.0-rp->protfrac)*(*y)(1,ngrid);
      sx[2]=rp->protfrac*(*y)(1,ngrid);
      rp->barn=(*y)(1,ngrid);
      rp->nd->msolve(2,sx,qqf);
      (*s)(1,jsf)=sx[1];
    } else {
      dx=(*x)[k]-(*x)[k-1];
      for(i=1;i<=3;i++) exy[i]=((*y)(i,k)+(*y)(i,k-1))/2.0;

      sx[1]=(1.0-rp->protfrac)*exy[1];
      sx[2]=rp->protfrac*exy[1];
      rp->barn=exy[1];
      rp->nd->msolve(2,sx,qqf);
      
      if (!finite(rp->n->mu)) {
      	rp->n->mu=rp->n->m;
      }
      if (!finite(rp->p->mu)) {
      	rp->p->mu=rp->p->m;
      }

      (*s)(1,jsf)=(*y)(1,k)-(*y)(1,k-1)-dx*exy[3]*exy[2];
      (*s)(2,jsf)=(*y)(2,k)-(*y)(2,k-1)-
	dx*exy[3]*(rp->n->mu-rp->mun)/rp->qnn;
      (*s)(3,jsf)=(*y)(3,k)-(*y)(3,k-1);
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
	(*s)(3,jsf)=(*y)(1,1)-(rp->nn0-delta*exp(big*expo));
	(*s)(4,jsf)=(*y)(2,1)-(rp->np0-epsi*exp(big*expo));
	(*s)(5,jsf)=(*y)(3,1)+delta*expo*exp(big*expo);
      } else {
	(*s)(2,jsf)=(*y)(1,1)-(rp->nn0-delta*exp(big*expo));
	(*s)(3,jsf)=(*y)(2,1)-(rp->np0-epsi*exp(big*expo));
	(*s)(4,jsf)=(*y)(3,1)+delta*expo*exp(big*expo);
	(*s)(5,jsf)=(*y)(4,1)+epsi*expo*exp(big*expo);
      }
      
    } else if (k>k2) {
      if (npdrip==0.0) {
	if (flatden) {
	  (*s)(1,jsf)=(*y)(2,ngrid);
	  (*s)(2,jsf)=(*y)(4,ngrid);
	} else {
	  (*s)(1,jsf)=(*y)(2,ngrid);
	}
      } else if (nndrip==0.0) {
	if (flatden) {
	  (*s)(1,jsf)=(*y)(1,ngrid);
	  (*s)(2,jsf)=(*y)(3,ngrid);
	} else {
	  (*s)(1,jsf)=(*y)(1,ngrid);
	}
      }
    } else {
      dx=(*x)[k]-(*x)[k-1];
      
      for(i=1;i<=5;i++) exy[i]=((*y)(i,k)+(*y)(i,k-1))/2.0;
      exx=((*x)[k]+(*x)[k-1])/2.0*exy[5];
      
      derivs(exx,exy,exdy,rp);
      exdy[5]=0.0;
      for(i=1;i<=5;i++) (*s)(i,jsf)=(*y)(i,k)-(*y)(i,k-1)-
			  dx*exy[5]*exdy[i];

    }

    //--------------------------------------------
    // RHS for both Qnn==Qnp and Qnn!=Qnp and APR

  } else {

    if (k==k1) {
      if (npdrip==0.0) (*s)(3,jsf)=(*y)(1,1)-nnrhs;
      else (*s)(3,jsf)=(*y)(1,1)-nprhs;
    } else if (k>k2) {
      if (npdrip==0.0) {
	(*s)(1,jsf)=(*y)(1,ngrid)-nndrip;
	if (nndrip>0.0) {
	  (*s)(2,jsf)=(*y)(3,ngrid)-rhslength;
	} else {
	  (*s)(2,jsf)=(*y)(2,ngrid);
	}
      } else {
	(*s)(1,jsf)=(*y)(1,ngrid)-npdrip;
	(*s)(2,jsf)=(*y)(2,ngrid);
      }
    } else {
      dx=(*x)[k]-(*x)[k-1];
    
      exy[1]=((*y)(1,k)+(*y)(1,k-1))/2.0;
      exy[2]=0.0;
      exy[3]=((*y)(2,k)+(*y)(2,k-1))/2.0;
      exy[4]=0.0;
      exy[5]=((*y)(3,k)+(*y)(3,k-1))/2.0;
      exx=((*x)[k]+(*x)[k-1])/2.0*exy[5];
      derivs(exx,exy,exdy,rp);

      // 10/6/04
      // I think this statement was missing in previous versions
      // and caused nasty things to happen. This should be
      // double-checked. FIXME!
      exdy[5]=0.0;

      (*s)(1,jsf)=(*y)(1,k)-(*y)(1,k-1)-dx*exy[5]*exdy[1];
      (*s)(2,jsf)=(*y)(2,k)-(*y)(2,k-1)-dx*exy[5]*exdy[3];
      (*s)(3,jsf)=(*y)(3,k)-(*y)(3,k-1)-dx*exy[5]*exdy[5];
    }
  }
  
  return 0;
}

double monfun(double lmonfact, int ngrid, void *pa) {
  relaxp *rp=(relaxp *)pa;
  bool guessdone, debug=false;
  double dx, xrhs, delta, epsi;
  ubvector y(5), dydx(5);
  int i, j, ilast=0, interpi;
  ofstream itout;
  char ch;

  monfact=lmonfact;

  //----------------------------------------------
  // Create object

  if (flatden) rel=new smnrelax(5,3,ngrid);
  else rel=new smnrelax(5,4,ngrid);
  rel->itmax=100;

  if (rp->iii==1) {
    //----------------------------------------------
    // Construct a guess by shooting
    
    xstor[1]=0.0;
    ystor(1,1)=rp->nn0;
    ystor(2,1)=rp->np0;
    ystor(3,1)=firstderiv;
    ystor(4,1)=firstderiv;
    dx=initialstep/((double)ngrid);
    
    if (debug) {
      cout << 1 << " " 
	   << xstor[1] << " " << ystor(1,1) << " " << ystor(2,1) << " "
	   << ystor(3,1) << " " << ystor(4,1) << endl;
    }
    
    guessdone=false;
    for(i=2;guessdone==false && i<=ngrid;i++) {
      xstor[i]=xstor[i-1]+dx;
      
      for(j=1;j<=4;j++) {
	y[j]=ystor(j,i-1);
      }
      derivs(xstor[i-1],y,dydx,rp);
      for(j=1;j<=4;j++) ystor(j,i)=ystor(j,i-1)+dx*dydx[j];
      
      if (ystor(1,i)<nndrip || ystor(2,i)<npdrip) {
	guessdone=true;
	ilast=i-1;
	i=ngrid+10;
      } else if (debug) {
	cout << i << " " 
	     << xstor[i] << " " << ystor(1,i) << " " << ystor(2,i) << " "
	     << ystor(3,i) << " " << ystor(4,i) << endl;
      }
    }
    
    //----------------------------------------------
    // Arrange guess into relaxation arrays and
    // stretch the solution to grid size=ngrid
    
    (*rel->x)[1]=xstor[1];
    (*rel->x)[ngrid]=xstor[ilast];
    for(i=2;i<ngrid;i++) {
      (*rel->x)[i]=(*rel->x)[1]+((double)(i-1))/((double)(ngrid-1))*
	((*rel->x)[ngrid]-(*rel->x)[1]);
    }
    interpi=1;
    for(i=1;i<=ngrid;i++) {
      while((*rel->x)[i]>xstor[interpi+1] && interpi<ilast-1) interpi++;
      while((*rel->x)[i]<xstor[interpi] && interpi>1) interpi--;
      for(j=1;j<=4;j++) {
	(*rel->y)(j,i)=ystor(j,interpi)+
	  (ystor(j,interpi+1)-ystor(j,interpi))*
	  ((*rel->x)[i]-xstor[interpi])/(xstor[interpi+1]-xstor[interpi]);
      }
      (*rel->y)(5,i)=((*rel->x)[ngrid]);
      (*rel->x)[i]=((double)(i-1))/((double)(ngrid-1));
    }

    if (false) {
      for(int iz=1;iz<=ngrid;iz++) {
	cout.precision(4);
	cout << iz << " " << (*rel->x)[iz] << " " 
	     << (*rel->y)(1,iz) << " " 
	     << (*rel->y)(2,iz) << " " 
	     << (*rel->y)(3,iz) << " " 
	     << (*rel->y)(4,iz) << " " 
	     << (*rel->y)(5,iz) << endl;
      }
    }
    
  } else {
    //----------------------------------------------
    // Use last solution for guess
    
    for(i=1;i<=ngrid;i++) {
      (*rel->x)[i]=xg1[i];
      for(j=1;j<=rel->ne;j++) {
	(*rel->y)(j,i)=yg1[j][i];
      }
    }

  }

  //----------------------------------------------
  // Solve

  rhsmode=false;
  relret1=rel->solve(relaxconverge,1.0);
  
  //----------------------------------------------
  // Copy solution for next guess

  for(i=1;i<=ngrid;i++) {
    xg1[i]=(*rel->x)[i];
    for(j=1;j<=rel->ne;j++) {
      yg1[j][i]=(*rel->y)(j,i);
    }
  }
  
  //----------------------------------------------
  // Rescale solution

  for(i=1;i<=rel->ngrid;i++) {
    (*rel->x)[i]=(*rel->x)[i]*(*rel->y)(5,i);
  }

  //----------------------------------------------
  // Store quantites at rhs for later use

  xrhs=(*rel->x)[ngrid];
  if (npdrip==0.0) {
    nnrhs=(*rel->y)(1,ngrid);
    cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
    if (nnrhs<rhsmin) {
      nnrhs=rhsmin;
      cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
    }
  } else {
    nprhs=(*rel->y)(2,ngrid);
    cout << "RHS, nprhs=" << nprhs << endl;
    if (nprhs<rhsmin) {
      nprhs=rhsmin;
      cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
    }
  }

  //----------------------------------------------
  // Store solution and delete smnrelax

  for(i=1;i<=ngrid;i++) {
    at->set(0,i-1,(*rel->x)[i]);
    for(j=1;j<=5;j++) {
      at->set(j,i-1,(*rel->y)(j,i));
    }
  }
  delete rel;
  at->set_nlines(ngrid);

  //----------------------------------------------
  // RHS
    
  if (fabs(rp->protfrac-0.5)>1.0e-4) {

    at->set_nlines(ngrid*2-1);

    rel=new smnrelax(3,1,ngrid);
    rel->itmax=100;

    if (rp->iii==1) {
      //----------------------------------------------
      // Construct a simple linear guess
      
      if (npdrip==0.0) {
	for(i=1;i<=ngrid;i++) {
	  (*rel->x)[i]=((double)(i-1))/((double)(ngrid-1));
	  (*rel->y)(1,i)=(nnrhs-nndrip)*(1.0-(*rel->x)[i])+nndrip;
	  (*rel->y)(2,i)=-0.08*(1.0-(*rel->x)[i]);
	  (*rel->y)(3,i)=0.001;
	}
      } else {
	for(i=1;i<=ngrid;i++) {
	  (*rel->x)[i]=((double)(i-1))/((double)(ngrid-1));
	  (*rel->y)(1,i)=(nprhs-npdrip)*(1.0-(*rel->x)[i])+npdrip;
	  (*rel->y)(2,i)=-0.08*(1.0-(*rel->x)[i]);
	  (*rel->y)(3,i)=0.001;
	}
      }
    } else {
      //----------------------------------------------
      // Use last solution for guess
      
      for(i=1;i<=ngrid;i++) {
	(*rel->x)[i]=xg2[i];
	for(j=1;j<=3;j++) {
	  (*rel->y)(j,i)=yg2[j][i];
	}
      }

    }

    //----------------------------------------------
    // Solve
  
    rhsmode=true;
    cout << "second solve" << endl;

    if (nndrip>0.0) {
      cout << "Rhs length: " << rhslength << endl;
      relret2=rel->solve(relaxconverge,1.0);

      while (-(*rel->y)(2,ngrid)>1.0e-4) {
	rhslength*=1.2;

	cout << "Rhs length: " << rhslength << endl;
	relret2=rel->solve(relaxconverge,1.0);
      }
    } else {
      relret2=rel->solve(relaxconverge,1.0);
    }
    
    //----------------------------------------------
    // Copy solution for next guess

    for(i=1;i<=ngrid;i++) {
      xg2[i]=(*rel->x)[i];
      for(j=1;j<=3;j++) {
	yg2[j][i]=(*rel->y)(j,i);
      }
    }
  
    //----------------------------------------------
    // Store solution and delete smnrelax

    for(i=1;i<=ngrid;i++) {
      at->set(0,i+ngrid-2,(*rel->x)[i]);
      for(j=1;j<=3;j++) {
	at->set(j,i+ngrid-2,(*rel->y)(j,i));
      }
    }
    delete rel;

    //----------------------------------------------
    // Rearrangment and rescaling
  
    if (npdrip==0.0) {
      for(i=ngrid-1;i<2*ngrid-1;i++) {
	at->set(5,i,at->get(3,i));
	at->set(4,i,0.0);
	at->set(3,i,at->get(2,i));
	at->set(2,i,0.0);
	at->set(0,i,at->get(0,i)*at->get(5,i)+xrhs);
      }
    } else {
      for(i=1;i<=rel->ngrid;i++) {
	at->set(5,i,at->get(3,i));
	at->set(4,i,at->get(2,i));
	at->set(3,i,0.0);
	at->set(2,i,at->get(1,i));
	at->set(1,i,0.0);
	at->set(0,i,at->get(0,i)*at->get(5,i)+xrhs);
      }
    }
    //    at->set_nlines(2*ngrid-1);
  }

  return 0.0;

}

double monfun2(double lmonfact, int ngrid, void *pa) {
  relaxp *rp=(relaxp *)pa;
  bool guessdone, debug=false;
  double dx, y[5], dydx[5], xrhs;
  int i, j, ilast=0, interpi;
  ofstream itout;
  char ch;
  ubvector sx(3), sy(3);

  monfact=lmonfact;

  //----------------------------------------------
  // Create object

  rel=new smnrelax(3,2,ngrid);
  rel->itmax=200;

  //----------------------------------------------
  // Construct a guess by shooting

  if (rp->iii==1) {
    xstor[1]=0.0;
    ystor(1,1)=rp->nn0+rp->np0;
    ystor(2,1)=firstderiv;
    dx=initialstep/((double)ngrid);
    
    guessdone=false;
    for(i=2;guessdone==false && i<=ngrid;i++) {
      xstor[i]=xstor[i-1]+dx;
      
      sx[1]=(1.0-rp->protfrac)*ystor(1,i-1);
      sx[2]=rp->protfrac*ystor(1,i-1);
      rp->barn=ystor(1,i-1);
      mm_funct qqf(qnnqnpfun);
      rp->nd->msolve(2,sx,qqf);
      
      ystor(1,i)=ystor(1,i-1)+dx*ystor(2,i-1);
      ystor(2,i)=ystor(2,i-1)+dx*(rp->n->mu-rp->mun)/(rp->qnn);
      
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
    
    (*rel->x)[1]=xstor[1];
    (*rel->x)[ngrid]=xstor[ilast];
    for(i=2;i<ngrid;i++) {
      (*rel->x)[i]=(*rel->x)[1]+((double)(i-1))/((double)(ngrid-1))*
	((*rel->x)[ngrid]-(*rel->x)[1]);
    }
    interpi=1;
    for(i=1;i<=ngrid;i++) {
      while((*rel->x)[i]>xstor[interpi+1] && interpi<ilast-1) interpi++;
      while((*rel->x)[i]<xstor[interpi] && interpi>1) interpi--;
      for(j=1;j<=2;j++) {
	(*rel->y)(j,i)=ystor(j,interpi)+
	  (ystor(j,interpi+1)-ystor(j,interpi))*
	  ((*rel->x)[i]-xstor[interpi])/(xstor[interpi+1]-xstor[interpi]);
      }
      (*rel->y)(3,i)=((*rel->x)[ngrid]);
      (*rel->x)[i]=((double)(i-1))/((double)(ngrid-1));
    }
  } else {
    //----------------------------------------------
    // Use last solution for guess
    
    for(i=1;i<=ngrid;i++) {
      (*rel->x)[i]=xg1[i];
      for(j=1;j<=rel->ne;j++) {
	(*rel->y)(j,i)=yg1[j][i];
      }
    }
  }

  //----------------------------------------------
  // Solve

  rhsmode=false;
  relret1=rel->solve(relaxconverge,1.0);
  
  //----------------------------------------------
  // Copy solution for next guess

  for(i=1;i<=ngrid;i++) {
    xg1[i]=(*rel->x)[i];
    for(j=1;j<=rel->ne;j++) {
      yg1[j][i]=(*rel->y)(j,i);
    }
  }
  
  //----------------------------------------------
  // Rescale solution and set neutron and proton
  // densities:

  at->set_nlines(ngrid);
  for(i=1;i<=rel->ngrid;i++) {

    sx[1]=(1.0-rp->protfrac)*(*rel->y)(1,i);
    sx[2]=rp->protfrac*(*rel->y)(1,i);

    rp->barn=(*rel->y)(1,i);
    mm_funct qqf(qnnqnpfun);
    rp->nd->msolve(2,sx,qqf);

    at->set(0,i-1,(*rel->x)[i]*(*rel->y)(3,i));
    at->set(1,i-1,sx[1]);
    at->set(2,i-1,sx[2]);
    at->set(3,i-1,(*rel->y)(2,i));
    at->set(4,i-1,0.0);
    at->set(5,i-1,(*rel->y)(3,i));
  }

  //----------------------------------------------
  // Store quantites at rhs for later use:
  
  xrhs=(*at)[0][ngrid-1];
  if (npdrip==0.0) {
    nnrhs=(*at)[1][ngrid-1];
    cout << "RHS, nnrhs=" << xrhs << " " << nnrhs << endl;
    if (nnrhs<rhsmin) {
      nnrhs=rhsmin;
      cout << "Adjusting nnrhs to nnrhs=" << nnrhs << endl;
    }
  } else {
    nprhs=(*at)[2][ngrid-1];
    cout << "RHS, nprhs=" << nprhs << endl;
    if (nprhs<rhsmin) {
      nprhs=rhsmin;
      cout << "Adjusting nprhs to nprhs=" << nprhs << endl;
    }
  }

  //----------------------------------------------
  // Delete smnrelax

  delete rel;

  //---------------------------------------------
  // RHS
    
  if (fabs(rp->protfrac-0.5)>1.0e-4) {
    rel=new smnrelax(3,1,ngrid);
    
    if (rp->iii==1) {
      //----------------------------------------------
      // Construct a simple linear guess
      
      if (npdrip==0.0) {
	for(i=1;i<=ngrid;i++) {
	  (*rel->x)[i]=((double)(i-1))/((double)(ngrid-1));
	  (*rel->y)(1,i)=(nnrhs-nndrip)*(1.0-(*rel->x)[i])+nndrip;
	  (*rel->y)(2,i)=-0.08*(1.0-(*rel->x)[i]);
	  (*rel->y)(3,i)=0.001;
	}
      } else {
	for(i=1;i<=ngrid;i++) {
	  (*rel->x)[i]=((double)(i-1))/((double)(ngrid-1));
	  (*rel->y)(1,i)=(nprhs-npdrip)*(1.0-(*rel->x)[i])+npdrip;
	  (*rel->y)(2,i)=-0.08*(1.0-(*rel->x)[i]);
	  (*rel->y)(3,i)=0.001;
	}
      }
    } else {
      //----------------------------------------------
      // Use last solution for guess
      
      for(i=1;i<=ngrid;i++) {
	(*rel->x)[i]=xg2[i];
	for(j=1;j<=rel->ne;j++) {
	  (*rel->y)(j,i)=yg2[j][i];
	}
      }
      
    }
    
    //----------------------------------------------
    // Solve
    
    rhsmode=true;
    relret2=rel->solve(relaxconverge,1.0);
    
    //----------------------------------------------
    // Copy solution for next guess
    
    for(i=1;i<=ngrid;i++) {
      xg2[i]=(*rel->x)[i];
      for(j=1;j<=rel->ne;j++) {
	yg2[j][i]=(*rel->y)(j,i);
      }
    }
    
    //----------------------------------------------
    // Store solution, rescale, and delete smnrelax
    
    at->set_nlines(2*ngrid-1);
    for(i=1;i<=rel->ngrid;i++) {
      at->set(5,i-2+ngrid,(*rel->y)(3,i));
      at->set(3,i-2+ngrid,(*rel->y)(2,i));
      at->set(4,i-2+ngrid,0.0);
      at->set(2,i-2+ngrid,0.0);
      at->set(1,i-2+ngrid,(*rel->y)(1,i));
      at->set(0,i-2+ngrid,(*rel->x)[i]*(*rel->y)(3,i)+xrhs);
    }
    delete rel;
  }
  
  return 0.0;

}

int qnnqnpfun(size_t sn, const ubvector &sx, ubvector &sy) {
  relaxp *rp;//=(relaxp *)pa;

  rp->n->n=sx[1];
  rp->p->n=sx[2];
  if (sx[1]<0.0 || sx[2]<0.0) return 1;
  
  rp->eos->calc_e(*rp->n,*rp->p,*rp->hb);

  sy[1]=rp->n->n+rp->p->n-rp->barn;
  sy[2]=rp->n->mu-rp->p->mu-rp->mun+rp->mup;

  return 0;
}

int ndripfun(size_t sn, const ubvector &sx, ubvector &sy) {
  relaxp *rp;//=(relaxp *)pa;
  double pleft, pright, munleft, munright;

  if (sx[1]<0.0 || sx[2]<0.0 || sx[3]<0.0) return 1;

  rp->n->n=sx[1];
  rp->p->n=sx[2];
  sy[1]=rp->p->n-rp->protfrac*(rp->p->n+rp->n->n);
  
  rp->eos->calc_e(*rp->n,*rp->p,*rp->hb);
  pleft=rp->hb->pr;
  munleft=rp->n->mu;

  rp->n->n=sx[3];
  rp->p->n=0.0;
  
  rp->eos->calc_e(*rp->n,*rp->p,*rp->hb);
  pright=rp->hb->pr;
  munright=rp->n->mu;

  sy[2]=pleft-pright;
  sy[3]=munleft-munright;

  return 0;
}

int pdripfun(size_t sn, const ubvector &sx, ubvector &sy) {
  relaxp *rp;//=(relaxp *)pa;
  double pleft, pright, mupleft, mupright;

  if (sx[1]<0.0 || sx[2]<0.0 || sx[3]<0.0) return 1;

  rp->n->n=sx[1];
  rp->p->n=sx[2];
  sy[1]=rp->p->n-rp->protfrac*(rp->p->n+rp->n->n);
  
  rp->eos->calc_e(*rp->n,*rp->p,*rp->hb);
  pleft=rp->hb->pr;
  mupleft=rp->p->mu;

  rp->n->n=0.0;
  rp->p->n=sx[3];
  
  rp->eos->calc_e(*rp->n,*rp->p,*rp->hb);
  pright=rp->hb->pr;
  mupright=rp->p->mu;

  sy[2]=pleft-pright;
  sy[3]=mupleft-mupright;

  return 0;
}

