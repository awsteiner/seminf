#include "relax.h"

using namespace std;

void relax::regrid(int new_ne, int new_nb, int new_grid) {
  int lo, hi, i, j;
  ubvector x2;
  ubmatrix y2;

  // Make the new space
  x2.resize(new_grid+1);
  y2.resize(new_ne+1,new_grid+1);

  // Set up boundaries properly
  x2[1]=x[1];
  x2[new_grid]=x[ngrid];
  for(j=1;j<=new_ne;j++) {
    if (j<=ne) {
      y2(j,1)=y(j,1);
      y2(j,new_grid)=y(j,ngrid);
    } else {
      y2(j,1)=0.0;
      y2(j,new_grid)=0.0;
    }
  }

  // Fill in the remaining
  for(i=2;i<=new_grid-1;i++) {
    x2[i]=x2[1]+((double)(i-1))/((double)(new_grid-1))*
      (x2[new_grid]-x2[1]);
  }
  lo=1;
  hi=2;
  for(i=2;i<=new_grid-1;i++) {
    while(x2[i]>x[hi]) {
      lo++;
      hi++;
    }
    for(j=1;j<=new_ne;j++) {
      if (j<=ne) {
	y2(j,i)=y(j,lo)+(x2[i]-x[lo])*(y(j,hi)-y(j,lo))/
	  (x[hi]-x[lo]);
      } else {
	y2(j,i)=0.0;
      }
    }
  }
  
  x=x2;
  y=y2;
  ngrid=new_grid;
  ne=new_ne;
  nb=new_nb;

  // Make space for a new s and c:
  s.resize(ne+1,2*ne+2);
  c.resize(ne+1);
  for(int ic=0;ic<ne+1;ic++) c[ic].resize(ne-nb+2,ngrid+2);
  
  return;
}

relax::relax(int tne, int tnb, int tngrid) {
  eps=1.0e-9; // Is this too small?
  ne=tne;
  nb=tnb;
  ngrid=tngrid;
  itmax=40;
  indexv.resize(ne+1);
  scalv.resize(ne+1);
  for(int i=1;i<=ne;i++) {
    indexv[i]=i;
    scalv[i]=1.0;
  }
  x.resize(ngrid+1);
  y.resize(ne+1,ngrid+1);
  s.resize(ne+1,2*ne+2);
  c.resize(ne+1);
  for(int ic=0;ic<ne+1;ic++) c[ic].resize(ne-nb+2,ngrid+2);
}

int relax::solve(double conv, double slowc) {
  int ic1,ic2,ic3,ic4,it,j,jj1,j2,j3,j4,j5,j6,j7,j8,j9;
  int jc1,jcf,jv,k,k1,k2,km,kp,nvars,ret;
  ubvector_int kmax;
  double err,errj,fac,vmax,vz;
  ubvector ermax;

  kmax.resize(ne+1);
  ermax.resize(ne+1);
  
  k1=1;
  k2=ngrid;
  nvars=ne*ngrid;
  jj1=1;
  j2=nb;
  j3=nb+1;
  j4=ne;
  j5=j4+jj1;
  j6=j4+j2;
  j7=j4+j3;
  j8=j4+j4;
  j9=j8+jj1;
  ic1=1;
  ic2=ne-nb;
  ic3=ic2+1;
  ic4=ne;
  jc1=1;
  jcf=ic3;
  for (it=1;it<=itmax;it++) {
    k=k1;
    
    ret=difeq(k,k1,k2,j9,ic3,ic4);
    if (ret!=0) {
      return ret;
    }
    ret=calcderiv(k,k1,k2,j9,ic1,ic4);
    if (ret!=0) {
      return ret;
    }
    ret=pinvs(ic3,ic4,j5,j9,jc1,k1);
    if (ret!=0) {
      return ret;
    }
    
    for (k=k1+1;k<=k2;k++) {

      kp=k-1;
      ret=difeq(k,k1,k2,j9,ic1,ic4);
      if (ret!=0) {
	return ret;
      }
      ret=calcderiv(k,k1,k2,j9,ic1,ic4);
      if (ret!=0) {
	return ret;
      }
      red(ic1,ic4,jj1,j2,j3,j4,j9,ic3,jc1,jcf,kp);
      ret=pinvs(ic1,ic4,j3,j9,jc1,k);
      if (ret!=0) {
	return ret;
      }
    }
    
    k=k2+1;
    ret=difeq(k,k1,k2,j9,ic1,ic2);
    if (ret!=0) {
      return ret;
    }
    ret=calcderiv(k,k1,k2,j9,ic1,ic4);
    if (ret!=0) {
      return ret;
    }
    red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2);
    ret=pinvs(ic1,ic2,j7,j9,jcf,k2+1);
    if (ret!=0) {
      return ret;
    }

    bksub(jcf,k1,k2);

    err=0.0;
    for (j=1;j<=ne;j++) {
      jv=indexv[j];
      errj=vmax=0.0;
      km=0;
      for (k=k1;k<=k2;k++) {
	vz=fabs(c[jv](1,k));
	if (vz > vmax) {
	  vmax=vz;
	  km=k;
	}
	errj += vz;
      }
      err += errj/scalv[j];
      ermax[j]=c[jv](1,km)/scalv[j];
      kmax[j]=km;
    }
    err /= nvars;
    fac=(err > slowc ? slowc/err : 1.0);
    fac=1.0;
    for (j=1;j<=ne;j++) {
      jv=indexv[j];
      for (k=k1;k<=k2;k++) {
	y(j,k) -= fac*(c[jv](1,k));
      }
    }

    ret=iter(it,err,fac,kmax,ermax);
    if (ret!=0) {
      return ret;
    }
    if (err<conv) {
      return 0;
    }
    
  }
  O2SCL_ERR("Exceeded maximum number of iterations in solve().",
	    o2scl::exc_efailed);
  return 10;
}

int relax::calcderiv(int k, int k1, int k2, int jsf, int is1, int isf) {
  double f,f2,h,tempd;
  int i,j,ret;
  
  if (k==k1) { 
    for(j=ne-nb+1;j<=ne;j++) {
      
      // go back to the old variables and calculate s[i][jsf]
      // for each run
      ret=difeq(k,k1,k2,jsf,is1,isf);
      if (ret!=0) return ret;
      f2=s(j,jsf);

      for(i=1;i<=ne;i++) {
	tempd=y(i,k);
	h=eps*fabs(tempd);
	if (fabs(h)<1.e-15) h=eps;
	y(i,k)=tempd+h;
	h=y(i,k)-tempd;
	ret=difeq(k,k1,k2,jsf,is1,isf);
	if (ret!=0) return ret;
	f=s(j,jsf);
	y(i,k)=tempd;

	s(j,ne+indexv[i])=(f-f2)/h;
	
      }
    }

  } else if (k>k2) {
    for(j=1;j<=ne-nb;j++) {
      
      // go back to the old variables and calculate s[i][jsf]
      // for each run
      ret=difeq(k,k1,k2,jsf,is1,isf);
      if (ret!=0) return ret;
      f2=s(j,jsf);
      
      for(i=1;i<=ne;i++) {
	tempd=y(i,ngrid);
	h=eps*fabs(tempd);
	if (fabs(h)<1.e-15) h=eps;
	y(i,ngrid)=tempd+h;
	h=y(i,ngrid)-tempd;
	ret=difeq(k,k1,k2,jsf,is1,isf);
	if (ret!=0) return ret;
	f=s(j,jsf);
	y(i,ngrid)=tempd;

	s(j,ne+indexv[i])=(f-f2)/h;
      }
    }
  } else {
    for(j=1;j<=ne;j++) {
      
      // go back to the old variables and calculate s[i][jsf]
      // for each run
      
      ret=difeq(k,k1,k2,jsf,is1,isf);
      if (ret!=0) return ret;
      f2=s(j,jsf);
      
      for(i=1;i<=ne;i++) {
	tempd=y(i,k-1);
	h=eps*fabs(tempd);
	if (fabs(h)<1.e-15) h=eps;
	y(i,k-1)=tempd+h;
	h=y(i,k-1)-tempd;
	ret=difeq(k,k1,k2,jsf,is1,isf);
	if (ret!=0) return ret;
	f=s(j,jsf);
	y(i,k-1)=tempd;
	s(j,indexv[i])=(f-f2)/h;
      
	tempd=y(i,k);
	h=eps*fabs(tempd);
	if (fabs(h)<1.e-15) h=eps;
	y(i,k)=tempd+h;
	h=y(i,k)-tempd;
	ret=difeq(k,k1,k2,jsf,is1,isf);
	if (ret!=0) return ret;
	f=s(j,jsf);
	y(i,k)=tempd;
	s(j,ne+indexv[i])=(f-f2)/h;
      }
    }
  }

  return 0;
}

void relax::bksub(int jf, int k1, int k2) {
  int nbf,im,kp,k,j,i;
  double xx;

  nbf=ne-nb;
  im=1;
  for (k=k2;k>=k1;k--) {
    if (k == k1) im=nbf+1;
    kp=k+1;
    for (j=1;j<=nbf;j++) {
      xx=c[j](jf,kp);
      for (i=im;i<=ne;i++) {
	c[i](jf,k)-=c[i](j,k)*xx;
      }
    }
  }
  for (k=k1;k<=k2;k++) {
    kp=k+1;
    for (i=1;i<=nb;i++) {
      c[i](1,k)=c[i+nbf](jf,k);
    }
    for (i=1;i<=nbf;i++) {
      c[i+nb](1,k)=c[i](jf,kp);
    }
  }
  return;
}

void relax::red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, 
		   int jmf, int ic1, int jc1, int jcf, int kc) {
  int loff,l,j,ic,i;
  double vx;

  loff=jc1-jm1;
  ic=ic1;
  for (j=jz1;j<=jz2;j++) {
    for (l=jm1;l<=jm2;l++) {
      vx=c[ic](l+loff,kc);
      for (i=iz1;i<=iz2;i++) s(i,l) -= s(i,j)*vx;
    }
    vx=c[ic](jcf,kc);
    for (i=iz1;i<=iz2;i++) s(i,jmf) -= s(i,j)*vx;
    ic += 1;
  }

  return;
}

int relax::pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k) {
  int js1,jpiv=0,jp=0,je2,jcoff,j,irow,ipiv=0,id,icoff,i;
  double pivinv,piv,dum,big;
  ubvector_int indxr(1+ie1+ie2);
  ubvector pscl(1+ie1+ie2);

  je2=je1+ie2-ie1;
  js1=je2+1;
  for (i=ie1;i<=ie2;i++) {
    big=0.0;
    for (j=je1;j<=je2;j++)
      if (fabs(s(i,j)) > big) big=fabs(s(i,j));
    if (big == 0.0) {
      O2SCL_ERR("Singular matrix in pinvs().",1);
      return 1;
    }
    pscl[i]=1.0/big;
    indxr[i]=0;
  }
  for (id=ie1;id<=ie2;id++) {
    piv=0.0;
    for (i=ie1;i<=ie2;i++) {
      if (indxr[i] == 0) {
	big=0.0;
	for (j=je1;j<=je2;j++) {
	  if (fabs(s(i,j)) > big) {
	    jp=j;
	    big=fabs(s(i,j));
	  }
	}
	if (big*pscl[i] > piv) {
	  ipiv=i;
	  jpiv=jp;
	  piv=big*pscl[i];
	}
      }
    }
    if (s(ipiv,jpiv) == 0.0) {
      O2SCL_ERR("Singular matrix in pinvs() (2).",2);
      return 2;
    }
    indxr[ipiv]=jpiv;
    pivinv=1.0/s(ipiv,jpiv);
    for (j=je1;j<=jsf;j++) s(ipiv,j) *= pivinv;
    s(ipiv,jpiv)=1.0;
    for (i=ie1;i<=ie2;i++) {
      if (indxr[i] != jpiv) {
	if (s(i,jpiv)) {
	  dum=s(i,jpiv);
	  for (j=je1;j<=jsf;j++) {
	    s(i,j) -= dum*s(ipiv,j);
	  }
	  s(i,jpiv)=0.0;
	}
      }
    }
  }
  jcoff=jc1-js1;
  icoff=ie1-je1;
  for (i=ie1;i<=ie2;i++) {
    irow=indxr[i]+icoff;
    if (irow<0) {
      O2SCL_ERR("irow<0 in pinvs().",3);
      return 3;
    } else {
      for (j=js1;j<=jsf;j++) {
	c[irow](j+jcoff,k)=s(i,j);
      }
    }
  }
  
  return 0;
}

int relax::iter(int k, double err, double fac, ubvector_int &kmax,
		ubvector &ermax) {
  
  if (k==1) {
    cout << "\tIter.\tError\t\tFac." << endl;
  }
  if (k>=1) {
    cout << "\t" << k << "  \t" << err << "\t" << fac << endl;
  }
  
  if (k>=itmax) {
    return 10;
  }
  
  return 0;

}

