/****************** planet-system ******************/

using namespace std;
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

/*************  important parameters  **************/
const int neqn = 4;           // number of equations
int schalter = 0;
bool abb=false;
const double tstep=0.001;
const double gen=0.001; //Genauigkeit
const double mu=0.05;
const double mu2=1-mu;
const double r1=0.2;
const double r2=0.01;
double omeg1,omeg2;
double z1[3], z2;	//fuer das Bisektionsverfahren
const double pi=4.*atan(1.);
//const double alpha=0.0995*pi;

/******************************************************/
int sgn(double h){
    if (h==0)	return 0;
    else		return (h>0) ? 1 : -1;
}
/******************** Runge-Kutta solver **********************************/
void rk(double y[],int n,double x,double h,
     void (derivs)(double, double[], double[]))
{
  	int j,i;
  	double h2,h6,xh, xhh,y1[n],k1[n],k2[n],k3[n];
  	h2=h*0.5;    h6=h/6.;
	//Simpsonregel
    xh=x+h2;	xhh=x+h;
    (derivs)(x,y,k1);
    for(i=0;i<n;i++) y1[i]=y[i]+h2*k1[i];
	(derivs)(xh,y1,k2);
	for(i=0;i<n;i++) y1[i]=y[i]-h*k1[i]+2*h*k2[i];
    (derivs)(xhh,y1,k3);
    x+=h;
    for(i=0;i<n;i++) y[i]+=h6*(k1[i]+4.*k2[i]+k3[i]);
}

void planet(double t,double x[],double dxdt[])
{
 double l1, l2, k1, k2;//x1=x[0], x2=x[1], x3=x[2], x4=x[3]
 k1=x[0]+mu;
 k2=k1-1.;
 l1=k1*k1+x[1]*x[1];
 l2=k2*k2+x[1]*x[1];
 dxdt[0] =x[2];
 dxdt[1] =x[3];
 dxdt[2] =2*x[3]+x[0]-(1-mu)*k1/pow(l1, 1.5)-mu*k2/pow(l2, 1.5);
 dxdt[3] =-2*x[2]+x[1]-(1-mu)*x[1]/pow(l1, 1.5)-mu*x[1]/pow(l2, 1.5);
}

/******************************************************/

double f(double y_0){
	double k1, k2,g;
 	k1=y_0+mu;
 	k2=k1-1.;
 	g=y_0-(1-mu)*sgn(k1)/(k1*k1)-mu*sgn(k2)/(k2*k2);
 	return g;
}


/******************************************************/

double Omega(double a){
	double u;
	u=a*a/2.+(1-mu)/fabs(a+mu)+mu/fabs(a-1.+mu)+mu*(1-mu)/2.;
	return u;
}
/*******************Bisection-Method********************/
/*finds zeros of equations*/
void bisec(int a, int m, double b, int i){
	if(a==m)z1[i]=b;
	else 	z2=b;
}


double lag(int j){
	double w,c;
	int m, g;
	w=fabs(z2-z1[j]);
	m=sgn(f(z1[j]));
  	while(w>gen){
    	c=(z2+z1[j])/2.;
    	g=sgn(f(c));
    	bisec(g,m, c,j);
    	w=fabs(z2-z1[j]);
	}
	return c;
}
/******************** Symplectic Integration **********************************/
void Hamilton(double t, double x[],double F[], bool Imp){
  	double l1, l2, k1, k2;//x1=x[0], x2=x[1], x3=x[2], x4=x[3]
  	if(Imp){
  		k1=x[0]+mu;
  		k2=k1-1.;
  		l1=k1*k1+x[1]*x[1];
  		l2=k2*k2+x[1]*x[1];
  		F[2]=-1*(x[0]-(1-mu)*k1/pow(l1, 1.5)-mu*k2/pow(l2, 1.5));
  		F[3]=-1*(x[1]-(1-mu)*x[1]/pow(l1, 1.5)-mu*x[1]/pow(l2, 1.5));
  	}
	else{
  		F[0]=x[2]+x[1];
  		F[1]=x[3]-x[0];
  	}
}

void Sympl_Integr(double x[],int n,double &t,double h,
     void (Ham_deriv)(double, double[], double[], bool))
{
  int i,k=int(n/2);
  double L[n], h2=h/2.;
  (Ham_deriv)(t,x,L,true);
  for(i=k;i<n;i++){//Halbschritt der Impulse
    x[i]+=-h2*L[i];
  }
  (Ham_deriv)(t,x,L,false);
  for(i=0;i<k;i++){//Schritt der Orte
    x[i]+=h*L[i];
  }
  (Ham_deriv)(t,x,L,true);
  for(i=k;i<n;i++){//Halbschritt der Impulse
    x[i]+=-h2*L[i];
  }
  t+=h;
}

/*********************Abbruchbedingungen***************/

double cond_kreis(double x,double a, double r){
	double abs_quad;
	abs_quad=r*r-(x-a)*(x-a);
	return abs_quad;
}

void cond(double x[]){
	double g;
	if (fabs(x[0])>2)abb=true;
	if (fabs(x[1])>2)abb=true;
	if(fabs(x[0]+mu)<r1){
		g=x[1]*x[1];
		if(g<cond_kreis(x[0],-1*mu,r1))abb=true;
	}
	if(fabs(x[0]-mu2)<r2){
		g=x[1]*x[1];
		if(g<cond_kreis(x[0],mu2,r2))abb=true;
	}
}

int  count(double x[],double l, int cnt){
	if(x[0]>mu2){
		if(sgn(l)!=sgn(x[1]))	cnt+=1;
	}
	return cnt;
}

double opt_wnk(double &wnk, double &wnk_step, int cnt, int &cnt1, int &j, int &i, int &max){
	if (cnt>cnt1){
		j=1;
		cnt1=cnt;
		wnk+=wnk_step;
		return j;
	}
	else{
		if (cnt==cnt1){
			if(cnt==0){
				wnk+=wnk_step;
				return 0;
			}
			else{
				if(j==max-1){
					wnk-=double(j)*wnk_step;
					max+=max;
					wnk_step=wnk_step/2.;
					i=max;
					return 0;
				}
				j+=1;
				wnk+=wnk_step;
				return j;
			}
		}
		else{
			if(j==1){
				wnk-=2*wnk_step;
				wnk_step=2*wnk_step/double(max);
			}
			else{
				wnk-=double(j)*wnk_step;
				wnk_step=wnk_step/double(max-j);
			}
			j=0;
			i=max;
			return j;
		}
	}
}

double energy(double x[]){
  return (x[0]*x[0]+x[1]*x[1])/2.+(1-mu)/sqrt((x[0]+mu)*(x[0]+mu)+x[1]*x[1])+mu/sqrt((x[0]-1.+mu)*(x[0]-1.+mu)+x[1]*x[1])+mu*(1-mu)/2.-(x[2]*x[2]+x[3]*x[3])/2.;
}

main( int argc, char *argv[] ){
 	int 
    cnt=0, 
    cnt1=0,
    k=0,
    max=10
    ;
 	bool 
    abb1=false
    ;
 	double 
    x[neqn], 
    t, 
    lag_punkt[3], 
    wnk_step=(pi/2.)/double(max), 
    wnk =0.38,
    w,
    l,
    m=0
    ;
  const char* 
    Dateiname
    ;
  schalter = atoi(argv[1]);
	z1[0]=-1.9;
  z1[1]=mu*(-1)+gen;
  z1[2]=1-mu+gen;
  if (schalter==1){
    Dateiname="DaT_Dateien//planet4.dat";
  }
  else{
    Dateiname="DaT_Dateien//planet5.dat";
  } 
  cout << Dateiname << endl;
  ofstream fout(Dateiname,ios::out); // file for output
  fout.setf(ios::scientific);fout.precision(6);
  	x[1]=0.; x[2]=0.; x[3]=0.; t=0.;  //x1=x[0], x2=x[1], x1^punkt=x[2], x2^punkt=x[3] Initial conditions
  for(int i=0;i<3;i++){
    z2=z1[i]+mu2;
    lag_punkt[i]=lag(i);
  }
	x[0]=lag_punkt[1];
  omeg1=Omega(x[0]);
	x[0]=0.15;	//Startpunkt fuer den Trek
	omeg2=Omega(x[0]);
	w=(sqrt(2*(omeg2-omeg1)))+0.078;
  ////////// Startparamter Ausgabe
  fout<<"#" << " x_1=x="<< x[0]<< ", x_2=y="<< x[1]<< ", x_3=x^dot="<< x[2]<< ", x_4=y^dot="<< x[3]<< ", Winkel="<< wnk<<"rad "<<", |x_vec|=" << w <<", schalter="<<schalter<<", tstep="<<tstep<< endl;
  fout<<"#" << " x_1=x     "<< "|     x_2=y  " <<"| energy      " << endl;
  //////////
	while(abb1==false){
		for(int i=0; i<max;i++){
			x[2]=cos(wnk)*w;  //Startgeschwindigkeit (x-Richtung) fuer den Trek
			x[3]=sin(wnk)*w;  //Startgeschwindigkeit (y-Richtung) fuer den Trek
      if (schalter==1){ //Symplectic Integration
        while(abb==false){
          l=x[1];
          x[2]=x[2]-x[1];
          x[3]=x[3]+x[0];
          Sympl_Integr(x,neqn,t,tstep,Hamilton);
          if (wnk_step<0.0001) fout<<x[0]<<" "<<x[1]<<" "<<energy(x)<<endl;
          x[2]=x[2]+x[1];
          x[3]=x[3]-x[0];
          cond(x);
          cnt=count(x,l,cnt);
          if (wnk_step<0.0001) cout<<cnt<<endl;
        }
      }
      else{ //Runge-Kutta Verfahren
        while(abb==false){
          l=x[1];
          rk(x,neqn,t,tstep,planet);
          if (wnk_step<0.0001) fout<<x[0]<<" "<<x[1]<<" "<<energy(x)<<endl;
          cond(x);
          count(x,l,cnt);
        }
      }
      if (wnk_step<0.0001){
        abb1=true;
        i=max;
      }
      k=opt_wnk(wnk, wnk_step, cnt, cnt1,k,i, max);
      x[0]=0.15; 
      x[1]=0.;
      abb=false; 
      cnt=0;
    }
    cout<<wnk_step<<endl;
    cnt1=0;
  }
  fout.close();
}
