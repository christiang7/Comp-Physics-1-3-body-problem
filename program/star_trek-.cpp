/****************** planet-system ******************/

using namespace std;
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>


/*************  important parameters  **************/
const int neqn = 4;           // number of equations
bool schalter = true;
bool abb=false;
bool abb4=false;
const double tstep=0.001;
const double gen=0.001; //Genauigkeit der Lagrangepunkte
const double eps=0.00001;
const double mu=0.05;
const double mu2=1-mu;
const double r1=0.2;
const double r2=0.01;
double omeg1,omeg2,omeg3;
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
 double l1, l2, k1, k2;//x1=y[0], x2=y[1], x3=y[2], x4=y[3]
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

double Omega(double a, double b){
	double u;
	u=(a*a+b*b)/2.+(1-mu)/sqrt((a+mu)*(a+mu)+b*b)+mu/sqrt((a-1.+mu)*(a-1.+mu)+b*b)+mu*(1-mu)/2.;
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
  	double l1, l2, k1, k2;//x1=y[0], x2=y[1], x3=y[2], x4=y[3]
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
	if (fabs(x[0])>2){
		abb=true;
		abb4=true;
	}
	if (fabs(x[1])>2){
		abb=true;
		abb4=true;
	}
	if(fabs(x[0]+mu)<r1){
		g=x[1]*x[1];
		if(g<cond_kreis(x[0],-1*mu,r1)){
			abb=true;
			//abb4=true;
		}
	}
	if(fabs(x[0]-mu2)<r2){
		g=x[1]*x[1];
		if(g<cond_kreis(x[0],mu2,r2))abb=true;
	}
}

int  count(double x[],double l, int cnt){
	if(x[1]>0){
		if(sgn(l-mu2)!=sgn(x[0]-mu2))	cnt+=1;
	}
	return cnt;
}

int  count1(double x[],double l, int cnt){
	if(x[1]<0){
		if(sgn(l+mu)!=sgn(x[0]+mu))	cnt+=1;
	}
	return cnt;
}

void opt_wnk(int j, double &wnk_step, double &wnk, double wnk_vec[], int i, int m){
	if(j==m)	wnk_step=eps;
	else{
		if(j==1){
			if (i==0){
				wnk=wnk_vec[i-1];
				wnk_step=wnk_step/(double(m));
			}
			if(i==m-1){
				wnk=wnk_vec[i-1];
				wnk_step=wnk_step/(double(m));
			}
			else{
				wnk=wnk_vec[i-1];
				wnk_step=2*wnk_step/(double(m));
			}
		}
		else{
			wnk=wnk_vec[i];
			wnk_step=double(j-1)*wnk_step/(double(m));
		}
	}
}

int list_max(int list1[], double list2[], int m){
	int y=0;
  	for(int i=0;i<m;i++)	if(list1[i]>y)	y=list1[i];
	return y;//2015
}

int list_koord(int list1[], double list2[], int m, int y){
  	for(int i=0;i<m;i++){
		if(list1[i]==y)	return i;
	}
}

int list_len(int list1[], double list2[], int m, int y){
	int j=0;
  	for(int i=0;i<m;i++){
		if(list1[i]==y)	j+=1;
	}
	return j;
}

double energy(double x[]){
  return (x[0]*x[0]+x[1]*x[1])/2.+(1-mu)/sqrt((x[0]+mu)*(x[0]+mu)+x[1]*x[1])+mu/sqrt((x[0]-1.+mu)*(x[0]-1.+mu)+x[1]*x[1])+mu*(1-mu)/2.-(x[2]*x[2]+x[3]*x[3])/2.;
}

main(){
 	int cnt=0, cnt1=0,k=0,max=0,m,m1=0, m2,len,len2, koord,koord2; //max =0; hier dazu ändern
 	int cnt_vec[max],cnt_vec2[max];
 	bool abb1=false,abb2=false;
 	double  t=0, wnk_step=(pi/2.)/double(max), wnk=0,w,w1,l,wnk_max, wnkwnk, ww;
 	double x[neqn] ,wnk_vec[max],w_vec[max],lag_punkt[3];
  	const char* Dateiname;
	z1[0]=-1.9;z1[1]=mu*(-1)+gen;z1[2]=1-mu+gen;
  	if (schalter)Dateiname="DaT_Dateien//planet4.dat";
  	else Dateiname="DaT_Dateien//planet5.dat";
 	 cout << Dateiname << endl;
  	x[1]=0.; x[2]=0.; x[3]=0.; t=0.;  //x1=x[0], x2=x[1], x1^punkt=x[2], x2^punkt=x[3] Initial conditions
  	for(int i=0;i<3;i++){
    	z2=z1[i]+mu2;
    	lag_punkt[i]=lag(i);
  	}
    omeg1=Omega(lag_punkt[2],0);
    omeg2=Omega(lag_punkt[0],0);
	x[0]=0.15;//gen;	//Startpunkt fuer den Trek
	omeg3=Omega(x[0],x[1]);
	w=(sqrt(2*(omeg3-omeg1)));
	w1=(sqrt(2*(omeg3-omeg2))+0.1);
	cout<<w<<" "<<w1<<endl;
	double w_step=(w1-w)/double(max);
  //abb2=false;
	while(abb2==false){
		for(int n=0; n<max;n++){
			while(abb1==false){
				for(int i=0; i<max;i++){
					x[2]=cos(wnk)*w-x[1];		//Startgeschwindigkeit (x-Richtung) fuer den Trek
					x[3]=sin(wnk)*w+x[0];		//Startgeschwindigkeit (y-Richtung) fuer den Trek
					if (schalter){ //Symplectic Integration
						while(abb==false){
							l=x[0];
							Sympl_Integr(x,neqn,t,tstep,Hamilton);
				  			cond(x);
						  	cnt=count(x,l,cnt);
						  	//cnt1=count1(x,l,cnt1);
			  			}
	  				}
					else{ 			//Runge-Kutta Verfahren
						while(abb==false){
							l=x[0];
					   		rk(x,neqn,t,tstep,planet);
			    			cond(x);
						    cnt=count(x,l,cnt);
					   }
	  				}
	  				//cout<<cnt<<" "<<n<<endl;
  					cnt_vec[i]=cnt;
  					wnk_vec[i]=wnk;
			  		if(abb4==true){
			  			if (cnt==4){
			  				//if(cnt1>=cnt2){
			  					cout<<cnt<<" "<<cnt1<<" "<<w<<" "<<wnk<<endl;
				  				ww=w;
				  				wnkwnk=wnk;
				  				n=max;
				  				abb2=true;
			  					abb1=true;
			  					i=max;
			  					//cnt1=cnt;
			  			//	}
			  			//	else	abb4=false;
			  			}
			  			else	abb4=false;
			  		}
			  		wnk+=wnk_step;
					  x[0]=0.15; x[1]=0.;abb=false;cnt=0;t=0;//cnt1=0;
	  				if(wnk_step<=eps){
						abb1=true;
						i=max;
					}
 	 			}
		  	m=list_max(cnt_vec,wnk_vec, max);
 			 	koord=list_koord(cnt_vec,wnk_vec, max, m);
  			len=list_len(cnt_vec,wnk_vec, max,m);
				opt_wnk(len,wnk_step,wnk,wnk_vec, koord, max);
			}
			if (m>m1){
				wnk_max=wnk_vec[koord];
				m1=m;
			}
			cnt_vec2[n]=m;
  		w_vec[n]=w;
  		w+=w_step;
  		wnk_step=(pi/2.)/double(max); wnk=0;abb1=false;
		}
		m2=list_max(cnt_vec2,w_vec, max);
 	 	koord2=list_koord(cnt_vec2,w_vec, max, m2);
  	len2=list_len(cnt_vec2,w_vec, max,m2);
		opt_wnk(len2,w_step,w,w_vec, koord2, max);
		//cout<<m2<<" "<<w_step<<" "<<w<<endl;
		m1=0;
		if(w_step<=eps){
			abb2=true;
		}
	}
	ofstream fout(Dateiname,ios::out); // file for output
	fout.setf(ios::scientific);fout.precision(6);
	w=(sqrt(2*(Omega(0.15,0)-Omega(lag_punkt[1],0) )))+0.06;//w_vec[koord2]; hier ändern von ww zu
	cout<<m2<<" "<<w<<" "<<wnkwnk<<endl;
	wnk=0.46*pi;//wnk_max; hier ändern von wnkwnk zu
	x[0]=0.15;
  x[1]=0.;
  x[2]=cos(wnk)*w-x[1];		//Startgeschwindigkeit (x-Richtung) fuer den Trek
	x[3]=sin(wnk)*w+x[0];		//Startgeschwindigkeit (y-Richtung) fuer den Trek
  x[2]=x[2]+x[1]; // hier auskommentieren
  x[3]=x[3]-x[0]; // hier auskommentieren
	if (schalter){ //Symplectic Integration
		fout<<x[0]<<" "<<x[1]<<" "<< energy(x)<<endl;
		while(abb==false){
			l=x[0];
      x[2]=x[2]-x[1];// hier auskommentieren
      x[3]=x[3]+x[0];// hier auskommentieren
			Sympl_Integr(x,neqn,t,tstep,Hamilton);
			x[2]=x[2]+x[1];// hier auskommentieren
      x[3]=x[3]-x[0];// hier auskommentieren
      fout<<x[0]<<" "<<x[1]<<" "<< energy(x)<<endl;
			cond(x);
			cnt=count(x,l,cnt);
		}
	}
	else{ 			//Runge-Kutta Verfahren
		while(abb==false){
			l=x[0];
			rk(x,neqn,t,tstep,planet);
		  fout<<x[0]<<" "<<x[1]<<" "<< energy(x)<<endl;
   		cond(x);
			cnt=count(x,l,cnt);
	   }
	}
	cout<<cnt<<endl;
	fout.close();
}
