#include<iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
using namespace std;

class gen{
public:
int isDG,qlim,type;
double mp,nq,V0,w0,Qmax;
gen();
};

gen::gen ()
{
isDG=0;
qlim=0;
V0=1.0;
w0=1.0;
Qmax=10000000;
type=1;
}

class load{
public:
int type;
double Pl0,Ql0,alpha,beta,kpf,kqf;
load();
};

//  Choose alpha=0 and beta = 0 for constant PQ type of load 
load::load()
{
type=0;
Pl0=0;
Ql0=0;
//alpha=0;
alpha=2.0;
kpf=0;
kqf=0;
//beta=0;
beta=2.0;
}
complex<double> sum(complex<double> * x,int n)
{
int i;
complex<double> y;
for(i=1;i<=n;i++)
y=y + x[i];
return y;
}

double max(double* x,int n)
{int i;
double y;
y=abs(x[1]);
for(i=2;i<=n;i++)
if(abs(x[i])>y)
y=abs(x[i]);

return y;
}

int min(complex<double>* x,int n)
{int i;
double y;
y=abs(x[1]);
int index=1;
for(i=2;i<=n;i++)
if(abs(x[i])<y)
{
y=abs(x[i]);
index=i;
}
return index;
}

int invshipley(double a[250][250],int n )
{
int i,j,k,precision_no=6;
for(k=1;k<=n;k++)
{
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
{
if(i==k || j==k)
continue;
a[i][j]=a[i][j] - (a[i][k]*a[k][j])/a[k][k];
}
for(i=1;i<=n;i++)
{if (i==k)
continue;
a[i][k]=(-1*a[i][k])/a[k][k];
}

for(j=1;j<=n;j++)
{if (j==k)
continue;
a[k][j]=(-1.0*a[k][j])/a[k][k];
}

a[k][k]=(-1.0/a[k][k]);
}
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
a[i][j]=-1.0*a[i][j];

return 0;
}

int main()
{

int i,j,l,from_bus,to_bus,k,no_iterations,no_buses,no_lines,temp,precision_no=4,no_PQ,no_pv,no_DG,linedata[120][3],no_tie;
double tolerance,base_V,base_VA,r,x,J[250][250],delP[250],delV[250],Vang,Vmag,Ploss,p,Qloss,Psys,Qsys;
complex<double> z,Z_base,SL[120],SG[120],Y[120][120],Scal[120],V[120],I[120],SLtotal,line_z[120],Sload;
string comment;
fstream fp;
double w0,PI=3.141592654,w=1.0,V0=1.0;
w0=2*PI*60;

gen DG[120];
load LD[120];

fp.open("islandip.txt",ios::in);

if(fp.is_open())
{
getline(fp,comment);
fp>>no_buses;
fp>>no_lines;
fp>>tolerance;
fp>>no_iterations;
fp>>base_V;
fp>>base_VA;
fp>>no_DG;
fp>>no_tie;
getline(fp,comment);
getline(fp,comment);
Z_base=(base_V*base_V)/base_VA;
for(i=1;i<=no_lines;i++)
{
fp>>temp;
fp>>from_bus;
fp>>to_bus;
fp>>r;
fp>>x;

z=complex<double>(r,(w0/1000.0)*x);

line_z[i]=z;

Y[from_bus][to_bus]=Y[from_bus][to_bus] - Z_base/z;
Y[from_bus][from_bus]=Y[from_bus][from_bus] + Z_base/z;
Y[to_bus][to_bus]=Y[to_bus][to_bus] + Z_base/z;
Y[to_bus][from_bus]=Y[from_bus][to_bus];

linedata[i][1]=from_bus;
linedata[i][2]=to_bus;
}
getline(fp,comment);
getline(fp,comment);

for(i=1;i<=no_buses;i++)
{
fp>>j;
fp>>r;
fp>>x;
fp>>temp;
if(temp==1)
{
LD[j].alpha=1.51;
LD[j].beta=3.4;
}
if(temp==2)
{
LD[j].alpha=0.18;
LD[j].beta=6.0;
}
if(temp==3)
{
LD[j].alpha=0.92;
LD[j].beta=4.04;
}

LD[j].Pl0=r*1000/base_VA;
LD[j].Ql0=x*1000/base_VA;
}
for(i=1;i<=no_buses;i++)
Sload+=SL[i];
getline(fp,comment);
getline(fp,comment);
for(i=1;i<=no_DG;i++)
{
fp>>temp;
fp>>r;
DG[temp].mp=r;
fp>>x;
DG[temp].nq=x;
fp>>x;
DG[temp].Qmax=x;
DG[temp].isDG=1;
fp>>r;
DG[temp].V0=r;
fp>>DG[temp].type;
}

getline(fp,comment);
getline(fp,comment);

for(i=1;i<=no_buses;i++)
getline(fp,comment);
getline(fp,comment);

for(i=2;i<=(no_pv +1);i++)
{
fp>>p;
fp>>Vmag;
V[i]=polar(Vmag,0.0);
SL[i]=(complex<double>(-p,0))/(base_VA);
}
}

else
{cout<<"Unable to open file\n";
return 0;
}


fp.close();

fp.open("nrlfop.txt",ios::out);

V[1]=1.0 + 0i;
for(i=2;i<=no_buses;i++)
V[i]= 1. + 0i;
for(k=1;k<=no_iterations;k++)
{
//Ybus
for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
Y[i][j]=complex<double>(0,0);

for(i=1;i<=no_lines;i++)
{
from_bus=linedata[i][1];
to_bus=linedata[i][2];
r=real(line_z[i]);
x=imag(line_z[i]);
z=complex<double>(r,w*x);
Y[from_bus][to_bus]=Y[from_bus][to_bus] - Z_base/z;

Y[from_bus][from_bus]=Y[from_bus][from_bus] + Z_base/z;

Y[to_bus][to_bus]=Y[to_bus][to_bus] + Z_base/z;

Y[to_bus][from_bus]=Y[from_bus][to_bus];
}

for(i=1;i<=no_buses;i++)
if(DG[i].isDG)
{
if(DG[i].type==1)
{
r=(1/DG[i].mp)*(DG[i].w0-w);
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==2)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==3)
{
double a,b;
a=(1/DG[i].mp)*(DG[i].w0-w);
b=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
r=(a+b)/2.0;
x=(b-a)/2.0;
}
DG[i].qlim=0;
if(x>DG[i].Qmax)
{
x=DG[i].Qmax;
DG[i].qlim=1;
}

if(x<-DG[i].Qmax)
{
x=-DG[i].Qmax;
DG[i].qlim=1;
}

SG[i]=complex<double>(r,x);
}

Sload=0+0i;
for(i=1;i<=no_buses;i++)
{
r=LD[i].Pl0*(pow(abs(V[i]),LD[i].alpha))*(1.0+LD[i].kpf*(w-1.0));
x=LD[i].Ql0*(pow(abs(V[i]),LD[i].beta))*(1.0+LD[i].kqf*(w-1.0));
SL[i]=complex<double>(r,x);
Sload+=SL[i];
}

for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
J[i][j]=0;

for(i=1;i<=no_buses;i++)
{complex<double> temp=0;
for(j=1;j<=no_buses;j++)
temp+=Y[i][j]*V[j];
I[i]=temp;
}

for(i=1;i<=no_buses;i++)
Scal[i]=(V[i]*conj(I[i]));

//    J11
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
{
if(i==j)
J[i-1][j-1]=-(imag(Scal[i]) + abs(V[i])*abs(V[i])*imag(Y[i][i]));
else
J[i-1][j-1]=-(abs(V[i])*abs(V[j])*abs(Y[i][j])*sin( arg(Y[i][j]) + arg(V[j]) - arg(V[i])));
}
//J21
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
J[i-2+no_buses][j-1]= real(Scal[i]) - abs(V[i])*abs(V[i])*real(Y[i][i]);
else
J[i-2+no_buses][j-1]= -abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J12
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
{
J[i-1][j-2+no_buses]= real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]) + real(SL[i])*LD[i].alpha;
if(DG[i].isDG)
{
if(DG[i].type==2)
J[i-1][j-2+no_buses]+=abs(V[i])/DG[i].nq;
if(DG[i].type==3)
J[i-1][j-2+no_buses]+=0.5*abs(V[i])/DG[i].nq;
}
}
else
J[i-1][j-2+no_buses]=abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J22
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
{
J[i-2+no_buses][j-2+no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]) + LD[i].beta*imag(SL[i]);
if(DG[i].isDG&&(DG[i].qlim==0))
{
if(DG[i].type==1)
J[i-2+no_buses][j-2+no_buses]+=(1/DG[i].nq)*abs(V[i]);
if(DG[i].type==3)
J[i-2+no_buses][j-2+no_buses]+=(0.5/DG[i].nq)*abs(V[i]);
}
}
else
J[i-2+no_buses][j-2+no_buses]=-abs(V[i])*abs(V[j])*abs(Y[i][j])*sin(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J13
j=2*no_buses - 1;
for(i=2;i<=no_buses;i++)
{
complex<double> temp=0+0i;
for(l=1;l<=no_buses;l++)
{
if(abs(Y[i][l]))
{
r=real(complex<double>(-1.0,0)/Y[i][l]);
temp+=(V[i]-V[l])*(Y[i][l]/w)*(complex<double>(1,0) + Y[i][l]*complex<double>(r,0))*conj(V[i]);
}
}
complex<double> Yii=0+0i;
int z;
for(z=1;z<=no_buses;z++)
Yii+=Y[i][z];
if(abs(Yii))
{
r=real(complex<double>(1.0,0)/Yii);
temp-=V[i]*(Yii/w)*(complex<double>(1,0)+complex<double>(r,0)*Yii)*conj(V[i]);
}
J[i-1][j]=real(temp) + LD[i].Pl0*pow(abs(V[i]),LD[i].alpha)*LD[i].kpf;
J[i-2+no_buses][j]=-imag(temp) + LD[i].Ql0*pow(abs(V[i]),LD[i].beta)*LD[i].kqf;
if(DG[i].isDG)
{
if(DG[i].type==1)
J[i-1][j]+=(1/DG[i].mp);
if(DG[i].type==2)
J[i-2+no_buses][j]-=(1/DG[i].mp);
if(DG[i].type==3)
{
J[i-1][j]+=(0.5/DG[i].mp);
J[i-2+no_buses][j]-=(0.5/DG[i].mp);
}
}
}

//J14
j=2*no_buses;
for(i=2;i<=no_buses;i++)
{
J[i-1][j]=abs(V[i])*abs(V[1])*abs(Y[i][1])*cos(arg(Y[i][1]) + arg(V[1]) - arg(V[i]));//abs(V[i])*abs(Y[i][1])*cos( arg(V[i])  - arg(Y[i][1]) );
}

//J24
j=2*no_buses;
for(i=2;i<=no_buses;i++)
{
J[i-2+no_buses][j]=-abs(V[i])*abs(V[1])*abs(Y[i][1])*sin(arg(Y[i][1]) + arg(V[1]) - arg(V[i]));//abs(V[i])*abs(Y[i][1])*sin( arg(V[i])  - arg(Y[i][1]) );
}

//J31
for(j=2;j<=no_buses;j++)
J[2*no_buses - 1][j-1]=-(abs(V[1])*abs(V[j])*abs(Y[1][j])*sin( arg(Y[1][j]) + arg(V[j]) - arg(V[1])));

//J32
for(j=2;j<=no_buses;j++)
{
J[2*no_buses - 1][j-2+no_buses]=abs(V[1])*abs(V[j])*abs(Y[1][j])*cos(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J41
i=2*no_buses;
for(j=2;j<=no_buses;j++)
{
J[i][j-1]=-abs(V[1])*abs(V[j])*abs(Y[1][j])*cos(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J42
i=2*no_buses;
for(j=2;j<=no_buses;j++)
{
J[i][j-2+no_buses]=-abs(V[1])*abs(V[j])*abs(Y[1][j])*sin(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J33 & J43
i=1;
{
complex<double> temp=0+0i;
for(l=1;l<=no_buses;l++)
{
if(abs(Y[i][l]))
{
r=real(complex<double>(-1.0,0)/Y[i][l]);
temp+=(V[i]-V[l])*(Y[i][l]/w)*(complex<double>(1,0) + Y[i][l]*complex<double>(r,0))*conj(V[i]);
}
}
complex<double> Yii=0+0i;
int z;
for(z=1;z<=no_buses;z++)
Yii+=Y[i][z];
if(abs(Yii))
{
r=real(complex<double>(1.0,0)/Yii);
temp-=V[i]*(Yii/w)*(complex<double>(1,0)+complex<double>(r,0)*Yii)*conj(V[i]);
}
J[2*no_buses - 1][2*no_buses - 1]=real(temp) + LD[i].Pl0*pow(abs(V[i]),LD[i].alpha)*LD[i].kpf*w;
J[2*no_buses][2*no_buses - 1]=-imag(temp) + LD[i].Ql0*pow(abs(V[i]),LD[i].beta)*LD[i].kqf*w;
if(DG[i].isDG)
{
if(DG[i].type==1)
J[2*no_buses - 1][2*no_buses - 1]+=(1/DG[i].mp);
if(DG[i].type==2)
J[2*no_buses][2*no_buses - 1]-=(1/DG[i].mp);
if(DG[i].type==3)
{
J[2*no_buses - 1][2*no_buses - 1]+=(0.5/DG[i].mp);
J[2*no_buses][2*no_buses - 1]-=(0.5/DG[i].mp);
}
}
}

//J34
i=1;
J[2*no_buses - 1][2*no_buses]=real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]) + real(SL[i])*LD[i].alpha;
if(DG[i].isDG)
{
if(DG[i].type==2)
J[2*no_buses - 1][2*no_buses]+=abs(V[i])/DG[i].nq;
if(DG[i].type==3)
J[2*no_buses - 1][2*no_buses]+=0.5*abs(V[i])/DG[i].nq;
}

//J44
i=1;
j=1;
J[2*no_buses][2*no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]);
if(DG[i].isDG&&(DG[i].qlim==0))
{
if(DG[i].type==1)
J[2*no_buses][2*no_buses]+=(1/DG[i].nq)*abs(V[i]);
if(DG[i].type==3)
J[2*no_buses][2*no_buses]+=(0.5/DG[i].nq)*abs(V[i]);
}

invshipley(J,2*no_buses);

for(i=2;i<=no_buses;i++)
delP[i-1]=-(real(Scal[i]) + real(SL[i]) -real(SG[i]));
for(i=no_buses;i<=2*(no_buses - 1);i++)
delP[i]=-(imag(Scal[i-no_buses+2]) + imag(SL[i-no_buses+2]) - imag(SG[i-no_buses+2]));


Ploss=0;
for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
Ploss+=real( Y[i][j]*(V[i]*conj(V[j]) + V[j]*conj(V[i])) );
Ploss=Ploss*0.5;

Qloss=0;
for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
Qloss-=imag( Y[i][j]*(V[i]*conj(V[j]) + V[j]*conj(V[i])) );
Qloss=Qloss*0.5;
i=1;
delP[2*no_buses - 1]=-(real(Scal[i]) + real(SL[i]) -real(SG[i]));
delP[2*no_buses]=-(imag(Scal[i]) + imag(SL[i]) - imag(SG[i]));

for(i=1;i<=2*(no_buses);i++)
{
double temp=0;
for(j=1;j<=2*(no_buses);j++)
temp+=J[i][j]*delP[j];
delV[i]=temp;
}

for(i=1;i<=(no_buses - 1);i++)
{
Vang=arg(V[i+1])+delV[i];
Vmag=abs(V[i+1])*(1.0 + delV[i+ no_buses - 1]);
V[i+1]=polar(Vmag,Vang);
}

w+=delV[2*no_buses - 1];
Vmag=abs(V[1])*(1.0 + delV[2*no_buses]);
V[1]=polar(Vmag,0.0);

if(tolerance>max(delP,(2*no_buses)))  
break;
}

fp<<std::fixed<<std::setprecision(precision_no);
if(k<=no_iterations)
fp<<"\n\nSolution converged after "<<k<<" iterations\n";
else
fp<<"\n\nSolution did not converge\n";
Psys=0;
Qsys=0;
fp<<"DG_bus\tPdg\tQdg\n";
for(i=1;i<=no_buses;i++)
if(DG[i].isDG)
{
if(DG[i].type==1)
{
r=(1/DG[i].mp)*(DG[i].w0-w);
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==2)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==3)
{
double a,b;
a=(1/DG[i].mp)*(DG[i].w0-w);
b=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
r=(a+b)/2.0;
x=(b-a)/2.0;
}
DG[i].qlim=0;

if(x>DG[i].Qmax)
{
x=DG[i].Qmax;
DG[i].qlim=1;
}

if(x<-DG[i].Qmax)
{
x=-DG[i].Qmax;
DG[i].qlim=1;
}
Qsys+=x;
Psys+=r;
fp<<i<<"\t"<<r<<"\t"<<x<<"\n";
}
fp<<"\nBus_no\t\tVoltage(magnitude)\t\tAngle(deg)\n";
fp<<"____________________________________________________________\n";
for(i=1;i<=no_buses;i++)
fp<<i<<"\t\t"<<abs(V[i])<<"\t\t\t\t"<<arg(V[i])*(180/PI)<<"\n";
fp<<"____________________________________________________________\n";

fp<<"\nScal\tSL\n";
for(i=1;i<=no_buses;i++)
fp<<Scal[i]<<"\t"<<SL[i]<<"\n";
SLtotal=sum(Scal,no_buses);
Ploss=real(SLtotal);
fp<<"\nPG="<<Psys*base_VA<<" W\tQG="<<Qsys*base_VA<<" VAR\nTotal P loss ="<<real(SLtotal)*base_VA<<" W\t\tTotal Q loss ="<<imag(SLtotal)*base_VA<<" VAR\n";
fp<<std::fixed<<std::setprecision(6);
fp<<"w="<<w<<endl;
fp<<"\n\nVoltage is minimum at bus "<<min(V,no_buses);
fp.close();
return 0;
}
