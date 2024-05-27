%This program calculates the relative error against variable L
clc;
clear all;
format long;
z1=2.627414291213963e+1;
z2= 3.239828844061399e+0;
z3= 1.446462700182916e+0;
z4= 1.039566130751374e+0;
L=1:1:10;       
N=15;
mu=1; 
kappa=1;
alpha=1;
Es_No_dB=0:10:30;
for t=1:4
ys(t)=10^(Es_No_dB(t)/10);%gamma bar
end
for t=1:10
%Calculation of Proposed Simpified Approximated SEP of BPSK    
I1(t)=0;
I2(t)=0;
I3(t)=0;
I4(t)=0;
for j=1:L(t)
for n=1:N
    I1(t)=I1(t)+alpha*mu^(mu+2*(j-1))*kappa^(j-1)*(1+kappa)^(mu+j-1)*1/gamma(mu+j-1)*1/factorial(j-1)*exp(-kappa*mu)*(-1)^(n-1)*mu^(n-1)*(1+kappa)^(n-1)*1/factorial(n-1)*(1/ys(1))^((n-1)*alpha/2+alpha/2*(mu+j-1))*(1/z1)^((n-1)*alpha/2+alpha/2*(mu+j-1))*gamma((n-1)*alpha/2+alpha/2*(mu+j-1));
    I2(t)=I2(t)+alpha*mu^(mu+2*(j-1))*kappa^(j-1)*(1+kappa)^(mu+j-1)*1/gamma(mu+j-1)*1/factorial(j-1)*exp(-kappa*mu)*(-1)^(n-1)*mu^(n-1)*(1+kappa)^(n-1)*1/factorial(n-1)*(1/ys(1))^((n-1)*alpha/2+alpha/2*(mu+j-1))*(1/z2)^((n-1)*alpha/2+alpha/2*(mu+j-1))*gamma((n-1)*alpha/2+alpha/2*(mu+j-1));
    I3(t)=I3(t)+alpha*mu^(mu+2*(j-1))*kappa^(j-1)*(1+kappa)^(mu+j-1)*1/gamma(mu+j-1)*1/factorial(j-1)*exp(-kappa*mu)*(-1)^(n-1)*mu^(n-1)*(1+kappa)^(n-1)*1/factorial(n-1)*(1/ys(1))^((n-1)*alpha/2+alpha/2*(mu+j-1))*(1/z3)^((n-1)*alpha/2+alpha/2*(mu+j-1))*gamma((n-1)*alpha/2+alpha/2*(mu+j-1));
    I4(t)=I4(t)+alpha*mu^(mu+2*(j-1))*kappa^(j-1)*(1+kappa)^(mu+j-1)*1/gamma(mu+j-1)*1/factorial(j-1)*exp(-kappa*mu)*(-1)^(n-1)*mu^(n-1)*(1+kappa)^(n-1)*1/factorial(n-1)*(1/ys(1))^((n-1)*alpha/2+alpha/2*(mu+j-1))*(1/z4)^((n-1)*alpha/2+alpha/2*(mu+j-1))*gamma((n-1)*alpha/2+alpha/2*(mu+j-1));
end
end
I(t)=1/16*(I1(t)+I2(t)+I3(t)+I4(t));    

%Calculation of Exact SEP of BPSK
SEP_exact(t)=0;
for j=1:L(t)
const=alpha.*mu.^(mu+2.*(j-1)).*kappa.^(j-1).*(kappa+1).^(mu+j-1).*1./2.*1./(gamma(mu+j-1)).*1./factorial(j-1).*exp(-kappa.*mu).*(1./ys(1)).^(alpha./2.*(mu+j-1));
fun = @(x)0.5.*const.*erfc(sqrt(x)).*x.^(alpha./2.*(mu+j-1)-1).*exp(-mu.*(1+kappa).*x.^(alpha./2).*(1./ys(1)).^(alpha./2));
SEP_exact(t)=SEP_exact(t)+integral(fun,0,Inf);
end
RE(t)=abs(SEP_exact(t)-I(t))/SEP_exact(t);
end
y=semilogy(L,RE,'-r' );
hold on;
