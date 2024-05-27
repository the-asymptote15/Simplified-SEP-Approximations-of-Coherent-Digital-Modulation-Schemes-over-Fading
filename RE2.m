%This program calculates the relative error in the upper bound of the truncation error
clc;
clear all;
format long;
z1=2.627414291213963e+1;
z2= 3.239828844061399e+0;
z3= 1.446462700182916e+0;
z4= 1.039566130751374e+0;
L=1:1:20;       
N=15;
mu=1; 
kappa=1;
alpha=1;
Es_No_dB=0:10:30;
for t=1:4
ys(t)=10^(Es_No_dB(t)/10);%gamma bar
end
for t=1:20
const_upperB(t)=1/8*(-1)^N*alpha*(mu*(1+kappa))^(mu+N)*exp(-kappa*mu)*1/2*1/factorial(N)*1/(ys(1))^(0.5*alpha*(mu+N))*mu^(2*L(t))*kappa^(L(t))*(1+kappa)^(L(t))*gamma(0.5*alpha*(mu+L(t)+N))*1/gamma(mu+L(t))*1/factorial(L(t));
c1=0.5*alpha*(mu+L(t)+N);
c2=-mu*(1+kappa)/(ys(1)*z1)^(alpha/2);
c3=-mu*(1+kappa)/(ys(1)*z2)^(alpha/2);
c4=-mu*(1+kappa)/(ys(1)*z3)^(alpha/2);
c5=-mu*(1+kappa)/(ys(1)*z4)^(alpha/2);
I(t)=const_upperB(t)*(z1^(-0.5*alpha*(mu+L(t)+N))*hypergeom([1, c1], N+1, c2)+z2^(-0.5*alpha*(mu+L(t)+N))*hypergeom([1, c1], N+1, c3)+z3^(-0.5*alpha*(mu+L(t)+N))*hypergeom([1, c1], N+1, c4)+z4^(-0.5*alpha*(mu+L(t)+N))*hypergeom([1, c1], N+1, c5));
%Calculation of Exact SEP of BPSK
SEP_exact(t)=0;
for j=1:L(t)
const=alpha.*mu.^(mu+2.*(j-1)).*kappa.^(j-1).*(kappa+1).^(mu+j-1).*1./2.*1./(gamma(mu+j-1)).*1./factorial(j-1).*exp(-kappa.*mu).*(1./ys(1)).^(alpha./2.*(mu+j-1));
fun = @(x)0.5.*const.*erfc(sqrt(x)).*x.^(alpha./2.*(mu+j-1)-1).*exp(-mu.*(1+kappa).*x.^(alpha./2).*(1./ys(1)).^(alpha./2));
SEP_exact(t)=SEP_exact(t)+integral(fun,0,Inf);
end
%RE in the upper bound of the truncation error
RE_upperB(t)=abs(I(t))/SEP_exact(t);
end
y=semilogy(L,RE_upperB,'-g' );
hold on;
