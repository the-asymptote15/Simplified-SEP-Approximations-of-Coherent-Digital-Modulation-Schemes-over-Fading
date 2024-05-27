%Fig.:8 of the manuscript is obtained using this code
tic
clear all;
alpha=1;
alpha2=alpha./2;
kappa=0.5;
mu=1;
M=16;
N=153;
N1=100;
sigma=sqrt(3./(M-1));
s=30;
z_1=0:5:s;
t1=10.^(z_1./10);

A=4.*(1-(1./sqrt(M)));
z=0:1:s;
t=10.^(z./10);
z1=26.27414291213963;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    n1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*sigma.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n1(i)=n1(i)+s1(j);
    end
end

z2=3.239828844061399;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    n2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*sigma.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n2(i)=n2(i)+s1(j);
    end
end

z3=1.446462700182916;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    n3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*sigma.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n3(i)=n3(i)+s1(j);
    end
end

z4=1.039566130751374;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    n4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*sigma.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n4(i)=n4(i)+s1(j);
    end
end
summ=0;
for i=1:s+1
    summ(i)=(n4(i)+n1(i)+n2(i)+n3(i));
end


%2nd term

A1=4.*(1-(1./sqrt(M))).^2;
%z1=26.27414291213963;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*sigma.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m1(i)=m1(i)+s1(j);
    end
end

%z2=3.239828844061399;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*sigma.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m2(i)=m2(i)+s1(j);
    end
end

%z3=1.446462700182916;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*sigma.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m3(i)=m3(i)+s1(j);
    end
end

%z4=1.039566130751374;
for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*sigma.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m4(i)=m4(i)+s1(j);
    end
end

for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m5(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z2).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m5(i)=m5(i)+s1(j);
    end
end


for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m6(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z3+z4).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m6(i)=m6(i)+s1(j);
    end
end

for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m7(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z3).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m7(i)=m7(i)+s1(j);
    end
end

for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m8(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z3).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m8(i)=m8(i)+s1(j);
    end
end

for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m9(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z4).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m9(i)=m9(i)+s1(j);
    end
end

for i=1:s+1
    y=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y=[y,coeff1];
    end
    m10(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z4).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m10(i)=m10(i)+s1(j);
    end
end
summ1=0;
for i=1:s+1
    summ1(i)=(m1(i)+m2(i)+m3(i)+m4(i)+m5(i)+m6(i)+m7(i)+m8(i)+m9(i)+m10(i));
end
summ2=0;
for i=1:s+1
    summ2(i)=summ(i)-summ1(i);
end
syms j x
% for i=1:31
%     f=@(x,j) (alpha2.*(1./2).*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sigma.*sqrt(x./2)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
%     f2= symsum(f,j,0,10);
%     ht = matlabFunction(f2);
%     I(i)=integral(ht,0,10);
% end
for i=1:s+1
    %f=@(x,j) ((1./4).*alpha2.*A1.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)))).*((A.*0.53.*erfc(sigma.*sqrt(x./2)))-(A1.*(1./4).*(erfc(sigma.*sqrt(x./2))).^2));
    f2= symsum(f,j,0,100);
    ht = matlabFunction(f2);
    I(i)=integral(ht,0,100);
    
end
% I2=0;
% for i=1:31
%     I2(i)=I(i)-I1(i);
% end
syms j x
% for i=1:31
%     f=@(x,j) (alpha2.*(1./2).*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sigma.*sqrt(x./2)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
%     f2= symsum(f,j,0,10);
%     ht = matlabFunction(f2);
%     I(i)=integral(ht,0,10);
% end
for i=1:7
    %f=@(x,j) ((1./4).*alpha2.*A1.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t1(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t1(i).^((alpha2).*(mu+j)))).*((A.*0.53.*erfc(sigma.*sqrt(x./2)))-(A1.*(1./4).*(erfc(sigma.*sqrt(x./2))).^2));
    f2= symsum(f,j,0,100);
    ht = matlabFunction(f2);
    I1(i)=integral(ht,0,100);
    
end

%legend('Exact','Proposed','Simulations')



%HQAM

% tic
% clear all;
% alpha=1;
% alpha2=alpha./2;
% kappa=1.5;
% mu=1;
% M=16;
% N=145;
% N1=10;
% s=30;
KCNN=6*(1-M^(-0.5))^2;
KNN=2*(3-(4*(M^(-0.5)))+M^(-1));
Z=24/(7*M-4);
% z=0:1:s;
% t=10.^(z./10);
% z_1=0:5:s;
% t1=10.^(z_1./10);
% z1=26.27414291213963;
% z2=3.239828844061399;
% z3=1.446462700182916;
% z4=1.039566130751374;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    o1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f=(y2(l).*B(j).*KNN.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*Z)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f;
      end
      o1(i)=o1(i)+s1(j);
    end
end


for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    o_2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*KNN*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*Z)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      o_2(i)=o_2(i)+s1(j);
    end
end


for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    o_3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*KNN.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*Z)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      o_3(i)=o_3(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    o_4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*KNN.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*Z)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      o_4(i)=o_4(i)+s1(j);
    end
end
summ3=0;
for i=1:s+1
    summ3(i)=(o_4(i)+o1(i)+o_2(i)+o_3(i));
end


%2nd term

A11=(2/3)*KCNN;
%z1=26.27414291213963;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1.*2.*Z)./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_1(i)=r_1(i)+s1(j);
    end
end

%z2=3.239828844061399;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2.*2.*Z)./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_2(i)=r_2(i)+s1(j);
    end
end

%z3=1.446462700182916;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((((z3.*2.*Z)./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_3(i)=r_3(i)+s1(j);
    end
end

%z4=1.039566130751374;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((((z4.*2.*Z)./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_4(i)=r_4(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_5(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z2).*(Z./3))).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_5(i)=r_5(i)+s1(j);
    end
end


for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_6(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z3+z4).*(Z./3))).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_6(i)=r_6(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_7(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z3).*(Z./3))).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_7(i)=r_7(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_8(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z3).*(Z./3))).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_8(i)=r_8(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_9(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z4).*(Z./3))).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_9(i)=r_9(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    r_10(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z4).*(Z./3))).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      r_10(i)=r_10(i)+s1(j);
    end
end
summe=0;
for i=1:s+1
    summe(i)=(r_1(i)+r_2(i)+r_3(i)+r_4(i)+r_5(i)+r_6(i)+r_7(i)+r_8(i)+r_9(i)+r_10(i));
end


%3 term

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((z1.*Z.*(2./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_1(i)=q_1(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z1+z2))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_2(i)=q_2(i)+s1(j);
    end
end
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z1+z3))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_3(i)=q_3(i)+s1(j);
    end
end

z4=1.039566130751374;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z1+z4))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_4(i)=q_4(i)+s1(j);
    end
end
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_5(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((z2.*Z.*(2./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_5(i)=q_5(i)+s1(j);
    end
end

%z2=3.239828844061399;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_6(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z2+z1))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_6(i)=q_6(i)+s1(j);
    end
end

%z3=1.446462700182916;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_7(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z2+z3))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_7(i)=q_7(i)+s1(j);
    end
end

%z4=1.039566130751374;
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_8(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z2+z4))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_8(i)=q_8(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_9(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z3+z1))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_9(i)=q_9(i)+s1(j);
    end
end


for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_10(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z3+z2))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_10(i)=q_10(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_11(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((z3.*Z.*(2./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_11(i)=q_11(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_12(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z3+z4))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_12(i)=q_12(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_13(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z4+z1))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_13(i)=q_13(i)+s1(j);
    end
end

for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_14(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z4+z2))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_14(i)=q_14(i)+s1(j);
    end
end
for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_15(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((Z.*(3.*z4+z3))./6).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_15(i)=q_15(i)+s1(j);
    end
end


for i=1:s+1
    y2=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y2=[y2,coeff1];
    end
    q_16(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y2(l).*B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./((z4.*Z.*(2./3)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q_16(i)=q_16(i)+s1(j);
    end
end

summ5=0;
for i=1:s+1
    summ5(i)=(q_1(i)+q_2(i)+q_3(i)+q_4(i)+q_5(i)+q_6(i)+q_7(i)+q_8(i)+q_9(i)+q_10(i)+q_11(i)+q_12(i)+q_13(i)+q_14(i)+q_15(i)+q_16(i));
end
sumf=0;
for i=1:s+1
    sumf(i)=summ3(i)+summe(i)-summ5(i);
end
syms j x
% for i=1:31
%     f=@(x,j) (alpha2.*(1./2).*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sigma.*sqrt(x./2)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
%     f2= symsum(f,j,0,10);
%     ht = matlabFunction(f2);
%     I(i)=integral(ht,0,10);
% end
for i=1:s+1
    %f=@(x,j) ((1./4).*alpha2.*A1.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)))).*((KNN.*(1./2).*erfc(sqrt((x.*Z)./2)))+((1./4).*A11.*(erfc(sqrt(Z.*x./3))).^2)-((1./4).*2.*KCNN.*erfc(sqrt(Z.*x./2)).*erfc(sqrt(Z.*x./6))));
    f2= symsum(f,j,0,10);
    ht = matlabFunction(f2);
    I3(i)=integral(ht,0,100);
    
end

syms j x

for i=1:7
    %f=@(x,j) ((1./4).*alpha2.*A1.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t1(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t1(i).^((alpha2).*(mu+j)))).*((KNN.*(1./2).*erfc(sqrt((x.*Z)./2)))+((1./4).*A11.*(erfc(sqrt(Z.*x./3))).^2)-((1./4).*2.*KCNN.*erfc(sqrt(Z.*x./2)).*erfc(sqrt(Z.*x./6))));
    f2= symsum(f,j,0,100);
    ht = matlabFunction(f2);
    I4(i)=integral(ht,0,100);
    
end
% I2=0;
% for i=1:31
%     I2(i)=I(i)-I1(i);
% end
%SQAQM



for i=1:s+1
    na(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1))))./(((z1.*sigma.^2)./2).^(alpha2.*(mu+(j-1))));
      na(i)=na(i)+f;
    end
end

z2=3.239828844061399;
for i=1:s+1
    nb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1))))./(((z2.*sigma.^2)./2).^(alpha2.*(mu+(j-1))));
      nb(i)=nb(i)+f;
    end
end

z3=1.446462700182916;
for i=1:s+1
    nc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1))))./(((z3.*sigma.^2)./2).^(alpha2.*(mu+(j-1))));
      nc(i)=nc(i)+f;
    end
end

z4=1.039566130751374;
for i=1:s+1
    nd(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(0.1325).*gamma(alpha2.*(mu+(j-1))))./(((z4.*sigma.^2)./2).^(alpha2.*(mu+(j-1))));
      nd(i)=nd(i)+f;
    end
end
summa=0;
for i=1:s+1
    summa(i)=(nd(i)+nc(i)+nb(i)+na(i));
end


%2nd term

A1=4.*(1-(1./sqrt(M))).^2;
%z1=26.27414291213963;
for i=1:s+1
    ma(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z1.*sigma.^2)).^(alpha2.*(mu+(j-1))));
      ma(i)=ma(i)+f;
    end
end

%z2=3.239828844061399;
for i=1:s+1
    mb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z2.*sigma.^2)).^(alpha2.*(mu+(j-1))));
      mb(i)=mb(i)+f;
    end
end

%z3=1.446462700182916;
for i=1:s+1
    mc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z3.*sigma.^2)).^(alpha2.*(mu+(j-1))));
      mc(i)=mc(i)+f;
    end
end

%z4=1.039566130751374;
for i=1:s+1
    md(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z4.*sigma.^2)).^(alpha2.*(mu+(j-1))));
      md(i)=md(i)+f;
    end
end

for i=1:s+1
    me(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z2).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1))));
      me(i)=me(i)+f;
    end
end


for i=1:s+1
    mf(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z3+z4).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1))));
      mf(i)=mf(i)+f;
    end
end

for i=1:s+1
    mg(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z3).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1))));
  
      mg(i)=mg(i)+f;
    end
end

for i=1:s+1
    mh(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z3).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1))));
      mh(i)=mh(i)+f;
    end
end

for i=1:s+1
    mi(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z4).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1))));
      mi(i)=mi(i)+f;
    end
end

for i=1:s+1
    mj(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z4).*(sigma.^2)./2)).^(alpha2.*(mu+(j-1))));
     
      mj(i)=mj(i)+f;
    end
end
summb=0;
for i=1:s+1
    summb(i)=(ma(i)+mb(i)+mc(i)+md(i)+me(i)+mf(i)+mg(i)+mh(i)+mi(i)+mj(i));
end
summc=0;
for i=1:s+1
    summc(i)=summa(i)-summb(i);
end





%HQAM


KCNN=6*(1-M^(-0.5))^2;
KNN=2*(3-(4*(M^(-0.5)))+M^(-1));
Z=24/(7*M-4);
% z=0:1:s;
% t=10.^(z./10);
% z_1=0:5:s;
% t1=10.^(z_1./10);
% z1=26.27414291213963;
% z2=3.239828844061399;
% z3=1.446462700182916;
% z4=1.039566130751374;
for i=1:s+1
    oa(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*KNN.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z1.*Z)./2).^(alpha2.*(mu+(j-1))));
 
      oa(i)=oa(i)+f;
    end
end


for i=1:s+1
    ob(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*KNN*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z2.*Z)./2).^(alpha2.*(mu+(j-1))));
      ob(i)=ob(i)+f;
    end
end


for i=1:s+1
    oc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*KNN.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z3.*Z)./2).^(alpha2.*(mu+(j-1))));
 
      oc(i)=oc(i)+f;
    end
end

for i=1:s+1
    od(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*KNN.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z4.*Z)./2).^(alpha2.*(mu+(j-1))));

      od(i)=od(i)+f;
    end
end
summd=0;
for i=1:s+1
    summd(i)=(od(i)+oa(i)+ob(i)+oc(i));
end


%2nd term

A11=(2/3)*KCNN;
%z1=26.27414291213963;
for i=1:s+1
    ra(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((((z1.*2.*Z)./3)).^(alpha2.*(mu+(j-1))));
     
      ra(i)=ra(i)+f;
    end
end

%z2=3.239828844061399;
for i=1:s+1
    rb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((((z2.*2.*Z)./3)).^(alpha2.*(mu+(j-1))));
 
      rb(i)=rb(i)+f;
    end
end

%z3=1.446462700182916;
for i=1:s+1
    rc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((((z3.*2.*Z)./3)).^(alpha2.*(mu+(j-1))));
      rc(i)=rc(i)+f;
    end
end

%z4=1.039566130751374;
for i=1:s+1
    rd(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((((z4.*2.*Z)./3)).^(alpha2.*(mu+(j-1))));
 
      rd(i)=rd(i)+f;
    end
end

for i=1:s+1
    re(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z2).*(Z./3))).^(alpha2.*(mu+(j-1))));

      re(i)=re(i)+f;
    end
end


for i=1:s+1
    rf(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z3+z4).*(Z./3))).^(alpha2.*(mu+(j-1))));
      rf(i)=rf(i)+f;
    end
end

for i=1:s+1
    rg(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z3).*(Z./3))).^(alpha2.*(mu+(j-1))));
      rg(i)=rg(i)+f;
    end
end

for i=1:s+1
    rh(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z3).*(Z./3))).^(alpha2.*(mu+(j-1))));
      rh(i)=rh(i)+f;
    end
end

for i=1:s+1
    ri(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z4).*(Z./3))).^(alpha2.*(mu+(j-1))));
      ri(i)=ri(i)+f;
    end
end

for i=1:s+1
    rj(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A11.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z4).*(Z./3))).^(alpha2.*(mu+(j-1))));
      rj(i)=rj(i)+f;
    end
end
summe=0;
for i=1:s+1
    summe(i)=(ra(i)+rb(i)+rc(i)+rd(i)+re(i)+rf(i)+rg(i)+rh(i)+ri(i)+rj(i));
end


%3 term

for i=1:s+1
    qa(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((z1.*Z.*(2./3)).^(alpha2.*(mu+(j-1))));
      qa(i)=qa(i)+f;
    end
end

for i=1:s+1
    qb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z1+z2))./6).^(alpha2.*(mu+(j-1))));
      qb(i)=qb(i)+f;
    end
end
for i=1:s+1
    qc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z1+z3))./6).^(alpha2.*(mu+(j-1))));
      qc(i)=qc(i)+f;
    end
end

z4=1.039566130751374;
for i=1:s+1
    qd(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z1+z4))./6).^(alpha2.*(mu+(j-1))));
      qd(i)=qd(i)+f;
    end
end
for i=1:s+1
    qe(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((z2.*Z.*(2./3)).^(alpha2.*(mu+(j-1))));
      qe(i)=qe(i)+f;
    end
end

%z2=3.239828844061399;
for i=1:s+1
    qf(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z2+z1))./6).^(alpha2.*(mu+(j-1))));
      qf(i)=qf(i)+f;
    end
end

%z3=1.446462700182916;
for i=1:s+1
    qg(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z2+z3))./6).^(alpha2.*(mu+(j-1))));
      qg(i)=qg(i)+f;
    end
end

%z4=1.039566130751374;
for i=1:s+1
    qh(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z2+z4))./6).^(alpha2.*(mu+(j-1))));
      qh(i)=qh(i)+f;
    end
end

for i=1:s+1
    qi(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z3+z1))./6).^(alpha2.*(mu+(j-1))));
      qi(i)=qi(i)+f;
    end
end


for i=1:s+1
    qj(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z3+z2))./6).^(alpha2.*(mu+(j-1))));
      qj(i)=qj(i)+f;
    end
end

for i=1:s+1
    qk(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((z3.*Z.*(2./3)).^(alpha2.*(mu+(j-1))));
      qk(i)=qk(i)+f;
    end
end

for i=1:s+1
    ql(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z3+z4))./6).^(alpha2.*(mu+(j-1))));
      ql(i)=ql(i)+f;
    end
end

for i=1:s+1
    qm(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z4+z1))./6).^(alpha2.*(mu+(j-1))));
      qm(i)=qm(i)+f;
    end
end

for i=1:s+1
    qn(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z4+z2))./6).^(alpha2.*(mu+(j-1))));
      qn(i)=qn(i)+f;
    end
end
for i=1:s+1
    qo(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((Z.*(3.*z4+z3))./6).^(alpha2.*(mu+(j-1))));
      qo(i)=qo(i)+f;
    end
end


for i=1:s+1
    qp(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*2.*KCNN.*(1./64).*gamma(alpha2.*(mu+(j-1))))./((z4.*Z.*(2./3)).^(alpha2.*(mu+(j-1))));
      qp(i)=qp(i)+f;
    end
end

summf=0;
for i=1:s+1
    summf(i)=(qa(i)+qb(i)+qc(i)+qd(i)+qe(i)+qf(i)+qg(i)+qh(i)+qi(i)+qj(i)+qk(i)+ql(i)+qm(i)+qn(i)+qo(i)+qp(i));
end
summg=0;
for i=1:s+1
    summg(i)=summd(i)+summe(i)-summf(i);
end
semilogy(z,I,'-r')
hold on;
semilogy(z,summ2,'+k')
hold on;
%HQAM
semilogy(z,I3,':c')
hold on;
semilogy(z,sumf,'+k')
hold on;
semilogy(z_1,I1,'diamondg')
hold on;
semilogy(z,summc,'--b')
hold on;
semilogy(z,summg,'--b')
semilogy(z_1,I4,'diamondg')
%legend('hide')
%legend('Exact(SQAM)','Proposed(SQAM)','Exact(HQAM)','Proposed(HQAM)')
legend('Exact(SQAM)','Proposed(SQAM)','Exact(HQAM)','Proposed(HQAM)','Simulations','Asymptotic')
toc
