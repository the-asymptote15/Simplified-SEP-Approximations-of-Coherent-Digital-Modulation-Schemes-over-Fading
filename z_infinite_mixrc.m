%Fig.:7 is obtained using this code
tic
clear all;
alpha=1;
alpha2=alpha./2;
kappa=0;
mu=2;
M=8;
N2=4;
N=171;
N1=4;
T=sqrt(96/(31*M*N2-32));
S=sqrt(12./((5.*M.*N2)-4));
A=4-(2./M)-(2./N2);
%A=4.*(1-(1./sqrt(M)));
z=0:1:30;
t=10.^(z./10);
z_1=0:5:30;
t1=10.^(z_1./10);
z1=26.27414291213963;
%RQAM
for i=1:31
    na(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z1.*S.^2)./2).^(alpha2.*(mu+(j-1))));
      na(i)=na(i)+f;
    end
end

z2=3.239828844061399;
for i=1:31
    nb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z2.*S.^2)./2).^(alpha2.*(mu+(j-1))));
      nb(i)=nb(i)+f;
    end
end

z3=1.446462700182916;
for i=1:31
    nc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z3.*S.^2)./2).^(alpha2.*(mu+(j-1))));
      nc(i)=nc(i)+f;
    end
end

z4=1.039566130751374;
for i=1:31
    nd(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z4.*S.^2)./2).^(alpha2.*(mu+(j-1))));
      nd(i)=nd(i)+f;
    end
end
summ=0;
for i=1:31
    summa(i)=(na(i)+nb(i)+nc(i)+nd(i));
end


%2nd term

A1=4.*(1-(1./M)).*(1-(1./N2));
%z1=26.27414291213963;
for i=1:31
    ma(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z1.*S.^2)).^(alpha2.*(mu+(j-1))));
      ma(i)=ma(i)+f;
    end
end

%z2=3.239828844061399;
for i=1:31
    mb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z2.*S.^2)).^(alpha2.*(mu+(j-1))));
      mb(i)=mb(i)+f;
    end
end

%z3=1.446462700182916;
for i=1:31
    mc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z3.*S.^2)).^(alpha2.*(mu+(j-1))));
      mc(i)=mc(i)+f;
    end
end

%z4=1.039566130751374;
for i=1:31
    md(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z4.*S.^2)).^(alpha2.*(mu+(j-1))));
      md(i)=md(i)+f;
    end
end

for i=1:31
    me(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z2).*(S.^2)./2)).^(alpha2.*(mu+(j-1))));
      me(i)=me(i)+f;
    end
end


for i=1:31
    mf(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z3+z4).*(S.^2)./2)).^(alpha2.*(mu+(j-1))));
      mf(i)=mf(i)+f;
    end
end

for i=1:31
    mg(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z3).*(S.^2)./2)).^(alpha2.*(mu+(j-1))));
      mg(i)=mg(i)+f;
    end
end

for i=1:31
    
    mh(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z3).*(S.^2)./2)).^(alpha2.*(mu+(j-1))));
      mh(i)=mh(i)+f;
    end
end

for i=1:31
    mi(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z4).*(S.^2)./2)).^(alpha2.*(mu+(j-1))));
      mi(i)=mi(i)+f;
    end
end

for i=1:31
    mj(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z4).*(S.^2)./2)).^(alpha2.*(mu+(j-1))));
      mj(i)=mj(i)+f;
    end
end
summb=0;
for i=1:31
    summb(i)=(ma(i)+mb(i)+mc(i)+md(i)+me(i)+mf(i)+mg(i)+mh(i)+mi(i)+mj(i));
end
summc=0;
for i=1:31
    summc(i)=summa(i)-summb(i);
end



%CQAM
for i=1:31
    oa(i)=0;
    for j=1:N1+1
        B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
        f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z1.*T.^2)./2).^(alpha2.*(mu+(j-1))));
        oa(i)=oa(i)+f;
    end
end

z2=3.239828844061399;
for i=1:31
    ob(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z2.*T.^2)./2).^(alpha2.*(mu+(j-1))));
      ob(i)=ob(i)+f;
    end
end

z3=1.446462700182916;
for i=1:31
    oc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z3.*T.^2)./2).^(alpha2.*(mu+(j-1))));
      oc(i)=oc(i)+f;
    end
end

z4=1.039566130751374;
for i=1:31
    od(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
       f=(B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1))))./(((z4.*T.^2)./2).^(alpha2.*(mu+(j-1))));

      od(i)=od(i)+f;
    end
end
summd=0;
for i=1:31
    summd(i)=(oa(i)+ob(i)+oc(i)+od(i));
end


%2nd term

A1=4-(4/M)-(4/N2)+(8/(M*N2));
%z1=26.27414291213963;
for i=1:31
    qa(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z1.*T.^2)).^(alpha2.*(mu+(j-1))));
      qa(i)=qa(i)+f;
    end
end

%z2=3.239828844061399;
for i=1:31
    qb(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z2.*T.^2)).^(alpha2.*(mu+(j-1))));

      qb(i)=qb(i)+f;
    end
end

%z3=1.446462700182916;
for i=1:31
    qc(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z3.*T.^2)).^(alpha2.*(mu+(j-1))));

      qc(i)=qc(i)+f;
    end
end

%z4=1.039566130751374;
for i=1:31

    qd(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1))))./(((z4.*T.^2)).^(alpha2.*(mu+(j-1))));
      qd(i)=qd(i)+f;
    end
end

for i=1:31

    qe(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z2).*(T.^2)./2)).^(alpha2.*(mu+(j-1))));

      qe(i)=qe(i)+f;
    end
end


for i=1:31

    qf(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z3+z4).*(T.^2)./2)).^(alpha2.*(mu+(j-1))));

      qf(i)=qf(i)+f;
    end
end

for i=1:31

    qg(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z3).*(T.^2)./2)).^(alpha2.*(mu+(j-1))));
 
      qg(i)=qg(i)+f;
    end
end

for i=1:31
    qh(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z3).*(T.^2)./2)).^(alpha2.*(mu+(j-1))));

      qh(i)=qh(i)+f;
    end
end

for i=1:31

    qi(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z1+z4).*(T.^2)./2)).^(alpha2.*(mu+(j-1))));
    qi(i)=qi(i)+f;
    end
end

for i=1:31

    qj(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      f=(B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1))))./((((z2+z4).*(T.^2)./2)).^(alpha2.*(mu+(j-1))));
       
      qj(i)=qj(i)+f;
    end
end
summe=0;
for i=1:31
    summe(i)=(qa(i)+qb(i)+qc(i)+qd(i)+qe(i)+qf(i)+qg(i)+qh(i)+qi(i)+qj(i));
end
summf=0;
for i=1:31
    summf(i)=summd(i)-summe(i);
end

%RQAM
for i=1:31
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
          f(l)=(y(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*S.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n1(i)=n1(i)+s1(j);
    end
end


for i=1:31
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
          f(l)=(y(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*S.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n2(i)=n2(i)+s1(j);
    end
end
for i=1:31
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
          f(l)=(y(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*S.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n3(i)=n3(i)+s1(j);
    end
end

for i=1:31
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
          f(l)=(y(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*S.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      n4(i)=n4(i)+s1(j);
    end
end
summ=0;
for i=1:31
    summ(i)=(n4(i)+n1(i)+n2(i)+n3(i));
end


%2nd term

A1=4.*(1-(1./M)).*(1-(1./N2));
%z1=26.27414291213963;
for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*S.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m1(i)=m1(i)+s1(j);
    end
end

%z2=3.239828844061399;
for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*S.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m2(i)=m2(i)+s1(j);
    end
end

%z3=1.446462700182916;
for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*S.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m3(i)=m3(i)+s1(j);
    end
end

%z4=1.039566130751374;
for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*S.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m4(i)=m4(i)+s1(j);
    end
end

for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z2).*(S.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m5(i)=m5(i)+s1(j);
    end
end


for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z3+z4).*(S.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m6(i)=m6(i)+s1(j);
    end
end

for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z3).*(S.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m7(i)=m7(i)+s1(j);
    end
end

for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z3).*(S.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m8(i)=m8(i)+s1(j);
    end
end

for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z4).*(S.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m9(i)=m9(i)+s1(j);
    end
end

for i=1:31
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
          f(l)=(y(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z4).*(S.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      m10(i)=m10(i)+s1(j);
    end
end
summ1=0;
for i=1:31
    summ1(i)=(m1(i)+m2(i)+m3(i)+m4(i)+m5(i)+m6(i)+m7(i)+m8(i)+m9(i)+m10(i));
end
summ2=0;
for i=1:31
    summ2(i)=summ(i)-summ1(i);
end
syms j x
% for i=1:21
%     f=@(x,j) (alpha2.*(1./2).*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sigma.*sqrt(x./2)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
%     f2= symsum(f,j,0,10);
%     ht = matlabFunction(f2);
%     I(i)=integral(ht,0,10);
% end
for i=1:31
    %f=@(x,j) ((1./4).*alpha2.*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)))).*((A.*(1./2).*erfc(S.*sqrt(x./2)))-(A1.*(1./4).*(erfc(S.*sqrt(x./2))).^2));
    f2= symsum(f,j,0,10);
    ht = matlabFunction(f2);
    I1(i)=integral(ht,0,100);
    
end
syms j x
for i=1:7
    %f=@(x,j) ((1./4).*alpha2.*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t1(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t1(i).^((alpha2).*(mu+j)))).*((A.*(1./2).*erfc(S.*sqrt(x./2)))-(A1.*(1./4).*(erfc(S.*sqrt(x./2))).^2));
    f2= symsum(f,j,0,10);
    ht = matlabFunction(f2);
    I2(i)=integral(ht,0,100);
    
end
% I2=0;
% for i=1:21
%     I2(i)=I(i)-I1(i);
% end


%CQAM
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    o1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f=(y1(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*T.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f;
      end
      o1(i)=o1(i)+s1(j);
    end
end

z2=3.239828844061399;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    o2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*T.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      o2(i)=o2(i)+s1(j);
    end
end

z3=1.446462700182916;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    o3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*T.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      o3(i)=o3(i)+s1(j);
    end
end

z4=1.039566130751374;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    o4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A.*(1./8).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*T.^2)./2).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      o4(i)=o4(i)+s1(j);
    end
end
summ3=0;
for i=1:31
    summ3(i)=(o4(i)+o1(i)+o2(i)+o3(i));
end


%2nd term

A1=4-(4/M)-(4/N2)+(8/(M*N2));
%z1=26.27414291213963;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q1(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z1.*T.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q1(i)=q1(i)+s1(j);
    end
end

%z2=3.239828844061399;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q2(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z2.*T.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q2(i)=q2(i)+s1(j);
    end
end

%z3=1.446462700182916;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q3(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z3.*T.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q3(i)=q3(i)+s1(j);
    end
end

%z4=1.039566130751374;
for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q4(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./64).*gamma(alpha2.*(mu+(j-1)+p)))./(((z4.*T.^2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q4(i)=q4(i)+s1(j);
    end
end

for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q5(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z2).*(T.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q5(i)=q5(i)+s1(j);
    end
end


for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q6(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z3+z4).*(T.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q6(i)=q6(i)+s1(j);
    end
end

for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q7(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z3).*(T.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q7(i)=q7(i)+s1(j);
    end
end

for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q8(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z3).*(T.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q8(i)=q8(i)+s1(j);
    end
end

for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q9(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z1+z4).*(T.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q9(i)=q9(i)+s1(j);
    end
end

for i=1:31
    y1=[1];
    for k=1:N-1
        coeff=(-mu.*(1+kappa))./(t(i).^(alpha./2));
        coeff1=(coeff.^k)./factorial(k);
        y1=[y1,coeff1];
    end
    q10(i)=0;
    for j=1:N1+1
      B(j)=(alpha2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(alpha2.*(mu+(j-1))));
      p=0;
      s1(j)=0;
      for l=1:N
          f(l)=(y1(l).*B(j).*A1.*(1./32).*gamma(alpha2.*(mu+(j-1)+p)))./((((z2+z4).*(T.^2)./2)).^(alpha2.*(mu+(j-1)+p)));
          p=p+1;
          s1(j)=s1(j)+f(l);
      end
      q10(i)=q10(i)+s1(j);
    end
end
summ4=0;
for i=1:31
    summ4(i)=(q1(i)+q2(i)+q3(i)+q4(i)+q5(i)+q6(i)+q7(i)+q8(i)+q9(i)+q10(i));
end
summ5=0;
for i=1:31
    summ5(i)=summ3(i)-summ4(i);
end
syms j x
% for i=1:31
%     f=@(x,j) (alpha2.*(1./2).*A.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sigma.*sqrt(x./2)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
%     f2= symsum(f,j,0,10);
%     ht = matlabFunction(f2);
%     I(i)=integral(ht,0,10);
% end
for i=1:31
    %f=@(x,j) ((1./4).*alpha2.*A1.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)))).*((A.*(1./2).*erfc(T.*sqrt(x./2)))-(A1.*(1./4).*(erfc(T.*sqrt(x./2))).^2));
    f2= symsum(f,j,0,100);
    ht = matlabFunction(f2);
    I3(i)=integral(ht,0,1000);
    
end
for i=1:7
    %f=@(x,j) ((1./4).*alpha2.*A1.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*(erfc(sigma.*sqrt(x./2))).^2)./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((alpha2).*(mu+j)));
    f=@(x,j) ((alpha2.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((alpha2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t1(i).^(alpha./2)).*(x.^(alpha./2))))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t1(i).^((alpha2).*(mu+j)))).*((A.*(1./2).*erfc(T.*sqrt(x./2)))-(A1.*(1./4).*(erfc(T.*sqrt(x./2))).^2));
    f2= symsum(f,j,0,100);
    ht = matlabFunction(f2);
    I4(i)=integral(ht,0,1000);
    
end
% I2=0;
% for i=1:31
%     I2(i)=I(i)-I1(i);
% end
%RQAM
semilogy(z,I1,'-r')
hold on;
semilogy(z,summ2,'+k')
hold on;

%CQAM
semilogy(z,I3,':c')
hold on;
semilogy(z,summ5,'+k')
hold on;
semilogy(z_1,I2,'diamondg')

hold on;
semilogy(z,summc,'--b')
hold on;
semilogy(z,summf,'--b')
hold on;
semilogy(z_1,I4,'diamondg')
%legend('hide')
legend('Exact(RQAM)','Proposed(RQAM)','Exact(XQAM)','Proposed(XQAM)','Simulations','Asymptotic')
toc

toc
