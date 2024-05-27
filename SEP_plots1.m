%Fig.:5(a),5(b),5(c),6 of the manuscript is obtained using this code.
tic
clear all
format long
a=1;
alpha=2;
a_2=alpha/2;
N1=3;
kappa=0;
mu=1;
N=145;
B=1;
%z1=14.92820266532045.*B;
z1=(26.27414291213963).*B;
%z1=1;
z=0:1:30;
t=10.^(z./10);
z_1=0:5:30;
t_1=10.^(z_1./10);
for i=1:31
    %A(i)=(alpha.*(mu.^mu).*(kappa+1).^(mu))./(2.*gamma(mu).*exp(kappa.*mu).*t(i).^(a_2.*mu));
    y=[1];
    for l=1:N-1
        k=(mu.*(1+kappa))./(t(i).^a_2);
        k1=((-k).^l)./factorial(l);
        y=[y k1];
    end
    n(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        p=0;
        s(j)=0;
        for j1=1:N
        
            f1=(vpa(A(j)).*y(j1).*(1./((vpa(z1)).^(a_2.*(mu+(j-1))+p.*a_2))).*gamma(a_2.*(mu+(j-1))+p.*a_2));
            %f1=gamma(a_2*(mu+j+p));
            p=p+1;
            s(j)=s(j)+vpa(f1);
        end
        n(i)=n(i)+vpa(s(j));
    end
end
% semilogy(z,n,'*r')
%z2=2.000000006712126.*B;
z2=(3.239828844061399).*B;
%z2=4/3;

for i=1:31
    %A(i)=(alpha.*mu.^(mu).*(kappa+1).^(mu))./(2.*gamma(mu).*exp(kappa.*mu).*t(i).^(a_2.*mu));
    y=[1];
    for l=1:N-1
        k=(mu.*(1+kappa))./(t(i).^a_2);
        k1=((-k).^l)./factorial(l);
        y=[y k1];
    end
    n1(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        p=0;
        s1(j)=0;
        for j1=1:N
            f1=(vpa(A(j)).*y(j1).*(1./(vpa(z2)).^(a_2.*(mu+(j-1))+p.*a_2)).*gamma(a_2.*(mu+(j-1))+p.*a_2));
            %f1=gamma(a_2*(mu+j+p));
            p=p+1;
            s1(j)=s1(j)+vpa(f1);
        end
        n1(i)=n1(i)+vpa(s1(j));
    end
end
%semilogy(z,s1,'-b')

%z3=1.071796761489149;
z3=(1.446462700182916).*B;
% 
for i=1:31
    %A(i)=(alpha.*mu.^(mu).*(kappa+1).^(mu))./(2.*gamma(mu).*(exp(kappa.*mu)).*(t(i).^(a_2.*mu)))
    y=[1];
    for l=1:N-1
        k=(mu.*(1+kappa))./(t(i).^a_2);
        k1=((-k).^l)./factorial(l);
        y=[y k1];
    end
    n2(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        p=0;
        s2(j)=0;
        for j1=1:N
        
            f1=(vpa(A(j)).*y(j1).*(1./(vpa(z3)).^(a_2.*(mu+(j-1))+p.*a_2)).*gamma(a_2.*(mu+(j-1))+p.*a_2));
            %f1=gamma(a_2*(mu+j+p));
            p=p+1;
            s2(j)=s2(j)+vpa(f1);
        end
        n2(i)=n2(i)+vpa(s2(j));
    end
end
%semilogy(z,s2,'-y')
% hold on;
z4=(1.039566130751374).*B;

for i=1:31
    %A(i)=(alpha.*mu.^(mu).*(kappa+1).^(mu))./(2.*gamma(mu).*(exp(kappa.*mu)).*(t(i).^(a_2.*mu)))
    y=[1];
    for l=1:N-1
        k=(mu.*(1+kappa))./(t(i).^a_2);
        k1=((-k).^l)./factorial(l);
        y=[y k1];
    end
    n3(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        p=0;
        s3(j)=0;
        for j1=1:N
        
            f1=(vpa(A(j)).*y(j1).*(1./(vpa(z4)).^(a_2.*(mu+(j-1))+p.*a_2)).*gamma(a_2.*(mu+(j-1))+p.*a_2));
            %f1=gamma(a_2*(mu+j+p));
            p=p+1;
            s3(j)=s3(j)+vpa(f1);
        end
        n3(i)=n3(i)+vpa(s3(j));
    end
end




% % %semilogy(z,s,'-r')
% 
summ=0;
for i=1:31
    summ(i)= (1/4).*a.*(n(i)+n1(i)+n2(i)+n3(i));
end

%assymptotic
for i=1:31
    %A(i)=(alpha.*(mu.^mu).*(kappa+1).^(mu))./(2.*gamma(mu).*exp(kappa.*mu).*t(i).^(a_2.*mu));
    m(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        f1=(vpa(A(j)).*(1./((vpa(z1)).^(a_2.*(mu+(j-1))))).*gamma(a_2.*(mu+(j-1))));
        m(i)=m(i)+vpa(f1);
    end
end

for i=1:31
    %A(i)=(alpha.*mu.^(mu).*(kappa+1).^(mu))./(2.*gamma(mu).*exp(kappa.*mu).*t(i).^(a_2.*mu));
    m1(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        f1=(vpa(A(j)).*(1./((vpa(z2)).^(a_2.*(mu+(j-1))))).*gamma(a_2.*(mu+(j-1))));
        m1(i)=m1(i)+vpa(f1);
    end
end
for i=1:31
    %A(i)=(alpha.*mu.^(mu).*(kappa+1).^(mu))./(2.*gamma(mu).*(exp(kappa.*mu)).*(t(i).^(a_2.*mu)));
    m2(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        f1=(vpa(A(j)).*(1./((vpa(z3)).^(a_2.*(mu+(j-1))))).*gamma(a_2.*(mu+(j-1))));
        m2(i)=m2(i)+vpa(f1);
    end
end

for i=1:31
    %A(i)=(alpha.*mu.^(mu).*(kappa+1).^(mu))./(2.*gamma(mu).*(exp(kappa.*mu)).*(t(i).^(a_2.*mu)))
    m3(i)=0;
    for j=1:N1+1
        A(j)=(a_2.*(mu.^(mu+2.*(j-1))).*(kappa.^(j-1)).*(kappa+1).^(mu+(j-1)))./(gamma(mu+(j-1)).*factorial(j-1).*exp(kappa.*mu).*t(i).^(a_2.*(mu+(j-1))));
        f1=(vpa(A(j)).*(1./(vpa(z4)).^(a_2.*(mu+(j-1)))).*gamma(a_2.*(mu+(j-1))));
        m3(i)=m3(i)+vpa(f1);
    end
end




% % %semilogy(z,s,'-r')
% 
summ2=0;
for i=1:31
    summ2(i)= (1/4).*(m(i)+m1(i)+m2(i)+m3(i));
end
%semilogy(z,summ,'*b')
% hold on
% 
% 
% 
% 
% 
%exact
syms j x
% alpha=4;
% mu=1;
% kappa=1;
% z=0:1:20;
% t=10.^(z./10);
for i=1:31
    f=@(x,j) (a_2.*a.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((a_2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sqrt(B.*x)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t(i).^((a_2).*(mu+j)));
    f2= symsum(f,j,0,3);
    ht = matlabFunction(f2);
    I(i)=integral(ht,0,100);
    
end
for i=1:7
    f=@(x,j) (a_2.*a.*(mu.^(mu+(2.*j))).*(kappa.^j).*((1+kappa).^(mu+j)).*(x.^(((a_2).*(mu+j))-1)).*exp((-mu.*(1+kappa)./t_1(i).^(alpha./2)).*(x.^(alpha./2))).*erfc(sqrt(B.*x)))./(gamma(mu+j).*factorial(j).*exp(kappa.*mu).*t_1(i).^((a_2).*(mu+j)));
    f2= symsum(f,j,0,3);
    ht = matlabFunction(f2);
    I1(i)=integral(ht,0,10);
    
end
semilogy(z,I,'-r')
hold on;
%semilogy(z,summ,'squarek')
semilogy(z,summ,'+k')
hold on;
semilogy(z_1,I1,'diamondg')
hold on;
semilogy(z,summ2,'--b')
%legend('hide')

legend('Exact','Proposed','Simulation','Asymptotic')
toc