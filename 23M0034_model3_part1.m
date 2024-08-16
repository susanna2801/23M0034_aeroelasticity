clc
clear; 

m = 3.3843;
S_theta = 0.0859;
I_theta = 0.0135;
c = 0.2540;
rho = 1.225;
K_h = 2818.8;
K_theta = 37.3;

a=0;
b= c/2;

%change a and K values to find flutter speed


r_theta = sqrt(I_theta/(m*b^2));
x_theta = S_theta/(m*b);
mu = m/(rho*pi*b^2);
w_h = sqrt(K_h/m);
w_theta = sqrt(K_theta/I_theta);
R=sqrt(w_h^2/w_theta^2);
Kbar=100;
su=Kbar/m*b*w_theta^2;
syms Omega
k_initial=0.1;
k_step=0.001;
k_final=0.2;
k_array=(k_initial:k_step:k_final)';
for j=1:size(k_array)
    k=k_array(j);
    C_k=(besselh(1,2,k))./(besselh(1,2,k)+1i.*besselh(0,2,k));
    Lh=-1+2i*C_k/k;
    Ltheta=(1/k)*(1i+2*C_k)+a+(2*C_k/k^2);
    Mh=0.5;
    Mtheta=(3/8)-(1i/k);
    lh=Lh;
    ltheta=Ltheta-(0.5+a)*Lh;
    mh=Mh+(0.5+a)*Lh;
    mtheta=Mtheta-(0.5+a)*(-Ltheta+Mh)-((0.5+a)^2)*Lh;
    
    D1=mu*(1-R^2*Omega)-lh;
    D2=mu*x_theta-ltheta+su*Omega;
    D3=mu*x_theta+mh+su*Omega;
    D4=mu*r_theta^2*(1-Omega)+mtheta;
    
    Det(j) = D1*D4-D2*D3;
    
    Root = vpa(solve(Det(j),Omega));
    RootR(j,:) =real(Root);
    RootI(j,:)=imag(Root);
    g(j,:)=imag(Root)./real(Root)
    w(j,:)=w_theta./sqrt(real(Root));
    U(j,:)=b*w(j,:)./k;
end
when_g_Zero_index = 0;
for iter = 1:size(k_array)
    
    sh(iter,2)=round(g(iter,2),3);
    if sh(iter,2)==0
        when_g_Zero_index=iter;
        break;
    end
end
disp(when_g_Zero_index)
s=w_theta/sqrt(RootR(when_g_Zero_index));
p=round(s,3)
U=p*b/k_array(when_g_Zero_index)