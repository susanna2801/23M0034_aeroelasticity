clc
clear; 

m = 3.3843;
S_theta = 0.0859;
I_theta = 0.0135;
c = 0.2540;
rho = 1.225;
K_h = 2818.8;
K_theta = 37.3;
a_mat=(-0.5:0.05:0)';
%a=0;
b= c/2;



r_theta = sqrt(I_theta/(m*b^2));
x_theta = S_theta/(m*b);
mu = m/(rho*pi*b^2);
w_h = sqrt(K_h/m);
w_theta = sqrt(K_theta/I_theta);
R=sqrt(w_h^2/w_theta^2);
syms Omega
k_initial=0.1;
k_step=0.01;
k_final=1.5;
k_array=(k_initial:k_step:k_final)';
for p=1:size(a_mat)
    a=a_mat(p)
    for j=1:size(k_array)
        k=k_array(j);
        Lh=-1+2i/k;
        Ltheta=(1/k)*(1i+2)+a+(2/k^2);
        Mh=0.5;
        Mtheta=(3/8)-(1i/k);
        lh=Lh;
        ltheta=Ltheta-(0.5+a)*Lh;
        mh=Mh+(0.5+a)*Lh;
        mtheta=Mtheta-(0.5+a)*(-Ltheta+Mh)-((0.5+a)^2)*Lh;
    
        D1=mu*(1-R^2*Omega)-lh;
        D2=mu*x_theta-ltheta;
        D3=mu*x_theta+mh;
        D4=mu*r_theta^2*(1-Omega)+mtheta;
    
        Det(j) = D1*D4-D2*D3;
    
        Root = vpa(solve(Det(j),Omega));
        RootR(j,:) =real(Root);
        RootI(j,:)=imag(Root);
        g(j,:)=imag(Root)./real(Root);
        w(j,:)=w_theta./sqrt(real(Root));
        U(j,:)=b*w(j,:)./k;
     end
     figure(1)
     plot(U,g)
     xlabel('Velocity(U)')
     ylabel('g')
     title('U-g plot')
     hold on
     figure(2)
     plot(U,w)
     xlabel('Velocity(U)')
     ylabel('frequency')
     title('U-freq plot')
     hold on
     

end

   