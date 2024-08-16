clc
clear; 

m = 3.3843;
S_theta = 0.0859;
I_theta = 0.0135;
c = 0.2540;
rho = 1.225;
K_h = 2818.8;
K_theta = 37.3;
a = (-0.5:0.05:0)';
b= c/2;

r_theta = sqrt(I_theta/(m*b^2));
x_theta = S_theta/(m*b);
mu = m/(rho*pi*b^2);
w_h = sqrt(K_h/m);
w_theta = sqrt(K_theta/I_theta);
w_ratio=sqrt(w_h^2/w_theta^2);

syms Omega

vbar_ini=0;
vbar_step=0.01;
vbar_fin=4;

vbar_mat=(vbar_ini:vbar_step:vbar_fin)';

for s=1:size(a)
    a1=a(s);
    for j = 1:size(vbar_mat)
        vbar= vbar_mat(j);
    
        D1 = Omega^2+w_ratio^2;
        D2 = Omega^2*(x_theta)+2*vbar^2/mu ;
        D3 = Omega^2*(x_theta);
        D4 = Omega^2*(r_theta^2)+r_theta^2-(2*(a(s)+0.5)*vbar^2)/mu;
    
        Det(j) = D1*D4-D2*D3;
        double(vbar);
        Root = vpa(solve(Det(j),Omega));
        RootR(j,:) = real(Root);
  
    end
    real_not_zero= 0;
    for i= 1: size(RootR)
        if RootR(i)~=0
            real_not_zero= i;
            break;
        end
    end
    Flutter_speed = double(vbar_mat(real_not_zero))*b*w_theta
    disp(a1)
    plot(vbar_mat,RootR,'k*')
    xlabel('ND Airpeed')
    ylabel('Real part of ND Frequency')
    title('Model 1 Flutter speed')
    hold on   
end  




