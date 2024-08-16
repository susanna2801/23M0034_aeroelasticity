clc
clear; 
m = 3.3843;
S_theta = 0.0859;
I_theta = 0.0135;
c = 0.2540;
rho = 1.225;
K_h = 2818.8;
K_theta = 37.3;
    
b= c/2;
r_theta = sqrt(I_theta/(m*b^2));
x_theta = S_theta/(m*b);
mu = m/(rho*pi*b^2);
w_h = sqrt(K_h/m);
w_theta = sqrt(K_theta/I_theta);
w_ratio=sqrt(w_h^2/w_theta^2);
a_mat= (-0.5:0.1:0)';
RootR=[];

syms Omega

vbar_ini=0;
vbar_step=0.1;
vbar_fin=8;

vbar_mat=(vbar_ini:vbar_step:vbar_fin)';

for s=1:size(a_mat)
    a=a_mat(s);
    RootR=real_root_calculation(a,vbar_mat,Omega,r_theta,x_theta,mu,w_ratio);
    real_is_not_zero= 0;
    for i= 1: size(RootR)
        if RootR(i) > 0
            real_is_not_zero= i;
            break;
        end
    end
    Flutter_speed = double(vbar_mat(real_is_not_zero))*b*w_theta

    disp(a)
    %c=['k*','g*','m*','c*','b*','r*'];
    plot(vbar_mat,RootR,'k*')
    xlabel('ND Airpeed')
    ylabel('Real part of ND Frequency')
    title('Model 1 (Flutter speed) for a1')
    if s<6
        hold on
    else
        hold off
    end    
end  


function RootR=real_root_calculation(a,vbar_mat,Omega,r_theta,x_theta,mu,w_ratio)

    for j = 1:size(vbar_mat)
        vbar= vbar_mat(j);
    
        D1 = Omega^2*(1+(1/mu))+w_ratio^2+Omega*(2*vbar/mu);
        D2 = Omega^2*(x_theta-(a/mu))+(2*vbar^2/mu)+Omega*((vbar/mu)+(2*vbar/mu)*(0.5-a)) ;
        D3 = Omega^2*(x_theta-(a/mu))-Omega*((2*vbar/mu)*(a+0.5));
        D4 = Omega^2*(r_theta^2+((0.125+a^2)/mu))+r_theta^2-(2*(a+0.5)*vbar^2)/mu-Omega*(((0.5-a)*vbar/mu)+(2*vbar*(a+0.5)*(0.5-a)/mu));
    
        Det(j) = D1*D4-D2*D3;
        double(vbar);
        Root = vpa(solve(Det(j),Omega));
        RootR(j,:) = real(Root);
  
    end


end    


