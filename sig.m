function [sigma1, sigma2, sigma3, sigma4, sigma5, mult] = sig(sigma1, sigma2, sigma3, sigma4, sigma5, r, s)
L=length(r)/5; 
mu=5; 
t_max = 2; 
t_inc=2; 
t_dec=2; 
% Adaptive t_inc , t_dec: 
t=sqrt(norm(r)/norm(s)); 
if ( (1 <= t) && (t < t_max) )
    t_inc = t;   
elseif((1/t_max < t) && (t < 1))
     t_inc = 1/t;
else
  t_inc = t_max;   
end 

t_dec = t_inc; 
mult=[1 1 1 1 1]; 
%% sigma 1

if( norm(r(1:L,1)) > mu*norm(s) )
    sigma1 = sigma1*t_inc; 
    mult(1)=1/t_inc; 
elseif( norm(s) > mu*norm(r(1:L,1)))
    sigma1=sigma1/t_dec;
    mult(1)=t_dec;
end 

%% sigma 2

if( norm(r(1+L:2*L,1)) > mu*norm(s) )
    sigma2 = sigma2*t_inc; 
    mult(2)=1/t_inc;
elseif( norm(s) > mu*norm(r(1+L:2*L,1)))
    sigma2=sigma2/t_dec;
    mult(2)=t_dec;
end 

%% sigma 3

if( norm(r(1+2*L:3*L,1)) > mu*norm(s) )
    sigma3 = sigma3*t_inc; 
    mult(3)=1/t_inc;
elseif( norm(s) > mu*norm(r(1+2*L:3*L,1)))
    sigma3=sigma3/t_dec;
    mult(3)=t_dec;
end 

%% sigma 4

if( norm(r(1+3*L:4*L,1)) > mu*norm(s) )
    sigma4 = sigma4*t_inc; 
    mult(4)=1/t_inc;
elseif( norm(s) > mu*norm(r(1+3*L:4*L,1)))
    sigma4=sigma4/t_dec;
    mult(4)=t_dec;
end 

%% sigma 5

if( norm(r(1+4*L:5*L,1)) > mu*norm(s) )
    sigma5 = sigma5*t_inc; 
    mult(5)=1/t_inc;
elseif( norm(s) > mu*norm(r(1+4*L:5*L,1)))
    sigma5=sigma5/t_dec;
    mult(5)=t_dec;
end 


end

