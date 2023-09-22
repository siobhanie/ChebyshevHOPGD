function [n,b,comp_size,c,ctilde,last2,tau_list,A0,A1,A2] = helm_setup2(ak2,N,b1,tau,f1,g)

load('param_est.mat')
n = length(b); b=b';  
btilde=randn(n,1); 
comp_size = n*(length(ak2)-1); 
c=zeros(comp_size,1); c(end-n+1:end)=b;
ctilde=zeros(comp_size,1); ctilde(end-n+1:end)=btilde;

ak2flip=fliplr(ak2);

last2=[];
for i=1:length(ak2flip)-1
    if i==1
       	last2=[last2 A2*ak2flip(i+1)]; %prev last i==2
    elseif i==2 
        last2=[last2 A2*ak2flip(i+1)]; %prev last i==3
    elseif i==length(ak2flip)-3 %prev last length - 1
        last2=[last2 A2*ak2flip(end-2)-A2*ak2flip(end)];     
    elseif i==length(ak2flip)-2 %new
        last2=[last2 A2*ak2flip(end-1) + (2*tau/b1)*(A2*ak2flip(end))]; 
    elseif i==length(ak2flip)-1
        last2=[last2 A0+A1*f1+A2*g(tau)];
    else
        last2=[last2 A2*ak2flip(i+1)];
    end  
end

t0 = 1; t1 = (1/b1)*tau; 
tau_list = [t0, t1]; 
for i=2:(N-2)
   tau_list = [tau_list (2/b1)*tau*tau_list(end) - tau_list(end-1)];  
end
tau_list=-tau_list(2:end); 