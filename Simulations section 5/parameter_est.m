clear all; clc; clf; close all; 
rng('default'); rng(1);

x_n = 1; x_0 = 0; n=10000; 
h = (x_n - x_0)/n; 
t0 = 0; T = 1; 
dt = .01; 
xs = [h:h:x_n-h]';

%initial condition 
u0 = sin(pi*xs); 
b = u0; 

diag1 = -(2/h^2)*ones(n-1,1); 
diag2 = (1/h^2)*ones(n-1,1); 
A1 = spdiags([diag2 diag1 diag2],-1:1,n-1,n-1); 

diag1b = (1/h)*ones(n-1,1); 
diag2b = -(1/h)*ones(n-1,1); 
A2 = spdiags([diag2b,diag1b],-1:0,n-1,n-1); 

dtA1 = -dt*A1; dtA2 = -dt*A2; 

f1 = @(mu1) 1+(sin(mu1)) ; 
f2 = @(mu2) 10 + cos(mu2) + pi*mu2;
A0 = speye(n-1); 
A_of_mu = @(mu) A0 + f1(mu(1))*dtA1 + f2(mu(2))*dtA2;

N1=5; N2=5;
a0=0; b0=.5; a1=0; b1=.5;  
mus=[]; mus2=[];
mu=linspace(0,0.5,N1); 
mu2=linspace(0,0.5,N2); 
load('a_melina.mat'); load('X1_melina.mat'); load('X2_melina.mat'); 
load('X3_melina.mat'); load('X4_melina.mat')

global n2; global a2 ; global X1global; global X2global; global X3global; global X4global
global muglobal; global mu2global; global guess1; 

X1global = X1; X2global = X2; X3global = X3; X4global = X4; n=length(X1(:,1)); 
muglobal = mu; mu2global = mu2; 
n2=n; a2=a; 

xmine = .4265; 
ymine = .3; 
f1 = @(mu1) 1+(sin(mu1)); 
f2 = @(mu2) 10 + cos(mu2) + pi*mu2;
A_of_mu = @(mu) A0 + f1(mu(1))*dtA1 + f2(mu(2))*dtA2;

%Perturb with random noise 
A = A_of_mu([xmine,ymine]); 
guess = A\u0;  
delta = (1.e-2)*randn(n,1); 
guess1 = guess + delta; 

minpunkt = fmincon(@make_int_sol,[.5 .5],[],[],[],[],[a0 a1],[b0 b1]);
relerr1 = norm(xmine-minpunkt(1))*100/norm(xmine);
relerr2 = norm(ymine-minpunkt(2))*100/norm(ymine);
sqrt(relerr1^2 + relerr2^2);

norm(guess-make_int_sol2(minpunkt))*100/norm(guess);


figure(1)
hold on
plot(mu,.25*ones(1,5),'ro')
plot(.25*ones(1,5),mu2,'ro')
plot(xmine,ymine,'*')
plot(minpunkt(1),minpunkt(2),'bo')







function approx1 = make_int_sol(val)
    global n2; global a2; global X1global; global X2global; global X3global; global X4global
    global muglobal; global mu2global; global guess1
    
    for k=1:length(a2)
        val3(k) = interp1(muglobal,X4global(:,k),val(1),'spline','extrap');
    end
    
    for k=1:length(a2)
        val4(k) = interp1(mu2global,X3global(:,k),val(2),'spline','extrap');
    end

    m=length(X1global(1,:)); 
    approx=zeros(n2,1); 
    for k=1:50
        approx = approx + a2(k)*X1global(:,k)*X2global(k)*val3(k)*val4(k); %interpolate in both dir.
    end

    approx1 = norm(approx-guess1);
end

function approx = make_int_sol2(val)
    global n2; global a2; global X1global; global X2global; global X3global; global X4global
    global muglobal; global mu2global; global guess1
    
    for k=1:length(a2)
        val3(k) = interp1(muglobal,X4global(:,k),val(1),'spline','extrap');
    end
    
    for k=1:length(a2)
        val4(k) = interp1(mu2global,X3global(:,k),val(2),'spline','extrap');
    end

    m=length(X1global(1,:)); 
    approx=zeros(n2,1); 
    for k=1:50
        approx = approx + a2(k)*X1global(:,k)*X2global(k)*val3(k)*val4(k); %interpolate in both dir.
    end
end


