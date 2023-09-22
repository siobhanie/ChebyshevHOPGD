clear; clc; clf; 

load('param_est.mat'); b=b';
load('a.mat'); load('X1.mat'); load('X2.mat'); load('X3.mat'); load('X4.mat')

N1=7; N2=N1; 
mus=[]; mus2=[];
a0=0; b0=1; a1=0; b1=1; 
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2);
n=length(X1(:,1)); 
A_of_mu = @(mu) A0 + (sin(mu(1))^2)*A1 + (cos(mu(2))^2)*A2;

global n2
global a2
global X1global
global X2global
global X3global
global X4global
global muglobal
global mu2global 
global guess1
n2=n; 
a2=a; 
X1global = X1;
X2global = X2; 
X3global = X3; 
X4global = X4; 
muglobal = mu; 
mu2global = mu2; 

xmine = .87; 
ymine = .12; 


A = A_of_mu([xmine,ymine]); 
guess = A\b;  
guess1 = guess; 

%minpunkt = fmincon(@make_int_sol,[0 0],[],[],[],[],[a0 a1],[b0 b1]);
minpunkt = fminsearch(@make_int_sol,[1,0]);

[xmine,ymine]
[minpunkt]

relerr1 = norm(xmine-minpunkt(1))*100/norm(xmine)
relerr2 = norm(ymine-minpunkt(2))*100 /norm(ymine)

figure(1)
hold on
plot(mu,.5*ones(7,1),'ro')
plot(.5*ones(7,1),mu,'ro')
plot(xmine,ymine,'r*')
plot(minpunkt(1),minpunkt(2),'bo')
xlabel('$\mu_1$','interpreter','latex')
ylabel('$\mu_2$','interpreter','latex')
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [15 15]); %Setthe paper to have width 5 and height 5. 


function approx1 = make_int_sol(val)
    global n2
    global a2
    global X1global
    global X2global
    global X3global
    global X4global
    global muglobal 
    global mu2global 
    global guess1
    %global ys1

    for k=1:length(a2)
        val3(k) = interp1(muglobal,X4global(:,k),val(1),'spline');
        %val3(k) = interp1(muglobal,X4global(:,k),xs1);
    end
    
    for k=1:length(a2)
        val4(k) = interp1(mu2global,X3global(:,k),val(2),'spline');
        %val4(k) = interp1(mu2global,X3global(:,k),ys1);
    end

    m=length(X1global(1,:)); 
    approx=zeros(n2,1); 
    for k=1:m
        approx = approx + a2(k)*X1global(:,k)*X2global(k)*val3(k)*val4(k); %interpolate in both dir.
    end

    approx1 = norm(approx-guess1);
end

