clc;clf;close all;

load('sparse_helm_eq.mat'); b=b';
load('a_sparse_grid.mat'); load('X1_sparse_grid.mat'); 
load('X2_sparse_grid.mat'); load('X3_sparse_grid.mat'); 
load('X4_sparse_grid.mat')

N1=7; N2=7; n=length(X1(:,1));
a0=1; b0=2; a1=1; b1=2; 
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 

f1 = @(mu) (cos(mu(1)) + mu(1)^4 + sin(mu(2)) + mu(2)); 
A_of_mu = @(mu) A0 + 2*pi^2*A1 + f1(mu)*A1;

xs1 = 1.6;
ys1 = 1.4; 

vals3=[]; vals4=[]; rel_err2=[]; 
%Evaluate the interpolation for each (x,y) of interest 
for i=1:length(xs1)
    for k=1:length(a)
        val3 = interp1(mu,X4(:,k),xs1(i),'spline','extrap');
        vals3(i,k) = val3; 
    end
end
    
for i=1:length(ys1)
    for k=1:length(a)
        val4 = interp1(mu2,X3(:,k),ys1(i),'spline','extrap');
        vals4(i,k) = val4; 
    end    
end

approx = []; sol_storage = []; 
%Compute error at the other locations, using the interpolated model 
for i=1:length(xs1)
    i
    tic;
    int_sol = make_int_sol(a,n,X1,X2,vals3(i,:),vals4(i,:));
    toc; 
    approx = [approx int_sol]; 
    A = A_of_mu([xs1(i),ys1(i)]); 
    tic;
    exact_sol = A\b; 
    toc;
    sol_storage = [sol_storage exact_sol]; 
    rel_err2(i) = norm(int_sol - exact_sol)*100/norm(exact_sol);
    rel_res(i) = norm(A*int_sol - b)/norm(b); 
end

% save('int_sol.mat','int_sol')
% rel_err2=rel_err2'

function approx = make_int_sol(a,n,X1,X2,val1,val2)
    m=length(X1(1,:)); 
    approx=zeros(n,1); 
    for k=1:m   
        approx = approx + a(k)*X1(:,k)*X2(k)*val1(k)*val2(k); %interpolate in both dir.
    end
end

