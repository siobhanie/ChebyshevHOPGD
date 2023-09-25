mus=[]; mus2=[]; N1=5; N2=5; 
mu=linspace(0,.5,N1); 
mu2=linspace(0,.5,N2); 
load('a_melina.mat'); load('X1_melina.mat'); load('X2_melina.mat'); 
load('X3_melina.mat'); load('X4_melina.mat')

global n2; global a2 ; global X1global; global X2global; global X3global; global X4global
global muglobal; global mu2global; global guess1; 
n2=n; a2=a; 
X1global = X1; X2global = X2; X3global = X3; X4global = X4; 
muglobal = mu; mu2global = mu2; 

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

f1 = @(mu1) 1+(sin(mu1)); 
f2 = @(mu2) 10 + cos(mu2) + pi*mu2;
A_of_mu = @(mu) A0 + f1(mu(1))*dtA1 + f2(mu(2))*dtA2;

xs1 = linspace(.2,.3,20); 
ys1 = linspace(.2,.3,20); 
[Xs,Ys]=meshgrid(xs1,ys1); 
xs1 = Xs(:); ys1=Ys(:); 
n=length(X1(:,1)); 

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
    int_sol = make_int_sol(a,n,X1,X2,vals3(i,:),vals4(i,:));
    approx = [approx int_sol]; 
    A = A_of_mu([xs1(i),ys1(i)]); 
    exact_sol = A\b; 
    sol_storage = [sol_storage exact_sol]; 
    rel_err2(i) = norm(int_sol - exact_sol)*100/norm(exact_sol);
    rel_res(i) = norm(A*int_sol - b)/norm(b); 
end

%save('int_sol.mat','int_sol')
rel_err2=rel_err2'

figure(1)
h=gca;
contourf(Xs,Ys,reshape(rel_err2,length(Xs),length(Ys)))
set(h,'zscale','log')
h.FontSize = 14; 
colorbar
xlabel('$\mu_1$','interpreter','latex','FontSize',24)
ylabel('$\mu_2$','interpreter','latex','FontSize',24)
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [15 15]); %Setthe paper to have width 5 and height 5. 

function approx = make_int_sol(a,n,X1,X2,val1,val2)
    m=length(X1(1,:)); 
    approx=zeros(n,1); 
    for k=1:50   
        approx = approx + a(k)*X1(:,k)*X2(k)*val1(k)*val2(k); %interpolate in both dir.
    end
end

