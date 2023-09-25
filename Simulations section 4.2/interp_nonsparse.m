clf; 

a0=1; b0=2; 
a1=1; b1=2; 
N1=11; N2=7;
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 
[Mu,Mu2]=meshgrid(mu,mu2); 
Mu = Mu(:); Mu2=Mu2(:);

load('nonsparse_helm_eq.mat'); b = b'; 
load('X1_full_grid.mat'); 
load('X2_full_grid.mat'); load('X3_full_grid.mat'); 

f1 = @(mu) cos(mu(1)) + mu(1)^3;
f2 = @(mu) sin(mu(2)) + mu(2)^2; 
A_of_mu = @(mu) A0 + f1(mu)*A2 + f2(mu)*A3;

n = length(X1(:,1)); 
xs1 = linspace(1,2,20); 
ys1 = linspace(1,2,20); 
[Xs,Ys]=meshgrid(xs1,ys1); 
xs1 = Xs(:); ys1=Ys(:); 


vals1=[]; vals2=[]; rel_err2=[]; 
for i=1:length(xs1)
    for k=1:length(X1(1,:))
        val1 = interp1(mu,X2(:,k),xs1(i),'spline','extrap');
        vals1(i,k) = val1; 
    end
end
    
for i=1:length(ys1)
    for k=1:length(X1(1,:))
        val2 = interp1(mu2,X3(:,k),ys1(i),'spline','extrap');
        vals2(i,k) = val2; 
    end    
end

approx = []; sol_storage = []; 
%Compute error at the other locations, using the interpolated model 
for i=1:length(xs1)
    i
    int_sol = make_int_sol(n,X1,vals1(i,:),vals2(i,:));
    approx = [approx int_sol]; 
    A = A_of_mu([xs1(i),ys1(i)]); 
    exact_sol = A\b; 
    sol_storage = [sol_storage exact_sol]; 
    rel_err2(i) = norm(int_sol - exact_sol)*100/norm(exact_sol);
end

% figure(2)
h=gca;
contourf(Xs,Ys,reshape(rel_err2,length(Xs),length(Ys)))
set(h,'zscale','log')
h.FontSize = 14; 
colorbar
xlabel('$\mu_1$','interpreter','latex','FontSize',24)
ylabel('$\mu_2$','interpreter','latex','FontSize',24)
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [15 15]); %Setthe paper to have width 5 and height 5. 

save('int_sol2.mat','int_sol')
rel_err2=rel_err2'

function approx = make_int_sol(n,X1,val1,val2)
    m=length(X1(1,:)); 
    approx=zeros(n,1); 
    for k=1:m   
        approx = approx + X1(:,k)*val1(k)*val2(k); %interpolate in both dir.
    end
end
