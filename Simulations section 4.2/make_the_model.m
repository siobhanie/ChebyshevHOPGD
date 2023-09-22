%clear all; clc; 
load('nonsparse_helm_eq.mat'); n = length(A0); b = b'; 
% f1 = @(mu) (cos(mu(1)) + mu(1)^2 + sin(mu(2)) + mu(2)); 
% A_of_mu = @(mu) A0 + 2*pi^2*A2 + f1(mu)*A2;

f1 = @(mu) cos(mu(1)) + mu(1)^3;
f2 = @(mu) sin(mu(2)) + mu(2)^2; 
A_of_mu = @(mu) A0 + f1(mu)*A2 + f2(mu)*A3;

% % %---Equidistant nodes 
a0=1; b0=2; 
a1=1; b1=2; 
N1=11; N2=7;
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 
[Mu,Mu2]=meshgrid(mu,mu2); 
Mu = Mu(:); Mu2=Mu2(:);

k=0; 
for j=1:length(mu2)
        
        for i=length(mu):-1:1
            k=k+1; 
            %solvec2 = A_of_mu([mu(i),mu2(j)])\b; 
            solvec = solution_mat1(:,k); 
            norm(solvec2-solvec);
            X(:,i,j) = solvec;
        end
end


% %---Build surrogate model using HOPGD ---------------
ec=1.e-4;
[ Model, e, eB] = HOPGD(X,ec);

%---Compute the relative error the model produces at some of the snapshot locations (mu1,mu2)
X1 = Model.F1; X2 = Model.F2; X3 = Model.F3;
save('X1_full_grid.mat','X1');
save('X2_full_grid.mat','X2');
save('X3_full_grid.mat','X3');




% k=length(X1(1,:)); rel_err=[]; 
% for i=1:N1
%     for j=1:N2
%         model_sol=zeros(n,1); 
%         for l=1:k
%             model_sol = model_sol + X1(:,l)*X2(i,l)*X3(j,l); 
%         end
%         exact_sol = X(:,i,j);
%         model_Tensor(:,i,j) = model_sol; 
%         perErr(i,j) = norm(exact_sol-model_sol)*100/norm(exact_sol); 
%     end
% end
