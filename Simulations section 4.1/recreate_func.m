%Seems like the best choice for accuracy and speed 

%clear all; 
clc; %close all; 
addpath '/Users/siobhanie/Dropbox/InterpolationProject/working codes/chebfun'
load('sparse_helm_eq.mat'); n = length(A0); b = b'; 

A_of_mu = @(mu) A0 + (2*pi^2 + cos(mu(1)) + mu(1)^4 + sin(mu(2)) + mu(2))*A1; 

% % %---Equidistant nodes 
a0=1; b0=2; a1=1; b1=2; 
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 

o=ones(1,N1); 
mus=[mu,((a0+b0)/2)*o]; mus=mus'; 
mus2=[((a1+b1)/2)*o,mu2]; mus2=mus2'; 
index=zeros(N1,N2); sparsegrid=[]; 


% %---Snapshots (generated using backslash for now)
X_T=ndSparse(zeros(n,1));  


for i=1:N1

    %solvec2 = A_of_mu([mus(i),mus2(i)])\b; 
    solvec = solution_mat1(:,N1-i+1);
    X_T(:,:,(N1-1)/2+1,i) = solvec;
    index((N1-1)/2+1,i)=1;
    sparsegrid=[sparsegrid; [mus(i),mus2(i)]]; 

end

for i=1:N2

    %solvec2 = A_of_mu([mus(end),mus2(i+N2)])\b;
    solvec = solution_mat2(:,N2-i+1);
    X_T(:,:,i,(N1-1)/2+1) = solvec;
    index(i,(N1-1)/2+1)=1;
    sparsegrid=[sparsegrid; [mus(end),mus2(i+N2)]];

end


%---Build surrogate model using HOPGD ---------------
ec=1.e-4;
index=logical(index); 
[ Model, e, eB] = HOPGD_axis( X_T,ec,index);


%---Compute the relative error the model produces at some of the snapshot locations (mu1,mu2)
%---Sanity check :) 
X1 = Model.F1; X2 = Model.F2; X3 = Model.F3; X4 = Model.F4; a = Model.alpha; 
k=length(a); rel_err=[]; 
save('a_sparse_grid.mat','a');
save('X1_sparse_grid.mat','X1');
save('X2_sparse_grid.mat','X2');
save('X3_sparse_grid.mat','X3');
save('X4_sparse_grid.mat','X4');



% for i=1:N1
%     model_sol=zeros(n,1); 
%     for j=1:k
%         model_sol = model_sol + a(j)*X1(:,j)*X2(j)*X4(i,j)*X3((N1-1)/2+1,j); 
%     end
%     %A = A_of_mu([mus(i),mus2(i)]); 
%     %exact_sol = A\b; 
%     exact_sol = f1([mus(i),mus2(i)]); 
%     rel_err(i)=norm(exact_sol - model_sol)/norm(exact_sol); 
% end
% for i=1:N2
%     model_sol=zeros(n,1); 
%     for j=1:k
%         model_sol = model_sol + a(j)*X1(:,j)*X2(j)*X4((N1-1)/2+1,j)*X3(i,j); 
%     end
%     %A = A_of_mu([mus(i+N1),mus2(i+N1)]); 
%     %exact_sol = A\b; 
%     exact_sol = f1([mus(i+N1),mus2(i+N1)]); 
%     rel_err(i+N1)=norm(exact_sol - model_sol)/norm(exact_sol); 
% end
% 
% rel_err=rel_err'; 
