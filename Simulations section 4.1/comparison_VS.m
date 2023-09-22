clear all; 
clc; close all; 

%A_of_mu = @(mu) 1 + mu(1) + mu(2) + mu(1)*mu(2) + (mu(1)^2)*(mu(2))^2; 
A_of_mu = @(mu) (1+mu(1)+mu(1)^2+mu(1)^3+mu(1)^4+mu(1)^5)*(1+mu(2)+mu(2)^2+mu(2)^3+mu(2)^4+mu(2)^5);   


% % %---Equidistant nodes 
a0=0; b0=1; a1=0; b1=1; N1=15; N2=15; n=10;
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 

o=ones(1,N1); 
mus=[mu,((a0+b0)/2)*o]; mus=mus'; 
mus2=[((a1+b1)/2)*o,mu2]; mus2=mus2'; 
index=zeros(N1,N2); sparsegrid=[]; 


% %---Snapshots (generated using backslash for now)
X_T=ndSparse(zeros(n,1));  


for i=1:N1

    %A_of_mu([mus(i),mus2(i)])
    X_T(:,:,(N1-1)/2+1,i) = ones(n,1)*A_of_mu([mus(i),mus2(i)]);
    index((N1-1)/2+1,i)=1;
    sparsegrid=[sparsegrid; [mus(i),mus2(i)]]; 

end

for i=1:N2

    X_T(:,:,i,(N1-1)/2+1) = ones(n,1)*A_of_mu([mus(end),mus2(i+N2)]);
    index(i,(N1-1)/2+1)=1;
    sparsegrid=[sparsegrid; [mus(end),mus2(i+N2)]];

end


%---Build surrogate model using HOPGD ---------------
ec=1.e-10;
index=logical(index); 
[ Model, e, eB] = HOPGD_axis( X_T,ec,index);


%---Compute the relative error the model produces at some of the snapshot locations (mu1,mu2)
%---Sanity check :) 
X1 = Model.F1; X2 = Model.F2; X3 = Model.F3; X4 = Model.F4; a = Model.alpha; 
k=length(a); rel_err=[]; 

for i=1:N1
    model_sol=zeros(n,1); 
    for j=1:k
        model_sol = model_sol + a(j)*X1(:,j)*X2(j)*X4(i,j)*X3((N1-1)/2+1,j); 
    end
    exact_sol = A_of_mu([mus(i),mus2(i)]); 
    rel_err(i)=norm(exact_sol - model_sol)/norm(exact_sol); 
end
for i=1:N2
    model_sol=zeros(n,1); 
    for j=1:k
        model_sol = model_sol + a(j)*X1(:,j)*X2(j)*X4((N1-1)/2+1,j)*X3(i,j); 
    end
    exact_sol = A_of_mu([mus(i+N1),mus2(i+N1)]); 
    rel_err(i+N1)=norm(exact_sol - model_sol)/norm(exact_sol); 
end

rel_err=rel_err'; 