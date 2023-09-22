
dtA1 = -dt*A1; dtA2 = -dt*A2; 
f1 = @(mu1) 1+(sin(mu1)) ; 
f2 = @(mu2) 10 + cos(mu2) + pi*mu2;
A_of_mu = @(mu) A0 + f1(mu(1))*dtA1 + f2(mu(2))*dtA2;


% % %---Equidistant nodes 
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 

o=ones(1,N1); 
mus=[mu,((a0+b0)/2)*o]; mus=mus'; 
mus2=[((a1+b1)/2)*o,mu2]; mus2=mus2'; 
index=zeros(N1,N2); sparsegrid=[]; 
n = length(A0); 

% %---Snapshots (generated using backslash for now)
X_T=ndSparse(zeros(n,1));  


for i=1:N1
    
   
    solvec1 = A_of_mu([mus(i),mus2(i)])\u0;
    solvec = solution_mat1(:,N1-i+1);
    norm(solvec-solvec1)

    X_T(:,:,(N1-1)/2+1,i) = solvec;
    index((N1-1)/2+1,i)=1;
    sparsegrid=[sparsegrid; [mus(i),mus2(i)]]; 

end

for i=1:N2
 
    solvec1 = A_of_mu([mus(end),mus2(i+N2)])\u0;
    solvec = solution_mat2(:,N2-i+1);
    norm(solvec1-solvec)
    
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

save('a_melina_c.mat','a');
save('X1_melina_c.mat','X1');
save('X2_melina_c.mat','X2');
save('X3_melina_c.mat','X3');
save('X4_melina_c.mat','X4');

