clear all; clc; close all; 
rng('default'); rng(1);
%Add path to Chebfun

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

%one
N1=5; N2=N1; mus=[]; mus2=[];
a0=.25; b0=.5; a1=.25; b1=.5;  
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 

mulist=linspace(a0,b0,N1);  
mu2hat = (a1+b1)/2; %Freeze mu2
x=chebfun('x');
g1 = f2(mu2hat);  
f2cheb = chebfun(f1,[-b0 b0]); %Chebyshev inter for func of mu1

N=max(length(f2cheb));
ak = chebpoly(f2cheb,N);  

maxit=100; 
mu2s=mulist; 

tau=.02; 
tol = 1.e-8; 

[n,b,comp_size,c,ctilde,last2,tau_list] = melina_setup1(ak,N,b0,tau,f1,g1,b,A0,dtA1,dtA2); 
 
Ptau = A0+dtA1*f1(tau)+dtA2*g1; 
[L,U,P,Q]=lu(Ptau); 
applyPtau = @(x) Q*(U\(L\(P*x))); 
applyPtauT = @(x) P'*(L'\(U'\(Q'*x))); 

eyen = speye(n); eyedn = speye(comp_size-n); 
zeros1 = sparse(n,comp_size-n); zeros2 = sparse(comp_size-n,n); 
pimat = [zeros1, eyen; eyedn, zeros2]; 
Pd=dtA1*ak(1); ak2=0; 
applyPreCon = @(x) pimat*(apply_Uinv(apply_Linv(x,tau,b0,N,applyPtau,last2,n),tau_list,n)); 
applyPreConT = @(x) apply_LinvT(apply_UinvT(pimat'*x,tau_list,n),tau,b0,N,applyPtauT,last2,n); 

rxold=c; rzold=ctilde;  
dx=0; dz=0;
kcount = 0; 
xi=1*ones(length(mu2s),1); xiold=1*ones(length(mu2s),1);
alphaold=1; rhoold=1;
dhatz=zeros(comp_size,length(mu2s)); dhatx=zeros(comp_size,length(mu2s)); 
xhat=zeros(comp_size,length(mu2s)); zhat=zeros(comp_size,length(mu2s)); 
    
for i=1:maxit
    i
    kcount=0; 
    rho=(rzold')*rxold;
    beta=-rho/rhoold;
    dx=rxold-beta*dx;
    dz=rzold-conj(beta)*dz;
    vec=apply_M(ak,ak2,n,A0,dtA1,dtA2,a0,b0,applyPreCon(dx),N,Pd);
    alpha=rho/((dz')*vec);

    %update the seed system 
    rx=rxold-alpha*(vec);
    rz=rzold-conj(alpha)*(applyPreConT(apply_M(ak,ak2,n,A0,dtA1,dtA2,a0,b0,dz,N,Pd')));

    for k=1:length(mu2s)
        mu(k)=-1./(-mu2s(k)+tau); 
        xinew(k)=(1-alpha*mu(k))*xi(k) + (alpha*beta/alphaold)*(xiold(k)-xi(k));
        alphahat(k)=-alpha*(xi(k)/xinew(k)); 
        betahat(k)=((xiold(k)/xi(k))^2)*beta;
        dhatx(:,k)=(1/xi(k))*rxold - betahat(k)*dhatx(:,k);
        xhat(:,k)=xhat(:,k)+alphahat(k)*dhatx(:,k);

        %%If you are interested in the adjoint solution also
        %dhatz(:,k)=(1/conj(xi(k)))*rzold - conj(beta)*dhatz(:,k);
        %zhat(:,k)=zhat(:,k)+conj(alphahat(k))*dhatz(:,k);
        
        xiold(k)=xi(k); xi(k)=xinew(k); 

        %Note, just for plotting purposes!!! The algorithm is cheaper without this.  
        [res,xreturn] = postprocess(xhat(:,k),k,b,applyPreCon,mu2s,A0,dtA1,dtA2,f1,g1,mu(k)); 
        res_mat(i+1,k)=res; 
        
        if k==length(mu2s)

            %Check residual of 'slowest'
            [res,xreturn] = postprocess(xhat(:,k),k,b,applyPreCon,mu2s,A0,dtA1,dtA2,f1,g1,mu(k)); 
            res_mat(i+1,k)=res; 

            %If the 'slowest' one converged 
            if res<=tol
                solution_mat1 = xreturn; 
                kcount = kcount + 1; 

                %Check residual of the others
                for l=length(mu2s)-1:-1:1 
                    [res,xreturn] = postprocess(xhat(:,l),l,b,applyPreCon,mu2s,A0,dtA1,dtA2,f1,g1,mu(l)); 
                    res_mat(i+1,l)=res; 
                    solution_mat1 = [solution_mat1 xreturn]; 
                    
                    if res<=tol
                        kcount = kcount + 1; 
                        
                        %If solutions have all converged, quit 
                        if kcount == length(mu2s)
                            res_mat(1,1)=1; res_mat(1,2)=1; res_mat(1,3)=1; res_mat(1,4)=1;
                            res_mat(1,5)=1; res_mat(1,6)=1; res_mat(1,7)=1; 
                            figure(1)
                            semilogy(res_mat(:,1))
                            hold on
                            semilogy(res_mat(:,2))
                            semilogy(res_mat(:,3))
                            semilogy(res_mat(:,4))
                            semilogy(res_mat(:,5))
                            semilogy(res_mat(:,6))
                            semilogy(res_mat(:,7))
                            xlabel('its')
                            ylabel('rel res')
                            return; 
                        end
                    
                    end
                
                end
            
            end        
        
        end
    
    end
    alphaold=alpha; 
    rxold=rx; rzold=rz; 
    rhoold=rho;
end

plot(solution_mat1(:,1))
hold on
plot(solution_mat1(:,2))
plot(solution_mat1(:,3))
plot(solution_mat1(:,4))
plot(solution_mat1(:,5))
plot(solution_mat1(:,6))
plot(solution_mat1(:,7))


function [res,xreturn] = postprocess(xhat,k,b,applyPreCon,mu2s,A0,dtA1,dtA2,f1,g1,mu)
    A=A0+dtA1*f1(mu2s(k))+dtA2*g1; 
    n=length(A0); 
    x2(:,k)=applyPreCon(xhat);
    x2long(:,k)=mu*x2(:,k); 
    xreturn = x2long(1:n,k); %xreturn is the corresponding solution 
    res=norm(A*x2long(1:n,k)-b)/norm(b); %returns norm of the rel res
end
