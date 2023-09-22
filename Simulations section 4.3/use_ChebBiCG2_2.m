%use this script second for generating snapshots in the parameter estimation simulation 

A_of_mu = @(mu) A0 + (sin(mu(1))^2)*A1 + (cos(mu(2))^2)*A2;

%one
N1=7; N2=N1; mus=[]; mus2=[];
a0=.5; b0=1; a1=0; b1=.5;
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 

mu2list=linspace(a1,b1,N2); 
ibound=b1; 
mu1hat = (a0+b0)/2; 
x=chebfun('x');
f1 = (sin(mu1hat))^2;
g = @(x) cos(x)^2;  
g2 = chebfun(g,[-ibound ibound]); 

N=max(length(f1),length(g2)); 
ak2 = chebpoly(g2,N); 
ak=ak2; 

maxit=40; 
mu2s=mu2list; 

tau=.4; 
tol = 1.e-10; 

[n,b,comp_size,c,ctilde,last2,tau_list,A0,A1,A2] = helm_setup2(ak2,N,ibound,tau,f1,g); 
 
Ptau = A0+A1*f1+A2*g(tau); 
[L,U,P,Q]=lu(Ptau); 
applyPtau = @(x) Q*(U\(L\(P*x))); 
applyPtauT = @(x) P'*(L'\(U'\(Q'*x))); 

eyen = speye(n); eyedn = speye(comp_size-n); 
zeros1 = sparse(n,comp_size-n); zeros2 = sparse(comp_size-n,n); 
pimat = [zeros1, eyen; eyedn, zeros2]; 
Pd=A2*ak2(1);
applyPreCon = @(x) pimat*(apply_Uinv(apply_Linv(x,tau,ibound,N,applyPtau,last2,n),tau_list,n)); 
applyPreConT = @(x) apply_LinvT(apply_UinvT(pimat'*x,tau_list,n),tau,ibound,N,applyPtauT,last2,n); 

rxold=c; rzold=ctilde;  
dx=0; dz=0;
kcount = 0; 
xi=1*ones(length(mu2s),1); xiold=1*ones(length(mu2s),1);
alphaold=1; rhoold=1;
dhatz=zeros(comp_size,length(mu2s)); dhatx=zeros(comp_size,length(mu2s)); 
xhat=zeros(comp_size,length(mu2s)); zhat=zeros(comp_size,length(mu2s)); 
res_mat=[];
 
for i=1:maxit
    i
    kcount=0; 
    rho=(rzold')*rxold;
    beta=-rho/rhoold;
    dx=rxold-beta*dx;
    dz=rzold-conj(beta)*dz;
    vec=apply_M(ak,ak2,n,A0,A1,A2,a1,ibound,applyPreCon(dx),N,Pd);
    alpha=rho/((dz')*vec);

    %update the seed system 
    rx=rxold-alpha*(vec);
    rz=rzold-conj(alpha)*(applyPreConT(apply_M(ak,ak2,n,A0,A1,A2,a1,ibound,dz,N,Pd')));

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
        [res,xreturn] = postprocess(xhat(:,k),k,b,applyPreCon,mu2s,A0,A1,A2,f1,g,mu(k)); 
        res_mat(i+1,k)=res; 
        
        if k==length(mu2s)

            %Check residual of 'slowest'
            [res,xreturn] = postprocess(xhat(:,k),k,b,applyPreCon,mu2s,A0,A1,A2,f1,g,mu(k)); 
            res_mat(i+1,k)=res; 

            %If the 'slowest' one converged 
            if res<=tol
                solution_mat2 = xreturn; 
                kcount = kcount + 1; 

                %Check residual of the others
                for l=length(mu2s)-1:-1:1 
                    [res,xreturn] = postprocess(xhat(:,l),l,b,applyPreCon,mu2s,A0,A1,A2,f1,g,mu(l)); 
                    res_mat(i+1,l)=res; 
                    solution_mat2 = [solution_mat2 xreturn]; 
                    
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

function [res,xreturn] = postprocess(xhat,k,b,applyPreCon,mu2s,A0,A1,A2,f1,g,mu)
    A=A0+A1*f1+A2*g(mu2s(k)); 
    n=length(A0); 
    x2(:,k)=applyPreCon(xhat);
    x2long(:,k)=mu*x2(:,k); 
    xreturn = x2long(1:n,k); %xreturn is the corresponding solution 
    res=norm(A*x2long(1:n,k)-b)/norm(b); %returns norm of the rel res
end
