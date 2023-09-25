function [ Model, e, eB ] = ALT_FIX_4D_axis( X,ec,index  )
nmax=1e6;  % Maximal mode number
vmax=400;  % Maximal iterations for each mode
ebr=1e-10; % Convergence criterion of each mode

K1=size(X,3);
K2=size(X,4);

%---- Index for the axis center --------
[Mc, index_c(1)]=max(sum(index,2));
[Mc, index_c(2)]=max(sum(index,1));


%---- Initiation -----------------------
f=ndSparse(zeros(size(X,1),size(X,2))); % N-D sparse matrix
bta=ndSparse(zeros(size(X,1),size(X,2)));
bt=ndSparse(zeros(size(X,1),size(X,2)));
   
alpha=1;
B.F1=ones(1,size(X,1))';
B.F2=ones(1,size(X,2))';
B.F3=ones(1,size(X,3))';
B.F4=ones(1,size(X,4))';

Mindex=[];

for k2=1:K2
    for k1=1:K1
        if index(k1,k2)==1
            f(:,:,k1,k2)=zeros(size(X,1),size(X,2));
            norm_X(k1,k2)=norm(full(X(:,:,k1,k2)));
            Mindex=[Mindex ; k1 k2 ];
        end
    end
end

fna=f;
R=f;
nex=0;

%----- (Optional) uncomment following to start from previous computations--
% if exist('calcdata.mat','file')==2       
%     load('calcdata')
%     clear e eB
%     for i=1:size(Mindex,1)
%         k1=Mindex(i,1);
%         k2=Mindex(i,2);        
%         D1=diag(Model.alpha.*Model.F3(k1,:)); 
%         D2=diag(Model.F4(k2,:));                                                               
%         fna(:,:,k1,k2)=Model.F1*D1*D2*Model.F2';
%         norm_X(k1,k2)=norm(full(X(:,:,k1,k2)));
%     end
%     nex=size(Model.alpha,2);
% end

%---- Altenating fix point algorithm --------
for n=(nex+1):nmax    
    for i=1:size(Mindex,1)
        k1=Mindex(i,1);
        k2=Mindex(i,2);   
        R(:,:,k1,k2)=(X(:,:,k1,k2))-(fna(:,:,k1,k2));  
        norm_R(k1,k2)=norm(full(R(:,:,k1,k2)));    
    end
 
    e(n)=max(max(norm_R(index)./norm_X(index))); 
    e(n)  %this displays (check the progress toward the stopping criteria!)
    if e(n)<ec
        break;
    end
    M=0;
    for i=1:size(Mindex,1)
        k1=Mindex(i,1);
        k2=Mindex(i,2);
        D1=diag(B.F3(k1,:)); 
        D2=diag(B.F4(k2,:));          
        M=M+full(R(:,:,k1,k2))*B.F2*D1*D2;  
    end
    B.F1=M/((B.F2'*B.F2)*(B.F3'*B.F3)*(B.F4'*B.F4)*alpha);
    B.F1=B.F1/norm(B.F1);

    %sanity check (writing purposes)
%     M2=0; 
%     for i=1:size(Mindex,1)
% 
%         k1=Mindex(i,1);
%         k2=Mindex(i,2);        
%         D1=diag(B.F3(k1,:)); 
%         D2=diag(B.F4(k2,:));  
%         M2=M2+full(R(:,:,k1,k2))*B.F2/(D1*D2);  
%     end
%     siobhan=M2/((B.F2'*B.F2)*alpha);
%     siobhan=siobhan/norm(siobhan);
% 
%     [B.F1,siobhan]

    M=0;
    for i=1:size(Mindex,1)
        k1=Mindex(i,1);
        k2=Mindex(i,2);
        D1=diag(B.F3(k1,:)); 
        D2=diag(B.F4(k2,:));                            
        M=M+full(R(:,:,k1,k2))'*B.F1*D1*D2; 
    end
    B.F2=M/((B.F1'*B.F1)*(B.F3'*B.F3)*(B.F4'*B.F4)*alpha);
    B.F2=B.F2/norm(B.F2);
    
    R_full=zeros(K1,K2);
    for i=1:size(Mindex,1)
        k1=Mindex(i,1);
        k2=Mindex(i,2);    
        R_full(k1,k2)=B.F1'*full(R(:,:,k1,k2))*B.F2;                                            
    end   

    for k1=1:K1    
        M=0;
        for k2=index_c(2):index_c(2)
            D2=diag(B.F4(k2,:));  
            M=M+R_full(k1,k2)*D2;    
        end
        B.F3(k1,:)=((B.F4'*B.F4)*(B.F2'*B.F2)*(B.F1'*B.F1)*alpha)\diag(M);
    end   
    B.F3=B.F3/norm(B.F3);

%     %sanity check (writing purposes)
%     for k1=1:K1    
%         M2=0;
%         for k2=index_c(2):index_c(2)
% 
%             D2=diag(B.F4(k2,:));  
%             M2=M2+R_full(k1,k2)/D2;    
%         end
%         siobhan(k1,:)=((B.F2'*B.F2)*(B.F1'*B.F1)*alpha)\diag(M2);
%     end   
%     siobhan=siobhan/norm(siobhan);
  
    for k2=1:K2   
        M=0;
        for k1=index_c(1):index_c(1)     
            D1=diag(B.F3(k1,:));
            M=M+R_full(k1,k2)*D1;                         
        end
        B.F4(k2,:)=((B.F3'*B.F3)*(B.F2'*B.F2)*(B.F1'*B.F1)*alpha)\diag(M);
    end  
    B.F4=B.F4/norm(B.F4); 

%     %sanity check (writing purposes)
%     for k2=1:K2   
%         M=0;
%         for k1=index_c(1):index_c(1)     
%             D1=diag(B.F3(k1,:));
%             M=M+R_full(k1,k2)/D1;                         
%         end
%         siobhan(k2,:)=((B.F2'*B.F2)*(B.F1'*B.F1)*alpha)\diag(M);
%     end  
%     siobhan=siobhan/norm(siobhan);
    
    M=0;
    for i=1:size(Mindex,1)
        k1=Mindex(i,1);
        k2=Mindex(i,2);       
        D1=diag(B.F3(k1,:)); 
        D2=diag(B.F4(k2,:));                                                                                         
        M=M+R_full(k1,k2)*D1*D2; 
    end
    alpha=M/((B.F1'*B.F1)*(B.F2'*B.F2)*(B.F3'*B.F3)*(B.F4'*B.F4)) ;
    
    for i=1:size(Mindex,1)
        k1=Mindex(i,1);
        k2=Mindex(i,2);       
        D1=diag(B.F3(k1,:)); 
        D2=diag(B.F4(k2,:));                                    
        bta(:,:,k1,k2)=alpha*B.F1*D1*D2*B.F2';             
    end
       
    for v=1:vmax    
                M=0;
                for i=1:size(Mindex,1)
                    k1=Mindex(i,1);
                    k2=Mindex(i,2);
                    D1=diag(B.F3(k1,:)); 
                    D2=diag(B.F4(k2,:));                            
                    M=M+full(R(:,:,k1,k2))*B.F2*D1*D2;                                              
                end
                B.F1=M/((B.F2'*B.F2)*(B.F3'*B.F3)*(B.F4'*B.F4)*alpha) ;
                B.F1=B.F1/norm(B.F1);

                M=0;
                for i=1:size(Mindex,1)
                    k1=Mindex(i,1);
                    k2=Mindex(i,2);
                    D1=diag(B.F3(k1,:)); 
                    D2=diag(B.F4(k2,:));                            
                    M=M+full(R(:,:,k1,k2))'*B.F1*D1*D2; 
                end
                B.F2=M/((B.F1'*B.F1)*(B.F3'*B.F3)*(B.F4'*B.F4)*alpha) ;
                B.F2=B.F2/norm(B.F2);

                R_full=zeros(K1,K2);
                for i=1:size(Mindex,1)
                    k1=Mindex(i,1);
                    k2=Mindex(i,2);
                    R_full(k1,k2)=B.F1'*full(R(:,:,k1,k2))*B.F2;                                            
                end   

                for k1=1:K1    
                    M=0;
                    for k2=index_c(2):index_c(2)       
                        D2=diag(B.F4(k2,:));                                       
                        M=M+R_full(k1,k2)*D2;                               
                    end
                    B.F3(k1,:)=((B.F4'*B.F4)*(B.F2'*B.F2)*(B.F1'*B.F1)*alpha)\diag(M);
                end   
                B.F3=B.F3/norm(B.F3);

                for k2=1:K2   
                    M=0;
                    for k1=index_c(1):index_c(1)     
                        D1=diag(B.F3(k1,:));
                        M=M+R_full(k1,k2)*D1;                         
                    end
                    B.F4(k2,:)=((B.F3'*B.F3)*(B.F2'*B.F2)*(B.F1'*B.F1)*alpha)\diag(M);
                end  
                B.F4=B.F4/norm(B.F4);                      

                M=0;
                for i=1:size(Mindex,1)
                    k1=Mindex(i,1);
                    k2=Mindex(i,2);
                    D1=diag(B.F3(k1,:)); 
                    D2=diag(B.F4(k2,:));                            
                    M=M+R_full(k1,k2)*D1*D2; 
                end
                alpha=M/((B.F1'*B.F1)*(B.F2'*B.F2)*(B.F3'*B.F3)*(B.F4'*B.F4)) ;

                for i=1:size(Mindex,1)
                    k1=Mindex(i,1);
                    k2=Mindex(i,2);
                    D1=diag(B.F3(k1,:)); 
                    D2=diag(B.F4(k2,:));                            
                    bt(:,:,k1,k2)=alpha*B.F1*D1*D2*B.F2';   
                    be(k1,k2)=norm(full(bt(:,:,k1,k2)-bta(:,:,k1,k2)));
                end
                 
                
                eB(v,n)=max(max(be))/max(max(norm_R)); 
                
           if eB(v,n)<ebr
                f=bt;                                                                              
                if n==1
                       Model.alpha(:,n)=alpha;   
                       Model.F1(:,n)=B.F1;
                       Model.F2(:,n)=B.F2;
                       Model.F3(:,n)=B.F3;
                       Model.F4(:,n)=B.F4;                      
                else
                    %---- Find out if it is a repetitive mode ------
                    for i=size(Model.alpha,2)
                           d_BF(1)=max(abs(B.F1-Model.F1(:,i))/max(abs(Model.F1(:,i))));
                           d_BF(2)=max(abs(B.F2-Model.F2(:,i))/max(abs(Model.F2(:,i)))); 
                           d_BF(3)=max(abs(B.F3-Model.F3(:,i))/max(abs(Model.F3(:,i))));
                           d_BF(4)=max(abs(B.F4-Model.F4(:,i))/max(abs(Model.F4(:,i)))); 
                           
                           if max(d_BF)<1e-2 %repetitive mode                              
                                   Model.alpha(i)=Model.alpha(i)+alpha;                               
                                   break
                           else   % not repetitive mode
                                   Model.alpha=[Model.alpha alpha];   
                                   Model.F1=[Model.F1 B.F1];
                                   Model.F2=[Model.F2 B.F2];
                                   Model.F3=[Model.F3 B.F3];
                                   Model.F4=[Model.F4 B.F4];                                   
                           end
                     end                    
                end
               fna=fna+f;                      

               bta=ndSparse(zeros(size(X,1),size(X,2)));
               bt=bta;

               %--- save compution data ---- 
               if (mod(n,2)==1)
                   save('calcdata','Model','eB','e');
               else 
                   save('calcdata1','Model','eB','e');
               end
               break   
           elseif v==100                   
                   disp(strcat('no convergence at iteration ',num2str(v)))                   
                   ebr=ebr*10;
                   disp(strcat('basis error:',num2str(ebr)))
           elseif v==200
                   disp(strcat('no convergence at iteration ',num2str(v)))                   
                   ebr=ebr*100;
                   disp(strcat('basis error:',num2str(ebr)))
           elseif v==400
                    disp(strcat('no convergence at iteration ',num2str(v),' save data and pass to next mode'))
                    f=bt;                                                                              
                    if n==1
                           Model.alpha(:,n)=alpha;   
                           Model.F1(:,n)=B.F1;
                           Model.F2(:,n)=B.F2;
                           Model.F3(:,n)=B.F3;
                           Model.F4(:,n)=B.F4;                           
                    else
                        %---- Find out if it is a repetitive mode ------
                        for i=size(Model.alpha,2)
                               d_BF(1)=max(abs(B.F1-Model.F1(:,i))/max(abs(Model.F1(:,i))));
                               d_BF(2)=max(abs(B.F2-Model.F2(:,i))/max(abs(Model.F2(:,i)))); 
                               d_BF(3)=max(abs(B.F3-Model.F3(:,i))/max(abs(Model.F3(:,i))));
                               d_BF(4)=max(abs(B.F4-Model.F4(:,i))/max(abs(Model.F4(:,i)))); 
                               
                               if max(d_BF)<1e-2 %repetitive mode                              
                                       Model.alpha(i)=Model.alpha(i)+alpha;                               
                                   break
                               else   % not repetitive mode
                                       Model.alpha=[Model.alpha alpha];   
                                       Model.F1=[Model.F1 B.F1];
                                       Model.F2=[Model.F2 B.F2];
                                       Model.F3=[Model.F3 B.F3];
                                       Model.F4=[Model.F4 B.F4];                                     
                               end
                         end                    
                    end
                   fna=fna+f;                      

                   bta=ndSparse(zeros(size(X,1),size(X,2)));
                   bt=bta;
                 
                   if (mod(n,2)==1)
                       save('calcdata','Model','eB','e');
                   else 
                       save('calcdata1','Model','eB','e');
                   end
                   break  
           elseif isnan(eB(v,n))
                   break  
           else 
                   bta=bt;                   
           end 
                   
            
    end
           
end
    

end




