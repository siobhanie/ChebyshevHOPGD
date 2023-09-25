function [ Model, e, eB ] = ALT_FIX_3D( X,ec )
    nmax=30; % Maximal mode number (Note, I changed this from 1000)
    vmax=100; % Maximal iterations for each mode
    ebr=1e-8;  % Convergence criterion of each mode
    ebr0=ebr;
    
    K=size(X,3); 
    Model.F1=[];
    Model.F2=[];
    Model.F3=[];

    f=zeros(size(X));
    bta1=0;
    bta2=0;
    bta3=0;
    B.F1=ones(1,size(X,1))';
    B.F2=ones(1,size(X,2))';
    B.F3=ones(1,size(X,3))';
    fna=0;
    fa=0;
       
    for k=1:K % change to parfor if needed
        norm_X(k)=norm(X(:,:,k));
    end
    norm_X=max(norm_X);
    %---- Altenating fix point algorithm --------
    for n=1:nmax
        R=X-fna;  
        for k=1:K % change to parfor if needed
            norm_R(k)=norm(R(:,:,k));          
        end
        e(n)=norm(norm_R./norm_X);
        norm(norm_R./norm_X)
        if e(n)<ec
            break;
        end
        
        M=0;     
%         for k=1:K 
%             D=diag(B.F3(k,:));
%             M=M+R(:,:,k)*B.F2*D;
%         end
        Rn=reshape(R,[size(R,1)*size(R,3) size(R,2)]);
        RnF2=Rn*B.F2;
        RnF2=reshape(RnF2,[size(R,1)  size(R,3)]);
        RnF2F3=RnF2*B.F3;
        M=RnF2F3;
        
        B.F1=M/((B.F2'*B.F2)*(B.F3'*B.F3)) ;
        B.F1=B.F1/norm(B.F1);

        M=0;
%         for k=1:K 
%             D=diag(B.F3(k,:));
%             M=M+R(:,:,k)'*B.F1*D;      
%         end
        Rn=reshape(R,[size(R,1) size(R,2)*size(R,3)]);
        Rn=Rn';
        RnF1=Rn*B.F1;
        RnF1=reshape(RnF1,[size(R,2)  size(R,3)]);
        RnF1F3=RnF1*B.F3;
        M=RnF1F3; 
        
        B.F2=M/((B.F1'*B.F1)*(B.F3'*B.F3)) ;
        B.F2=B.F2/norm(B.F2);

%         for k=1:K
%             B.F3(k,:)=((B.F2'*B.F2)*(B.F1'*B.F1))\diag(B.F1'*R(:,:,k)*B.F2);
%         end         
        Rn=reshape(R,[size(R,1)*size(R,3) size(R,2)]);
        RnF2=Rn*B.F2;      
        RnF2=reshape(RnF2,[size(R,1) size(R,3)]);
        RnF2F1=B.F1'*RnF2;
        RnF2F1=RnF2F1';
        B.F3=((B.F2'*B.F2)*(B.F1'*B.F1))\RnF2F1;        
        B.F3=B.F3/norm(B.F3);
        
%             br1=norm(B.F1);
%             br2=norm(B.F2);
%             br3=norm(B.F3);
            
        RnF2F1F3=B.F3'*RnF2F1;
        alpha=RnF2F1F3/((B.F3'*B.F3)*(B.F2'*B.F2)*(B.F1'*B.F1));

        for v=1:vmax     
                   M=0;
%                     for k=1:K
%                         D=diag(B.F3(k,:));
%                         M=M+R(:,:,k)*B.F2*D;
%                     end
                    Rn=reshape(R,[size(R,1)*size(R,3) size(R,2)]);
                    RnF2=Rn*B.F2;
                    RnF2=reshape(RnF2,[size(R,1)  size(R,3)]);
                    RnF2F3=RnF2*B.F3;
                    M=RnF2F3;
                    
                    B.F1=M/((B.F2'*B.F2)*(B.F3'*B.F3)) ;
                    B.F1=B.F1/norm(B.F1);
                    

                    M=0;
%                     for k=1:K  
%                         D=diag(B.F3(k,:));
%                         M=M+R(:,:,k)'*B.F1*D;      
%                     end
                    Rn=reshape(R,[size(R,1) size(R,2)*size(R,3)]);
                    Rn=Rn';
                    RnF1=Rn*B.F1;
                    RnF1=reshape(RnF1,[size(R,2)  size(R,3)]);
                    RnF1F3=RnF1*B.F3;
                    M=RnF1F3; 
                    B.F2=M/((B.F1'*B.F1)*(B.F3'*B.F3)) ;
                    B.F2=B.F2/norm(B.F2);

%                     for k=1:K
%                         B.F3(k,:)=((B.F2'*B.F2)*(B.F1'*B.F1))\diag(B.F1'*R(:,:,k)*B.F2);
%                     end 
                    Rn=reshape(R,[size(R,1)*size(R,3) size(R,2)]);
                    RnF2=Rn*B.F2;      
                    RnF2=reshape(RnF2,[size(R,1) size(R,3)]);
                    RnF2F1=B.F1'*RnF2;
                    RnF2F1=RnF2F1';
                    B.F3=((B.F2'*B.F2)*(B.F1'*B.F1))\RnF2F1;
                    B.F3=B.F3/norm(B.F3);
                    
                    RnF2F1F3=B.F3'*RnF2F1;
                    alpha=RnF2F1F3/((B.F3'*B.F3)*(B.F2'*B.F2)*(B.F1'*B.F1));
                    
%                        bt1=norm(B.F1);
%                        bt2=norm(B.F2);
%                        bt3=norm(B.F3);       
%                        eB.e1(v,n)=norm(bt1-bta1)/br1;
%                        eB.e2(v,n)=norm(bt2-bta2)/br2;
%                        eB.e3(v,n)=norm(bt3-bta3)/br3;
                   for k=1:K % change to parfor if needed
                      D=diag(B.F3(k,:));    
                      f(:,:,k)=B.F1*D*B.F2'*alpha;
                   end
                   df=f-fa;
                   eB.df(v,n)=mean(abs(df),'all')/mean(norm_X);
%                if (eB.e1(v,n)<ebr&&eB.e2(v,n)<ebr&&eB.e3(v,n)<ebr)||eB.df(v,n)<ebr
               if  eB.df(v,n)<ebr
%                    for k=1:K
%                    D=diag(B.F3(k,:));    
%                    f(:,:,k)=B.F1*D*B.F2';  
%                    end
                   alpha=alpha^(1/3);
                   Model.F1=[Model.F1 B.F1*alpha];
                   Model.F2=[Model.F2 B.F2*alpha];
                   Model.F3=[Model.F3 B.F3*alpha];
                   fna=fna+f;
%                    bta1=0;
%                    bta2=0;
%                    bta3=0;
                   B.F2=ones(1,size(X,2))';
                   B.F3=ones(1,size(X,3))';
                   ebr=ebr0;
                   break     
               else 
                   fa=f;
%                    bta1=bt1;
%                    bta2=bt2;
%                    bta3=bt3;

               end

        end
              
    end

end

