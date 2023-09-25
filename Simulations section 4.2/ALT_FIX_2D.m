function [ Model, e, eB ] = ALT_FIX_2D( X,ec )
    nmax=1000; % Maximal mode number
    vmax=2400; % Maximal iterations for each mode
    ebr=1e-6;  % Convergence criterion of each mode
    
    f=zeros(size(X));
    bta1=0;
    bta2=0;
    B.F1=ones(1,size(X,1))';
    B.F2=ones(1,size(X,2))';
    fna=0;
       
    norm_X=norm(X(:,:));

    %---- Altenating fix point algorithm --------
    for n=1:nmax
        R=X-fna;  
        norm_R=norm(R(:,:));          
        e(n)=norm(norm_R./norm_X);
        if e(n)<ec
            break;
        end
        M=0;
        M=M+R(:,:)*B.F2;
        B.F1=M/((B.F2'*B.F2)) ;

        M=0;
        M=M+R(:,:)'*B.F1;      
        B.F2=M/((B.F1'*B.F1)) ;
    
        br1=norm(B.F1);
        br2=norm(B.F2);

        for v=1:vmax     
               M=0;
               M=M+R(:,:)*B.F2;
               B.F1=M/((B.F2'*B.F2)) ;

               M=0;
               M=M+R(:,:)'*B.F1;      
               B.F2=M/((B.F1'*B.F1)) ;

               bt1=norm(B.F1);
               bt2=norm(B.F2);     
               eB.e1(v,n)=norm(bt1-bta1)/br1;
               eB.e2(v,n)=norm(bt2-bta2)/br2;
                   
               if eB.e1(v,n)<ebr&&eB.e2(v,n)<ebr   
                   f(:,:)=B.F1*B.F2';                    
                   Model.F1(:,n)=B.F1;
                   Model.F2(:,n)=B.F2;
                   fna=fna+f;
                   bta1=0;
                   bta2=0;
                   B.F2=ones(1,size(X,2))';
                   break     
               else 
                   bta1=bt1;
                   bta2=bt2;
               end

        end
              
    end

end

