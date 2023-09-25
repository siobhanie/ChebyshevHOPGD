function [ Model, e, eB ] = ALT_FIX_4D( X,ec )
        nmax=1000; % Maximal mode number
        vmax=2400; % Maximal iterations for each mode
        ebr=1e-6;  % Convergence criterion of each mode
        
        K1=size(X,3);
        K2=size(X,4);

        f=zeros(size(X));
        bta1=0;
        bta2=0;
        bta3=0;
        bta4=0;

        B.F1=ones(1,size(X,1))';
        B.F2=ones(1,size(X,2))';
        B.F3=ones(1,size(X,3))';
        B.F4=ones(1,size(X,4))';
        fna=0;
        
        for k2=1:K2
            for k1=1:K1
            norm_X(k1,k2)=norm(X(:,:,k1,k2));
            index(k1,k2) = ~all(norm_X(k1,k2)==0,2); % Non-zero values
            end
        end
        %---- Altenating fix point algorithm --------
        for n=1:nmax
            R=X-fna;  
            for k2=1:K2
                for k1=1:K1
                norm_R(k1,k2)=norm(R(:,:,k1,k2));
                end
            end
            e(n)=max(max(norm_R(index)./norm_X(index)));
            if e(n)<ec
                break;
            end
            M=0;
            for k2=1:K2
                D2=diag(B.F4(k2,:));
                for k1=1:K1
                    if index(k1,k2)==1
                        D1=diag(B.F3(k1,:));       
                        M=M+R(:,:,k1,k2)*B.F2*D1*D2;
                    end
                end
            end
            B.F1=M/((B.F2'*B.F2)*(B.F3'*B.F3)*(B.F4'*B.F4)) ;
            M=0;
            for k2=1:K2
                D2=diag(B.F4(k2,:));
                for k1=1:K1 
                    if index(k1,k2)==1
                        D1=diag(B.F3(k1,:));       
                        M=M+R(:,:,k1,k2)'*B.F1*D1*D2; 
                    end
                end
            end
            B.F2=M/((B.F1'*B.F1)*(B.F3'*B.F3)*(B.F4'*B.F4)) ;
            for k1=1:K1    
                M=0;
                for k2=1:K2      
                    if index(k1,k2)==1
                        D2=diag(B.F4(k2,:));
                        M=M+B.F1'*R(:,:,k1,k2)*B.F2*D2;   
                    end
                end
                B.F3(k1,:)=((B.F4'*B.F4)*(B.F2'*B.F2)*(B.F1'*B.F1))\diag(M);
            end   
            for k2=1:K2   
                M=0;
                for k1=1:K1      
                    if index(k1,k2)==1
                        D1=diag(B.F3(k1,:));
                        M=M+B.F1'*R(:,:,k1,k2)*B.F2*D1;        
                    end
                end
                B.F4(k2,:)=((B.F3'*B.F3)*(B.F2'*B.F2)*(B.F1'*B.F1))\diag(M);
            end  
                br1=norm(B.F1);
                br2=norm(B.F2);
                br3=norm(B.F3);
                br4=norm(B.F4);
            for v=1:vmax
                       M=0;
                        for k2=1:K2
                            D2=diag(B.F4(k2,:));
                            for k1=1:K1
                                if index(k1,k2)==1
                                    D1=diag(B.F3(k1,:));                
                                    M=M+R(:,:,k1,k2)*B.F2*D1*D2;
                                end
                            end
                        end
                        B.F1=M/((B.F2'*B.F2)*(B.F3'*B.F3)*(B.F4'*B.F4)) ;

                        M=0;
                        for k2=1:K2
                            D2=diag(B.F4(k2,:));
                            for k1=1:K1 
                                if index(k1,k2)==1
                                    D1=diag(B.F3(k1,:));                  
                                    M=M+R(:,:,k1,k2)'*B.F1*D1*D2; 
                                end
                            end
                        end
                        B.F2=M/((B.F1'*B.F1)*(B.F3'*B.F3)*(B.F4'*B.F4)) ;

                        for k1=1:K1    
                            M=0;
                            for k2=1:K2      
                                if index(k1,k2)==1
                                    D2=diag(B.F4(k2,:));
                                    M=M+B.F1'*R(:,:,k1,k2)*B.F2*D2;    
                                end
                            end
                            B.F3(k1,:)=((B.F4'*B.F4)*(B.F2'*B.F2)*(B.F1'*B.F1))\diag(M);
                        end   
                        for k2=1:K2   
                            M=0;
                            for k1=1:K1     
                                if index(k1,k2)==1
                                    D1=diag(B.F3(k1,:));
                                    M=M+B.F1'*R(:,:,k1,k2)*B.F2*D1;  
                                end
                            end
                            B.F4(k2,:)=((B.F3'*B.F3)*(B.F2'*B.F2)*(B.F1'*B.F1))\diag(M);
                        end  
                           bt1=norm(B.F1);
                           bt2=norm(B.F2);
                           bt3=norm(B.F3);     
                           bt4=norm(B.F4); 
                           eB.e1(v,n)=norm(bt1-bta1)/br1;
                           eB.e2(v,n)=norm(bt2-bta2)/br2;
                           eB.e3(v,n)=norm(bt3-bta3)/br3;
                           eB.e4(v,n)=norm(bt4-bta4)/br4;
                   if eB.e1(v,n)<ebr&&eB.e2(v,n)<ebr&&eB.e3(v,n)<ebr&&eB.e4(v,n)<ebr
                       for k2=1:K2  
                           D2=diag(B.F4(k2,:)); 
                           for k1=1:K1
                               D1=diag(B.F3(k1,:));                     
                               f(:,:,k1,k2)=B.F1*D1*D2*B.F2'; 
                           end
                       end
                       Model.F1(:,n)=B.F1;
                       Model.F2(:,n)=B.F2;
                       Model.F3(:,n)=B.F3;
                       Model.F4(:,n)=B.F4;
                       fna=fna+f;
                       bta1=0;
                       bta2=0;
                       bta3=0;
                       bta4=0;
                       B.F2=ones(1,size(X,2))';
                       B.F3=ones(1,size(X,3))';
                       B.F4=ones(1,size(X,4))';
                       break     
                   else 
                       bta1=bt1;
                       bta2=bt2;
                       bta3=bt3;
                       bta4=bt4;
                   end
            end
        end


end

