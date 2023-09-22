mus=[]; mus2=[]; a0=.25; a1=a0; b0=.5; b1=b0; 
mu=linspace(a0,b0,N1); 
mu2=linspace(a1,b1,N2); 
load('a_melina_b.mat'); load('X1_melina_b.mat'); load('X2_melina_b.mat'); 
load('X3_melina_b.mat'); load('X4_melina_b.mat')

global n2; global a2 ; global X1global; global X2global; global X3global; global X4global
global muglobal; global mu2global; global guess1; 
n2=n; a2=a; 
X1global = X1; X2global = X2; X3global = X3; X4global = X4; 
muglobal = mu; mu2global = mu2; 

f1 = @(mu1) 1+(sin(mu1)); 
f2 = @(mu2) 10 + cos(mu2) + pi*mu2;
A_of_mu = @(mu) A0 + f1(mu(1))*dtA1 + f2(mu(2))*dtA2;

minpunkt = fmincon(@make_int_sol,[minpunkt],[],[],[],[],[a0 a1],[b0 b1]);
relerr1 = norm(xmine-minpunkt(1))*100/norm(xmine);
relerr2 = norm(ymine-minpunkt(2))*100 /norm(ymine);
sqrt(relerr1^2 + relerr2^2);

norm(guess-make_int_sol2(minpunkt))*100/norm(guess)

figure(2)
hold on
plot(linspace(0,.5,5),.25*ones(1,5),'ro')
plot(.25*ones(1,5),linspace(0,.5,5),'ro')
plot(mu,.375*ones(1,5),'ro')
plot(.375*ones(1,5),mu2,'ro')
plot(xmine,ymine,'*')
plot(minpunkt(1),minpunkt(2),'bo')


% 
% data1 = [mu',.375*ones(1,5)'];
% data2 = [.375*ones(1,5)',mu2'];
% data3 = [xmine,ymine];
% data4 = [minpunkt(1),minpunkt(2)];
% 
% 
% header = {};
% header{1} = 'column 1'; header{2} = 'column 2'; 
% header = strjoin(header, ',');
% 
% fid = fopen('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest5_melina.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest5_melina.csv',data1,'-append');
% 
% 
% fid = fopen('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest6_melina.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest6_melina.csv',data2,'-append');
% 
% 
% fid = fopen('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest7_melina.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest7_melina.csv',data3,'-append');
% 
% 
% fid = fopen('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest8_melina.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/experiment/doc/raw data/parest8_melina.csv',data4,'-append');
% 



function approx1 = make_int_sol(val)
    global n2; global a2; global X1global; global X2global; global X3global; global X4global
    global muglobal; global mu2global; global guess1
    
    for k=1:length(a2)
        val3(k) = interp1(muglobal,X4global(:,k),val(1),'spline','extrap');
    end
    
    for k=1:length(a2)
        val4(k) = interp1(mu2global,X3global(:,k),val(2),'spline','extrap');
    end

    m=length(X1global(1,:)); 
    approx=zeros(n2,1); 
    for k=1:25
        approx = approx + a2(k)*X1global(:,k)*X2global(k)*val3(k)*val4(k); %interpolate in both dir.
    end

    approx1 = norm(approx-guess1);
end

function approx = make_int_sol2(val)
    global n2; global a2; global X1global; global X2global; global X3global; global X4global
    global muglobal; global mu2global; global guess1
    
    for k=1:length(a2)
        val3(k) = interp1(muglobal,X4global(:,k),val(1),'spline','extrap');
    end
    
    for k=1:length(a2)
        val4(k) = interp1(mu2global,X3global(:,k),val(2),'spline','extrap');
    end

    m=length(X1global(1,:)); 
    approx=zeros(n2,1); 
    for k=1:25
        approx = approx + a2(k)*X1global(:,k)*X2global(k)*val3(k)*val4(k); %interpolate in both dir.
    end
end

