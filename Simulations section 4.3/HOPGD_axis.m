function [ Model, e, eB] = HOPGD_axis( X,ec,index )
%%%
% HOPGD MODEL - Sparse implementation 
% ----------------------INPUT----------------------
% X       The n-dimensional sparse array, which can be up to 5D in the
%         present version, eg. X(p1,p2,p3,p4,p5)
%         The first two dimension should be the natural dimension of X
%         eg. Given a vector X1 associted with parameters p3*,p4*,p5*, use
%         the following to assigne the components of X
%         X(:,:,p3*,p4*,p5*)=X1;
%         Note: present version only accounts for the components in the
%         axes of the parameter space [p3 p4 p5]
%         eg. if the axis center is p3*,p4*,p5*, one only changes one parameter
%         at a time to get the axes.
% ec      Trancation error, eg. ec=0.1 (percentage)
% index   The logical indicator matrix. eg. if a component of X is known
%         for p3*,p4*,p5*, index(p3*,p4*,p5*)=logical(1), 
%         Otherwise, index(p3*,p4*,p5*)=logical(0)
% ----------------------OUTPUT---------------------
% Model   Separated modes
% e       Evolution of trancation error with respect to mode number
% eB      Convergence of each mode
% -------------------------------------------------
% Copyright (C) 2020  Ye Lu
% Northwestern University, Evanston, Illinois, US, ye.lu@northwestern.edu
%%%
    nd=length(size(X)); % Matrix dimension
    switch nd         
        case 4
        [ Model, e, eB ] = ALT_FIX_4D_axis( X,ec,index  );            
        case 5
        [ Model, e, eB ] = ALT_FIX_5D_axis( X,ec,index  );       
    end
    
end

