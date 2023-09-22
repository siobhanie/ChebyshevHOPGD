function [ Model, e, eB] = HOPGD( X,ec ) 
%%%
% HOPGD MODEL - Standard implementation
% ----------------------INPUT----------------------
% X       The n-dimensional input array, which can be up to 4D in the
%         present version, eg. X(p1,p2,p3,p4)
% ec      Trancation error, eg. ec=0.1 (percentage)
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
        case 2
        [ Model, e, eB ] = ALT_FIX_2D( X,ec ); 
        case 3
        [ Model, e, eB ] = ALT_FIX_3D( X,ec );            
        case 4
        [ Model, e, eB ] = ALT_FIX_4D( X,ec );            
    end  
end

