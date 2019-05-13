function [dFX dFdX] = calc_dFnl_Double_Pendule_numerique(X,dX,t)

global m1 m2 l1 l2 g mu;

dFX=zeros(length(X));
dFdX=zeros(length(X));

dFX(1,1)=