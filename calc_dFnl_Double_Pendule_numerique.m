function [dFX, dFdX] = calc_dFnl_Double_Pendule_numerique(X,dX,Fnl)

global m1 m2 l1 l2 g delta;

dFX=zeros(length(X));
dFdX=zeros(length(X));
% On stock au d�but du calcul les variables que l'on a rentr� car dans la
% suite on va changer leurs valeurs et on souhaiterai les garder
X_stock=zeros(length(X));
dX_stock=zeros(length(X));
X_stock(1)=X(1);
X_stock(2)=X(2);
dX_stock(1)=dX(1);
dX_stock(2)=dX(2);

% On utilise ensuite la formule d'int�gration d'ordre 1 et de pr�cision
% d'ordre h : (f(x+h)-f(x))/h. Pour se faire on recalcule Fnl grace � notre
% fonction calc_Fnl mais en changeant les variables comme ci dessous
X(1)=X(1)+delta;
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl1X1=Fnl_temp(1);

X(1)=X_stock(1);
X(2)=X(2)+delta;
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl1X2=Fnl_temp(1);

X(1)=X(1)+delta;
X(2)=X_stock(2);
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl2X1=Fnl_temp(2);

X(1)=X_stock(1);
X(2)=X(2)+delta;
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl2X2=Fnl_temp(2);

X(2)=X_stock(2);
dX(1)=dX(1)+delta;
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl1dX1=Fnl_temp(1);

dX(1)=dX_stock(1);
dX(2)=dX(2)+delta;
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl1dX2=Fnl_temp(1);

dX(1)=dX(1)+delta;
dX(2)=dX_stock(2);
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl2dX1=Fnl_temp(2);

dX(1)=dX_stock(1);
dX(2)=dX(2)+delta;
Fnl_temp=calc_Fnl_Double_Pendule(X,dX);
deltaFnl2dX2=Fnl_temp(2);

% Application des formules d'int�gration num�rique finale 
dFX(1,1)=(deltaFnl1X1-Fnl(1))/delta;
dFX(1,2)=(deltaFnl1X2-Fnl(1))/delta;
dFX(2,1)=(deltaFnl2X1-Fnl(2))/delta;
dFX(2,2)=(deltaFnl2X2-Fnl(2))/delta;

dFdX(1,1)=(deltaFnl1dX1-Fnl(1))/delta;
dFdX(1,2)=(deltaFnl1dX2-Fnl(1))/delta;
dFdX(2,1)=(deltaFnl2dX1-Fnl(2))/delta;
dFdX(2,2)=(deltaFnl2dX2-Fnl(2))/delta;