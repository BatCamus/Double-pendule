function [tt,Xt,dXt,Tr]=Temps_de_retournement_N(X0,dX0,t_init,dt)

global m1 m2 l1 l2 g mu;
% global knl f0 OMEG
global M NUM

precNR=1.e-5;
% sol init
Tr=t_init;
n=1;
X=X0;dX=dX0; % X et dX sont des matrices (2,1) avec X(1) et X(2) représente les angles theta1 et theta2  
%% et dX(1) et dX(2) les vitesses associés

% Fnl=zeros(size(X));
dFX=zeros(length(X));
dFdX=zeros(length(X));
Fnl=calc_Fnl_Double_Pendule(X,dX);
ddX=M\(-Fnl);
tt(n)=t;
Xt(:,n)=X;
dXt(:,n)=dX;

y1=-l1*cos(X(1,1));
y2=y1-l2*cos(X(2,1));

w=sqrt(l1/g);
% integration temporelle
while y2<0 && Tr<w*1e3   %Boucle sur les pas de temps
    n=n+1;
    % prediction
    iter=0;
    X=X+dt*dX+(dt^2/2)*ddX;
    dX=dX+dt*ddX;
    ddX=ddX;
    % Calcul du residu
    P=0;
%   P=f0*cos(OMEG*t);
    Fnl=calc_Fnl_Double_Pendule(X,dX);
    res=P-M*ddX-Fnl;
    normres=norm(res);
    while (normres>precNR);    %Newton Raphson
        iter=iter+1;
    
        if(NUM==0)
        % Calcul de la Jacobienne calculé à la main
        [dFX, dFdX]=calc_dFnl_Double_Pendule(X,dX);
        end
        if(NUM==1)
        % Calcul de la Jacobienne numériquement 
        [dFX, dFdX]=calc_dFnl_Double_Pendule_numerique(X,dX,Fnl);
        end
        
        K_eff=(4/dt^2)*M+(2/dt)*dFdX+dFX;
        % Calcul de la correction
        deltaX=K_eff\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
        % Calcul du residu
         Fnl=calc_Fnl_Double_Pendule(X,dX);
        res=P-M*ddX-Fnl;
%         normres=norm(res)/norm(P);
        normres=norm(deltaX)/norm(X);
    end
    Tr=Tr+dt;
    tt(n)=t;
    Xt(:,n)=X;
    dXt(:,n)=dX;
    y1=-l1*cos(X(1,1));
    y2=y1-l2*cos(X(2,1));
end
