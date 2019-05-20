function [tt,Xt,dXt]=newmark_Double_Pendule(X0,dX0,t_init,dt,t_tot)

global m1 m2 l1 l2 g mu;
% global knl f0 OMEG
global M NUM

precNR=1.e-11;
% sol init
t=t_init;
n=1;
X=X0;dX=dX0; % X et dX sont des matrices (2,1) avec X(1) et X(2) représente les angles theta1 et theta2  
%% et dX(1) et dX(2) les vitesses associés

Fnl=zeros(size(X))
dFX=zeros(length(X))
dFdX=zeros(length(X))
P=zeros(size(X))
%P=f0*cos(OMEG*t);
Fnl=calc_Fnl_Double_Pendule(X,dX)
ddX=M\(P-Fnl);
tt(n)=t;
Xt(:,n)=X;
dXt(:,n)=dX;

% integration temporelle
for t=t_init+dt:dt:t_tot;   %Boucle sur les pas de temps
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
    tt(n)=t;
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
