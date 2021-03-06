function [tt,Xt,dXt]=newmark_Double_Pendule(X0,dX0,t_init,dt,t_tot)

global m1 m2 l1 l2 g mu;
% global knl f0 OMEG
global M NUM

precNR=1.e-11;
% solution initiale
t=t_init;
n=1;
X=X0;dX=dX0; % X et dX sont des matrices (2,1) avec X(1) et X(2) represente les angles theta1 et theta2  
% et dX(1) et dX(2) les vitesses associes

Fnl=zeros(size(X)); % Matrice representant les efforts non lin�aires appliques au pendule
dFX=zeros(length(X)); % Derive de Fnl par rapport a theta1/2
dFdX=zeros(length(X)); % Derive de Fnl par rapport a theta point 1/2
Fnl=calc_Fnl_Double_Pendule(X,dX)
ddX=M\(-Fnl);
tt(n)=t;
Xt(:,n)=X;
dXt(:,n)=dX;

% integration temporelle
for t=t_init+dt:dt:t_tot   %Boucle sur les pas de temps
    n=n+1;
    % prediction des positions, vitesse et acc�l�ration 
    iter=0;
    X=X+dt*dX+(dt^2/2)*ddX;
    dX=dX+dt*ddX;
    ddX=ddX;
    % Calcul du residu
    Fnl=calc_Fnl_Double_Pendule(X,dX); % Calcul des efforts non lineaire
    res=-M*ddX-Fnl;
    normres=norm(res);
    while (normres>precNR)    %Newton Raphson : on r�soud le syst�me pour un pas de temps jusqu'a atteindre une certaine precision
        iter=iter+1;
    
        if (NUM==0)
        % Calcul de la Jacobienne calcule a la main
        [dFX, dFdX]=calc_dFnl_Double_Pendule(X,dX);
        end
        if NUM
        % Calcul de la Jacobienne numeriquement 
        [dFX, dFdX]=calc_dFnl_Double_Pendule_numerique(X,dX,Fnl);
        end

        K_eff=(4/dt^2)*M+(2/dt)*dFdX+dFX; % Matrice de raideur effective 
        % Calcul de la correction : deltaX
        deltaX=K_eff\res;
        % On corrige les positions, vitesses et acceleration avec la
        % corection calculee
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
        % Calcul du residu
        Fnl=calc_Fnl_Double_Pendule(X,dX);
        res=-M*ddX-Fnl;
        normres=norm(deltaX)/norm(X);
    end
    tt(n)=t;
    Xt(:,n)=mod(X+pi,2*pi)-pi;
    dXt(:,n)=dX;
end
