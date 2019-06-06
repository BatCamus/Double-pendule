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
%         Fnl=calc_Fnl_Double_Pendule(X,dX);
%         Non appel à la fonction pour gagner du temps
                ODE=0;
                Fnl=zeros(size(X));
                %
                c=cos(X(1)-X(2));
                s=sin(X(1)-X(2));

                if(ODE==0)
                Fnl(1)=(m2*l1*c*s*(dX(1))^2+m2*l2*s*(dX(2))^2+m1*g*sin(X(1))+m2*g*s*cos(X(2)))/(m1*l1+m2*l1*(s)^2);
                Fnl(2)=-((m1+m2)*l1*s*dX(1)^2+m2*l2*c*s*dX(2)^2+(m1+m2)*g*c*sin(X(1))-(m1+m2)*g*sin(X(2)))/(m1*l2+m2*l2*s^2);
                end

                % Formule trouvé sur internet
                if(ODE==1)
                Fnl(1)=-((-m2*l1*sin(X(1)-X(2))*cos(X(1)-X(2))*dX(1)^2)-(m2*l2*sin(X(1)-X(2))*dX(2)^2)-(m1*g*sin(X(1)))-(m2*g*sin(X(1)-X(2))*cos(X(2))))/(m1*l1+m2*l1*sin(X(1)-X(2))^2);
                Fnl(2)=-(((m1+m2)*l1*sin(X(1)-X(2))*dX(1)^2)+(m2*l2*sin(X(1)-X(2))*cos(X(1)-X(2))*dX(2)^2)+((m1+m2)*g*sin(X(1)-X(2))*cos(X(1))))/(m1*l2+m2*l2*sin(X(1)-X(2))^2);
                end
                
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
