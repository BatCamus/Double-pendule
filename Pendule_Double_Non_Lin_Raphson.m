%% pendule double non linéaire methode Newton Raphson

clear all
close all
global m1 m2 l1 l2 g mu

%% Declaration variable
g = 9.81;         % gravitÃ© terrestre
m1 = 1;           % masse du pendule 1
m2 = 0.9;           % masse du pendule 2
l1 = 1;           % longueur du pendule 1
l2 = 0.9;           % longueur du pendule 2
theta10 =0.8;      % angle formÃ© par le pendule 1 avec la verticale
theta20 = 0.1;        % angle formÃ© par le pendule 2 avec la verticale
theta10p= 0;         %vitesse angulaire initiale du pendule 1
theta20p= 0;         % vitesse angulaire initiale du pendule 1
theta10pp = 0;     % accÃ©lÃ©ration angulaire initiale du pendule 1
theta20pp = 0;     % accÃ©lÃ©ration angulaire initiale du pendule 2
mu = m2/m1;       % rapport des masses : utile pour simplifier l'Ã©quation


w1 = sqrt((g*(1+mu)*(l1+l2)+g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
w2 = sqrt((g*(1+mu)*(l1+l2)-g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
A1 = (1+mu)/mu-l1*w1^2/(mu*g);
A2 = (1+mu)/mu-l1*w2^2/(mu*g);
C1 = (theta20-A2*theta10)/(A1-A2);
C2 = (A1*theta10-theta20)/(A1-A2);
phi1 = asin((theta20p-A2*theta10p)/(C1*w1*(A2-A1)));
phi2 = asin((A1*theta10p-theta20p)/(C2*w2*(A2-A1)));

Niter= 20000; % Nombre d'itÃ©rations
dt = 0.001; % Intervalle de temps
tf = Niter * dt; %Temps de modÃ©lisation 
t =0:dt:tf ; %Matrice temps
delta=1e-6; %pas pour la dérivation

theta_NL0=[theta10 , theta10p ,theta20, theta20p];
z=zeros(1,6);
z(1,[1:4])=theta_NL0;
z(1,[5:6])=[0 0];

x=zeros(Niter,1);

%% Solution analytique
aTheta=zeros(Niter+1,2);
aTheta(:,1)=C1*cos(w1*t+phi1)+C2*cos(w2*t+phi2);
aTheta(:,2)=C1*A1*cos(w1*t+phi1)+C2*A2*cos(w2*t+phi2);

%% solution par résolution numérique (Newton Raphson) du systéme non linéaire 

for j= 1:Niter+1  
%   initialisation du test de convergence et calcule de la valeur courante
%   de theta
    ecart=1;
    %
       while ecart > 1E-8
%      définition du systéme à résoudre
       f1=(m1+m2)*l1*z(3)+m2*l2*cos(z(1)-z(4))*z(6)+m2*l2*sin(z(1)-z(4))*z(5)^2+(m1+m2)*g*sin(z(1));
       f2=m2*l2*z(6)+m2*l1*cos(z(1)-z(4))*z(3)-m2*l1*sin(z(1)-z(4))*z(2)^2+m2*g*sin(z(4));
       f=[f1;f2];
%      calcul de la matrice Jacobienne (soit analytique soit
%      differentiation numérique)
       % analytique
%        J1=[cos(alpha) 0];
%        J2=[sin(theta)*sin(alpha) x*cos(theta)*sin(alpha)-h*cos(psi)*cos(theta)*cos(alpha)+h*sin(psi)*sin(theta)];
%        J=[J1;J2];
       % differentiation numérique
       J=Jacobienne_double_pendule(z,delta,f1,f2);
%      calcul de l'increment et de la nouvelle solution approximée
       dz=-J\f;
       z=z+dz;
%      calcul de l'ecart par rapport à la solution précédente
       ecart=abs(dz);
%      affectation de la nouvelle solution 
       x(j,1)=z(1); x(j,3)=z(4);
       end
end
%% Calcul erreur 
Erreur=zeros(Niter+1,2);
Erreur(:,1)=abs((x(:,1)-aTheta(:,1)));
Erreur(:,2)=abs((x(:,3)-aTheta(:,2)));

figure(1)
plot(t,Erreur(:,1),'.b');
hold on
plot(t,Erreur(:,2),'-r');
legend('Erreur sur theta1','Erreur sur theta2');
hold off

P1=zeros(Niter+1,2);
P1(:,1)=l1.*sin(x(:,1));
P1(:,2)=l1.*cos(x(:,1));

P2=zeros(Niter+1,1);
P2(:,1)=l2.*sin(x(:,3));
P2(:,2)=l2.*cos(x(:,3));

figure(2);
axis([-(l1+l2) (l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]); %// freeze axes
title('Double pendule')
pendule_masse1=plot(P1(1,1),-P1(1,2),'k.','MarkerSize',40,'Color','red');
hold on
pendule_tige1=plot([0,P1(1,1)],[0,-P1(1,2)],'LineWidth',1);
hold on
pendule_masse2=plot(P2(1,1),-P2(1,2),'k.','MarkerSize',40,'Color','red');
hold on
pendule_tige2=plot([P1(1,1),P2(1,1)],[-P1(1,2),-P2(1,2)],'LineWidth',1);
hold on

longueur1=sqrt(P1(:,1).^2+P1(:,2).^2);
longueur2=sqrt(P2(:,1).^2+P2(:,2).^2);

P2(:,1)=l2.*sin(x(:,3))+P1(:,1);
P2(:,2)=l2.*cos(x(:,3))+P1(:,2);

pendule_traj=plot(P2(1,1),-P2(1,2),'.b','Markersize',5);
hold on
axis([-(l1+l2) (l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]); %// freeze axes
for j = 1:Niter
    hold on
    set(pendule_masse1,'XData',P1(j,1),'YData',-P1(j,2));
    set(pendule_tige1,'XData',[0,P1(j,1)],'YData',[0,-P1(j,2)]);
    set(pendule_masse2,'XData',P2(j,1),'YData',-P2(j,2));
    set(pendule_tige2,'XData',[P1(j,1),P2(j,1)],'YData',[-P1(j,2),-P2(j,2)]);
    
    
    if j> 100
     set(pendule_traj,'XData',P2(j-100:j,1),'YData',-P2(j-100:j,2));
    end
    if j<100
        set(pendule_traj,'XData',P2(1:j,1),'YData',-P2(1:j,2));
    end 
%     plot(P2(j,1),-P2(j,2),'.b','Markersize',5)
    drawnow
end