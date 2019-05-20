%% Déclaration des variables et initalisation des constantes 
clear all
close all
global m1 m2 l1 l2 g mu
g = 9.81;         % gravité terrestre
m1 = 1;           % masse du pendule 1
m2 = 0.9;           % masse du pendule 2
l1 = 1;           % longueur du pendule 1
l2 = 0.9;           % longueur du pendule 2
theta10 =0.5;      % angle formé par le pendule 1 avec la verticale
theta20 = 0.1;        % angle formé par le pendule 2 avec la verticale
theta10p= 0;         %vitesse angulaire initiale du pendule 1
theta20p= 0;         % vitesse angulaire initiale du pendule 1
theta10pp = 0;     % accélération angulaire initiale du pendule 1
theta20pp = 0;     % accélération angulaire initiale du pendule 2
mu = m2/m1;       % rapport des masses : utile pour simplifier l'équation


w1 = sqrt((g*(1+mu)*(l1+l2)+g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
w2 = sqrt((g*(1+mu)*(l1+l2)-g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
A1 = (1+mu)/mu-l1*w1^2/(mu*g);
A2 = (1+mu)/mu-l1*w2^2/(mu*g);
C1 = (theta20-A2*theta10)/(A1-A2);
C2 = (A1*theta10-theta20)/(A1-A2);
phi1 = asin((theta20p-A2*theta10p)/(C1*w1*(A2-A1)));
phi2 = asin((A1*theta10p-theta20p)/(C2*w2*(A2-A1)));

Niter= 2000; % Nombre d'itérations
dt = 0.01; % Intervalle de temps
tf = Niter * dt; %Temps de modélisation 
t =0:dt:tf ; %Matrice temps

theta_NL0=[theta10 , theta10p ,theta20, theta20p];
    
[tt,x]=ode45(@Pendule_Double_Non_Lin, t ,theta_NL0);

aTheta=zeros(Niter+1,2);
aTheta(:,1)=C1*cos(w1*t+phi1)+C2*cos(w2*t+phi2);
aTheta(:,2)=C1*A1*cos(w1*t+phi1)+C2*A2*cos(w2*t+phi2);

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


