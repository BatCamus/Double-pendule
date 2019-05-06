%% Linéaire
%%%%%%%%%%%%%%%%%%%%%%


clear all
%% Déclaration des variables et initalisation des constantes 
g = 9.81;         % gravité terrestre
m1 = 2;           % masse du pendule 1
m2 = 5;           % masse du pendule 2
l1 = 3;           % longueur du pendule 1
l2 = 2;           % longueur du pendule 2
theta10 = 0.01;      % angle formé par le pendule 1 avec la verticale
theta20 = 0.01;       % angle formé par le pendule 2 avec la verticale
theta10p= 0.1;         %vitesse angulaire initiale du pendule 1
theta20p= 0.2;         % vitesse angulaire initiale du pendule 1
mu = m2/m1;       % rapport des masses : utile pour simplifier l'équation

Niter= 1000; % Nombre d'itérations
dt = 0.01; % Intervalle de temps
% tf = Niter * dt; %Temps de modélisation 
% t = [0:dt:tf] ; %Matrice temps

%% Constantes simplificatrices

w1 = sqrt((g*(1+mu)*(l1+l2)+g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
w2 = sqrt((g*(1+mu)*(l1+l2)-g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
A1 = (1+mu)/mu-l1*w1^2/(mu*g);
A2 = (1+mu)/mu-l1*w2^2/(mu*g);
C1 = (theta20-A2*theta10)/(A1-A2);
C2 = (A1*theta10-theta20)/(A1-A2);
phi1 = asin((theta20p-A2*theta10p)/(C1*w1*(A2-A1)));
phi2 = asin((A1*theta10p-theta20p)/(C2*w2*(A2-A1))); 

%% Déclaration et initialisation des matrices

P1=zeros(Niter+1,2); %Matrice position du premier pendule
P2=zeros(Niter+1,2); %Matrice position du deuxième pendule
theta=zeros(Niter+1,6); %Matrice angle et vitesses de rotation


theta0 = [theta10 , theta10p ,theta20, theta20p];
theta(1,1) = theta0(1);
theta(1,2) = theta0(2);
theta(1,4) = theta0(3);
theta(1,5) = theta0(4);


P1(1,1) = l1*sin(theta10); 
P1(1,2) = l1*cos(theta10);
P2(1,1) = l2*sin(theta20)+P1(1,1); 
P2(1,2) = l2*cos(theta20)+P1(1,2);

%% Boucle Euler explicite 

for i=1:Niter
   
    theta(i,3) = ((mu*g*theta(i,4))-((1+mu)*g*theta(i,1)))/l1;     % thetapp pendule 1
    theta(i,6) = ((1+mu)*g*theta(i,1)-(1+mu)*g*theta(i,4))/l2;   % thetapp pendule 2
    
    theta(i+1,1) = theta(i,1) + dt * theta(i,2);  % theta pendule 1
    theta(i+1,4) = theta(i,4) + dt * theta(i,5);  % theta pendule 2
    theta(i+1,2) = theta(i,2) + dt * theta(i,3);  % thetap pendule 1
    theta(i+1,5) = theta(i,5) + dt * theta(i,6);  % thetap pendule 2
    
   
end
 P1(:,1)=l1*sin(theta(:,1));
 P1(:,2)=l1*cos(theta(:,1));
    
 P2(:,1)=l2*sin(theta(:,4))+P1(:,1);
 P2(:,2)=l2*cos(theta(:,4))+P1(:,2);

 
%% Boucle Verlet 

for i=1:Niter
   
    theta(i,3) = ((mu*g*theta(i,4))-((1+mu)*g*theta(i,1)))/l1;     % thetapp pendule 1
    theta(i,6) = ((1+mu)*g*theta(i,1)-(1+mu)*g*theta(i,4))/l2;   % thetapp pendule 2
    
  
    theta(i+1,1) = theta(i,1) + dt * theta(i,2)+ ((dt^2)/2) * theta(i,3) ; % theta pendule 1
    theta(i+1,4) = theta(i,4) + dt * theta(i,5)+ ((dt^2)/2) * theta(i,6); % theta pendule 2

    theta(i+1,2) = theta(i,2)+ (dt/2) * theta(i,3) + (dt/2) * (((mu*g*theta(i+1,4))-((1+mu)*g*theta(i+1,1)))/l1);  % thetapp pendule 1
    theta(i+1,5) = theta(i,5)+ (dt/2) * theta(i,6) + (dt/2) * (((1+mu)*g*theta(i+1,1)-(1+mu)*g*theta(i+1,4))/l2);  % thetapp pendule 2
    
    
end
 P1(:,1)=l1*sin(theta(:,1));
 P1(:,2)=l1*cos(theta(:,1));
    
 P2(:,1)=l2*sin(theta(:,4))+P1(:,1);
 P2(:,2)=l2*cos(theta(:,4))+P1(:,2);


%% Affichage graphique


figure(1);
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

for j = 1:Niter
    set(pendule_masse1,'XData',P1(j,1),'YData',-P1(j,2));
    set(pendule_tige1,'XData',[0,P1(j,1)],'YData',[0,-P1(j,2)]);
    set(pendule_masse2,'XData',P2(j,1),'YData',-P2(j,2));
    set(pendule_tige2,'XData',[P1(j,1),P2(j,1)],'YData',[-P1(j,2),-P2(j,2)]);
    axis([-(l1+l2) (l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]); %// freeze axes
    pause(0.005)
end
