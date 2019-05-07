%% Linéaire
%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all
%% Déclaration des variables et initalisation des constantes 
g = 9.81;         % gravité terrestre
m1 = 2;           % masse du pendule 1
m2 = 5;           % masse du pendule 2
l1 = 3;           % longueur du pendule 1
l2 = 2;           % longueur du pendule 2
theta10 = 0.02;      % angle formé par le pendule 1 avec la verticale
theta20 = 0.02;       % angle formé par le pendule 2 avec la verticale
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

%% Energies

       %Energies cinétiques
       
       Ec1=0.5*m1*(l1^2)*(theta(:,2).^2); %Energie cinétique pendule 1
       Ec2=0.5*m2*((l1^2)*(theta(:,2).^2)+(l2^2)*(theta(:,5).^2)+(2*l1*l2).*(cos(theta(:,1)-theta(:,4)).*theta(:,2).*theta(:,5))); %Energie cinétique pendule 2
        
       %Energies potentielles
       
       Ep1=(-m1*g*l1).*cos(theta(:,1)); %Energie cinétique pendule 1
       Ep2=(-m2*g)*(l1.*cos(theta(:,1))+l2.*cos(theta(:,4))); %Energie cinétique pendule 2
       
       
       









%% Affichage graphique
tf=dt*Niter;
max1=max(Ec1); min1=min(Ec1);  %Max et min de l'énergie cinétique sur le pendule 1
max2=max(Ec2); min2=min(Ec2);  %Max et min de l'énergie cinétique sur le pendule 2
maxt=max(Ec1+Ec2); mint=min(Ec1+Ec2); %Max et min de l'énergie cinétique totale

max3=max(Ep1); min3=min(Ep1);  %Max et min de l'énergie potentielle sur le pendule 1
max4=max(Ep2); min4=min(Ep2);  %Max et min de l'énergie potentielle sur le pendule 2
max34=max(max3,max4); %Max entre Ep1 et Ep2
maxp=max(Ep1+Ep2); minp=min(Ep1+Ep2); %Max et min de l'énergie potentielle totale

maxtot=max(Ep1+Ep2+Ec1+Ec2) ; mintot=min(Ep1+Ep2+Ec1+Ec2) ; 
    %% détermination de la position initiale

    context_graph=1; % tracé de la position initial
    Graph_Pendule(context_graph,P1(1,1),P1(1,2),P2(1,1),P2(1,2),l1,l2,0,tf,Ec1(1),Ec2(1),Ep1(1),Ep2(1),maxt,mint,max34,minp,maxtot,mintot);

    %% actualisation position
    for j = 1:Niter
        t=dt*j;
        context_graph=2; % reactualisation du tracé pour afficher la position courante
        Graph_Pendule(context_graph, P1(j,1),P1(j,2),P2(j,1),P2(j,2),l1,l2,t,tf,Ec1(j),Ec2(j),Ep1(j),Ep2(j),maxt,mint,max34,minp,maxtot,mintot);
        drawnow;
    end
