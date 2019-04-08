%% Solution analytique
%%%%%%%%%%%%%%%%%%%%%%

%% Constantes

g = 9.81;

m1=2; %masse pendule 1
m2=3; % masse pendule 2
mu = m1/m2;

l1 = 3; %longueur pendule 1
l2 = 2; %longueur pendule 2

Niter= ?    ; % Nombre d'itérations
tau= ?      ; % Intervalle de temps

%% Constantes simplificatrices

w1 = sqrt((g*(1+mu)*(l1+l2)+g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
w2 = sqrt((g*(1+mu)*(l1+l2)-g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
A1 = (1+mu)/mu-l1*w1^2/(mu*g);
A2 = (1+mu)/mu-l1*w2^2/(mu*g);
C1 = (theta20-A2*theta10)/(A1-A2);
C2 = (A1*theta10-theta20)/(A1-A2);
phi1 = asin((theta20p-A2*theta10p)/(C1*w1*(A2-A1));
phi2 = asin((A1*theta10p-theta20p)/(C2*w2*(A2-A1)); 

%% Déclaration et initialisation des matrices

X=zeros(Niter+1,1); %Matrice position n°1
Y=zeros(Niter+1,1); %Matrice position n°2
theta=zeros(Niter+1,4); %Matrice angle et vitesses de rotation

theta10 = ;
theta10p = ;
theta20 = ;
theta20p = ;
theta0 = [theta10 ,theta10p,theta20,theta20p];
theta(1,:) = theta0;
theta1=theta(:,1) % theta1
theta1p=theta(:,2) % theta point 1
theta2=theta(:,3) % theta 2
theta2p=theta(:,4) % theta point 2 


