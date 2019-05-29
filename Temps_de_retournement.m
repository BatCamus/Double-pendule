%% Déclaration des variables et initalisation des constantes 
clear all
%close all
global m1 m2 l1 l2 g mu
g = 9.81;           % gravité terrestre
m1 = 2;             % masse du pendule 1
m2 = 3;           % masse du pendule 2
l1 = 3;             % longueur du pendule 1
l2 = 2;           % longueur du pendule 2
theta10d =80;       % angle formé par le pendule 1 avec la verticale
theta20d =80;      % angle formé par le pendule 2 avec la verticale
theta10=theta10d*pi/180;
theta20=theta20d*pi/180;
theta10p= 0;        %vitesse angulaire initiale du pendule 1
theta20p= 0;        % vitesse angulaire initiale du pendule 1
theta10pp = 0;      % accélération angulaire initiale du pendule 1
theta20pp = 0;      % accélération angulaire initiale du pendule 2
mu = m2/m1;         % rapport des masses : utile pour simplifier l'équation

anim=0;
grille=1;



w1 = sqrt((g*(1+mu)*(l1+l2)+g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
w2 = sqrt((g*(1+mu)*(l1+l2)-g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
A1 = (1+mu)/mu-l1*w1^2/(mu*g);
A2 = (1+mu)/mu-l1*w2^2/(mu*g);
C1 = (theta20-A2*theta10)/(A1-A2);
C2 = (A1*theta10-theta20)/(A1-A2);
phi1 = asin((theta20p-A2*theta10p)/(C1*w1*(A2-A1)));
phi2 = asin((A1*theta10p-theta20p)/(C2*w2*(A2-A1)));



Niter= 10000;      % Nombre d'itérations
dt = 0.005;          % Intervalle de temps
tf = Niter * dt;    %Temps de modélisation 
t =0:dt:tf ;        %Matrice temps


delta=0.1; %Pas d'intégration de Fnl

M=eye(2); %Constante pour Newmark

NUM=0; % Choix de jacobienne numérique ou analytique 0 pour analytique 1 pour numérique
%% Solution numérique
 theta_NL0=[theta10 , theta10p ,theta20, theta20p];
 [tt,x]=ode45(@Pendule_Double_Non_Lin, t ,theta_NL0);

%% Calul position en cartésien
 
P1=zeros(Niter+1,2);
P1(:,1)=l1.*sin(x(:,1));
P1(:,2)=l1.*cos(x(:,1));

P2=zeros(Niter+1,1);
P2(:,1)=l2.*sin(x(:,3))+P1(:,1);
P2(:,2)=l2.*cos(x(:,3))+P1(:,2);

% Solution analytique
 aTheta=zeros(Niter+1,2);
 aTheta(:,1)=C1*cos(w1*t+phi1)+C2*cos(w2*t+phi2);
 aTheta(:,2)=C1*A1*cos(w1*t+phi1)+C2*A2*cos(w2*t+phi2);
 
 
%  %% Energies
% 
%        %Energies cinétiques
%        
%        Ec1=0.5*m1*(l1^2)*(x(:,2).^2); %Energie cinétique pendule 1
%        Ec2=0.5*m2*((l1^2)*(x(:,2).^2)+(l2^2)*(x(:,4).^2)+(2*l1*l2).*(cos(x(:,1)-x(:,3)).*x(:,2).*x(:,4))); %Energie cinétique pendule 2
%         
%        %Energies potentielles
%        
%        Ep1=(-m1*g*l1).*cos(x(:,1)); %Energie cinétique pendule 1
%        Ep2=(-m2*g)*(l1.*cos(x(:,1))+l2.*cos(x(:,3))); %Energie cinétique pendule 2
%        
% % Affichage graphique
% tf=dt*Niter;
% max1=max(Ec1); min1=min(Ec1);  %Max et min de l'énergie cinétique sur le pendule 1
% max2=max(Ec2); min2=min(Ec2);  %Max et min de l'énergie cinétique sur le pendule 2
% maxt=max(Ec1+Ec2); mint=min(Ec1+Ec2); %Max et min de l'énergie cinétique totale
% 
% max3=max(Ep1); min3=min(Ep1);  %Max et min de l'énergie potentielle sur le pendule 1
% max4=max(Ep2); min4=min(Ep2);  %Max et min de l'énergie potentielle sur le pendule 2
% max34=max(max3,max4); %Max entre Ep1 et Ep2
% maxp=max(Ep1+Ep2); minp=min(Ep1+Ep2); %Max et min de l'énergie potentielle totale
% 
% maxtot=max(Ep1+Ep2+Ec1+Ec2) ; mintot=min(Ep1+Ep2+Ec1+Ec2) ; 
%     % détermination de la position initiale
% 
%     context_graph=1; % tracé de la position initial
%     Graph_Pendule(context_graph,P1(1,1),P1(1,2),P2(1,1),P2(1,2),l1,l2,0,tf,Ec1(1),Ec2(1),Ep1(1),Ep2(1),maxt,mint,max34,minp,maxtot,mintot);
% 
%     % actualisation position
%     for j = 1:Niter
%         t=dt*j;
%         context_graph=2; % reactualisation du tracé pour afficher la position courante
%         Graph_Pendule(context_graph, P1(j,1),P1(j,2),P2(j,1),P2(j,2),l1,l2,t,tf,Ec1(j),Ec2(j),Ep1(j),Ep2(j),maxt,mint,max34,minp,maxtot,mintot);
%         drawnow;
%     end

%% solution par résolution numérique Newmark
% X0=[theta10 ; theta10p];
% dX0=[theta20; theta20p];
% 
% tic
% [tt,xt,dxt]=newmark_Double_Pendule(X0,dX0,t0,dt,tf);
% toc
% xt=xt';

%% Grille  
if grille

    dtheta=1;
    Range=180;
    
    t1d=-Range:dtheta:Range;
    t2d=-Range:dtheta:Range;
    t1=t1d*pi/180;
    t2=t2d*pi/180;

    W1=zeros(length(t1),length(t2));
    W2=zeros(length(t1),length(t2));

    for j=1:length(t1)
        for i=1:length(t2)
            %ode45
            tic
            theta_NL0=[t1(j) , 0 ,t2(i),0];
            [tt,x]=ode45(@Pendule_Double_Non_Lin, t ,theta_NL0);

            %Analytique
%             C1 = (t2(i)-A2*t1(j))/(A1-A2);
%             C2 = (A1*t1(j)-(t2(i)))/(A1-A2);
%             phi1 = 0;
%             phi2 = 0;
% 
%             aTheta=zeros(Niter+1,2);
%             aTheta(:,1)=C1*cos(w1*t+phi1)+C2*cos(w2*t+phi2);
%             aTheta(:,2)=C1*A1*cos(w1*t+phi1)+C2*A2*cos(w2*t+phi2);

            %Erreur
%             Erreur=zeros(Niter+1,2);
%             Erreur(:,1)=abs((x(:,1)-aTheta(:,1)));
%             Erreur(:,2)=abs((x(:,3)-aTheta(:,2)));
% 
%             ErreurRel1=abs(Erreur(:,1)/max(aTheta(:,1)));
%             ErreurRel2=abs(Erreur(:,2)/max(aTheta(:,2)));
% 
%             erreur1=mean(ErreurRel1)*100;
%             erreur2=mean(ErreurRel2)*100;
%             erreur=max(erreur1,erreur2);

            P1(:,2)=l1.*cos(x(:,1));
            P2(:,2)=l2.*cos(x(:,3))+P1(:,2);
            k=1;
            z=false;
            while z==false && k<Niter
                if P2(k,2)<0
                    z=true;
                    Tr=k*dt;
                end
                k=k+1;
            end
            if z==true
                W1(j,i)=Tr; 
            else
                W1(j,i)=NaN;
            end
%             W2(j,i)=erreur2;

            toc
        end
%             W1(0,0)=0;
%             W2(0,0)=0;
            figure(25)
            h=pcolor(t1d,t2d,log10(W1))
            set(h, 'EdgeColor', 'none');
%             caxis([0 5])
            caxis('auto')
            xlabel('Theta1 initial en degres');
            ylabel('Theta2 initial en degres');
            title('Temps de retournement');
            colorbar('EastOutside')
            drawnow
           
%             figure(26)
%             pcolor(t1d,t2d,W2)
%             caxis([0 5])
%             %caxis('auto')
%             xlabel('Theta1 en degres');
%             ylabel('Theta2 en degres');
%             title('Erreur relative sur theta2');
%             colorbar('EastOutside')
%             drawnow
    end
    %% Position
    figure(2);
    plot(t,x(:,1),'Color','red');
    hold on
    plot(t,aTheta(:,1),'Color','blue')
    hold off

    figure(3)
    plot(t,x(:,3),'Color','red');
    hold on
    plot(t,aTheta(:,2),'Color','blue')
end
%% Erreur
Erreur=zeros(Niter+1,2);

Erreur(:,1)=abs((x(:,1)-aTheta(:,1)));
Erreur(:,2)=abs((x(:,3)-aTheta(:,2)));

ErreurRel1=abs(Erreur(:,1)/max(aTheta(:,1)));
ErreurRel2=abs(Erreur(:,2)/max(aTheta(:,2)));

erreur1=mean(ErreurRel1);
erreur2=mean(ErreurRel2);

erreur=max(erreur1,erreur2);

%% Animation
if anim 
    figure(5);
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

    axis([-1.2*(l1+l2) 1.2*(l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]); % freeze axes

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
end


