%% pendule double non linéaire methode Newmark

clear all
close all
global m1 m2 l1 l2 g mu 
global M NUM delta Niter

%% Declaration variable
g = 9.81;         % gravitÃ© terrestre
m1 = 2;           % masse du pendule 1
m2 = 3;           % masse du pendule 2
l1 = 3;           % longueur du pendule 1                                                            
l2 =2 ;% longueur du pendule 2
theta10=30*pi/180;      % angle forme par le pendule 1 avec la verticale
theta20= 5*pi/180;        % angle forme par le pendule 2 avec la verticale
theta10p= 0;         %vitesse angulaire initiale du pendule 1
theta20p= 1;         % vitesse angulaire initiale du pendule 1
theta10pp = 0;     % accÃ©lÃ©ration angulaire initiale du pendule 1
theta20pp = 0;     % accÃ©lÃ©ration angulaire initiale du pendule 2
mu = m2/m1;       % rapport des masses : utile pour simplifier l'Ã©quation

scz=get(0,'screensize'); %Taille écran


w1 = sqrt((g*(1+mu)*(l1+l2)+g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
w2 = sqrt((g*(1+mu)*(l1+l2)-g*sqrt((1+mu)^2*(l1+l2)^2-4*(1+mu)*l1*l2))/(2*l1*l2));
A1 = (1+mu)/mu-l1*w1^2/(mu*g);
A2 = (1+mu)/mu-l1*w2^2/(mu*g);
C1 = (theta20-A2*theta10)/(A1-A2);
C2 = (A1*theta10-theta20)/(A1-A2);
phi1 = asin((theta20p-A2*theta10p)/(C1*w1*(A2-A1)));
phi2 = asin((A1*theta10p-theta20p)/(C2*w2*(A2-A1)));

Niter= 10000; % Nombre d'itÃ©rations
dt = 0.01; % Intervalle de temps
tf = Niter * dt; %Temps de modÃ©lisation 
t0=0;
t =t0:dt:tf ; %Matrice temps
delta=0.1; %Pas d'intégration de Fnl

M=eye(2); %Constante pour Newmark

NUM=0; % Choix de jacobienne numérique ou analytique 0 pour analytique 1 pour numérique
ERR_petit_angle=0; % Affichage erreur petit angle
ERR_ODE_45=0; % Comparaison ODE 45
ANIM=0; %Animation
POINCARE=1; %Graphe poincare
Ener_Newmark=0; %Graphe Energie Newmark
Ener_ODE45=0; %Graphe Energie ODE45


%% Solution analytique (VERIFICATION PETITS ANGLE)
aTheta=zeros(Niter+1,2);
aTheta(:,1)=C1*cos(w1*t+phi1)+C2*cos(w2*t+phi2);
aTheta(:,2)=C1*A1*cos(w1*t+phi1)+C2*A2*cos(w2*t+phi2);

%% solution par résolution numérique Newmark
X0=[theta10 ; theta10p];
dX0=[theta20; theta20p];

tic
[tt,xt,dxt]=newmark_Double_Pendule(X0,dX0,t0,dt,tf);
toc
xt=xt';
dxt=dxt';
t=t';
tt=tt';

%% Calcul erreur petits angles
if ERR_petit_angle
    Erreur=zeros(Niter+1,2);
    Erreur(:,1)=abs((xt(:,1)-aTheta(:,1)));
    Erreur(:,2)=abs((xt(:,2)-aTheta(:,2)));

    figure(1)
    plot(t,Erreur(:,1),'.b');
    hold on
    plot(t,Erreur(:,2),'-r');
    legend('Erreur sur theta1','Erreur sur theta2');
    hold off
end 

%% Comparaison ODE 45 Newamark
if ERR_ODE_45

    theta_NL0=[theta10 , theta10p ,theta20, theta20p];
    options = odeset('AbsTol',1e-11,'RelTol',1e-11); 
    [tt,x]=ode45(@Pendule_Double_Non_Lin, t ,theta_NL0,options);

    ERR_ODE_45=zeros(Niter+1,2);
    ERR_ODE_45(:,1)=abs(xt(:,1)-x(:,1));
    ERR_ODE_45(:,2)=abs(xt(:,2)-x(:,3));

    figure(2)
    set(figure(2),'position',[3*scz(3)/4-10 0 scz(3)/4 scz(4)/2.2-40]);
    plot(tt,xt(:,1))
    hold on
    plot(tt,x(:,1),'r-')
    title('Theta 1 pour ODE 45 et Newmark')
    hold off

    figure(3)
    set(figure(3),'position',[3*scz(3)/4-10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
    plot(tt,xt(:,2))
    hold on
    plot(tt,x(:,3))
    title('Theta 2 pour ODE 45 et Newmark')
    hold off

    figure(4)
    set(figure(4),'position',[100 300 scz(3)/2 scz(4)/2-40]);
    subplot(1,2,1)
    plot(tt,ERR_ODE_45(:,1))
    title('Erreur ODE 45-Newmark theta 1')


    subplot(1,2,2)
    plot(tt,ERR_ODE_45(:,2))
    title('Erreur ODE 45-Newmark theta 2')

end

%% Affichage

P1=zeros(Niter+1,2);
P1(:,1)=l1.*sin(xt(:,1));
P1(:,2)=l1.*cos(xt(:,1));

P2=zeros(Niter+1,1);
P2(:,1)=l2.*sin(xt(:,2));
P2(:,2)=l2.*cos(xt(:,2));

if ANIM
    figure(5);
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

    P2(:,1)=l2.*sin(xt(:,2))+P1(:,1);
    P2(:,2)=l2.*cos(xt(:,2))+P1(:,2);

    pendule_traj=plot(P2(1,1),-P2(1,2),'.b','Markersize',5);
    hold on
    axis([-(l1+l2) (l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]); %// freeze axes


    for j = 1:Niter
         hold on
         set(pendule_masse1,'XData',P1(j,1),'YData',-P1(j,2));
         set(pendule_tige1,'XData',[0,P1(j,1)],'YData',[0,-P1(j,2)]);
         set(pendule_masse2,'XData',P2(j,1),'YData',-P2(j,2));
         set(pendule_tige2,'XData',[P1(j,1),P2(j,1)],'YData',[-P1(j,2),-P2(j,2)]);

% 
%             if j> 100
%              set(pendule_traj,'XData',P2(j-100:j,1),'YData',-P2(j-100:j,2));
%             end
%             if j<100
%                 set(pendule_traj,'XData',P2(1:j,1),'YData',-P2(1:j,2));
%             end 
            plot(P2(j,1),-P2(j,2),'.b','Markersize',5)
            drawnow
        end
end 

%% Section de poincaré

if POINCARE

    n1=max(size(xt(:,1)),size(xt(:,2)));
    %set the index of poincare points to 1
    np1=1;
    np2=1; 
    %Création matrices ps 
    ps1=zeros(n1(1),2);
    ps2=zeros(n1(1),2);
    
   %Précision des points d'intersections
   Atol=10^(-1);

    

%



    for i=2:n1(1)
            %Trouver les points passant par un plan défini (ici theta1=0)
            if (xt(i,1)*xt(i-1,1)<0 && dxt(i,1)>0 &&  abs(xt(i,1))<2)
                
                %choix du point le plus proche entre celui à gauche et à droite du plan theta1=0
                if(abs(xt(i,1))<abs(xt(i-1,1)))
                    
                    % Sauvegarde des points d'intersection en theta1=0
                    ps1(np1,1)=xt(i,1);
                    ps1(np1,2)=dxt(i,1);
                else 
                    % Sauvegarde des points d'intersection en theta1=0
                    ps1(np1,1)=xt(i-1,1);
                    ps1(np1,2)=dxt(i-1,1);
                end 
                % Incrémentation 
                np1=np1+1;
            end
            %Trouver les points passant par un plan défini (ici theta2=0)
            if (xt(i,2)*xt(i-1,2)<0 && dxt(i,2)>0 && abs(xt(i,2))<2) 
                
                 %choix du point le plus proche entre celui à gauche et à droite du plan theta2=0
                 if(abs(xt(i,2))<abs(xt(i-1,2)))
                     
                    % Sauvegarde des points d'intersection en theta2=0
                    ps2(np2,1)=xt(i,2);
                    ps2(np2,2)=dxt(i,2);
                 else 
                    % Sauvegarde des points d'intersection en theta2=0
                    ps2(np2,1)=xt(i-1,2);
                    ps2(np2,2)=dxt(i-1,2);
                 end  
                %  Incrémentation 
                np2=np2+1;
            end
    end


    figure(6) 
    set(figure(6),'position',[10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
    plot(xt(:,1),dxt(:,1),'c-','Markersize',2)
    xlabel('theta1')
    ylabel('d(theta1)/dt (rad/s)') 
    title('Portrait de phase en theta1') 
    hold on 
    axis([min(xt(:,1))-0.2 max(xt(:,1))+0.2 min(dxt(:,1))-0.2 max(dxt(:,1))+0.2]);

    %Boucle affichage de la section de poincaré
    for i=1:np1-1
        plot(ps1(i,1),ps1(i,2),'r+','markersize', 5)
        % Possibilité de faire un affichage en temps réel
        %pause(2);
    end



                    %%%%%%%%%%%%

    figure(7)
    set(figure(7),'position',[10 0 scz(3)/4 scz(4)/2.2-40]);
    x1=zeros(np1-1,2);
    for i=1:np1-1
        x1(i,1)=ps1(i,1);
        x1(i,2)=ps1(i,2);
    end 
    plot(x1(:,1),x1(:,2),'r+','markersize', 5)
    xlabel('theta1')
    ylabel('d(theta1)/dt (rad/s)')
    title('Section de poincaré en theta1')
    axis([min(xt(:,1))-0.2 max(xt(:,1))+0.2 min(dxt(:,1))-0.2 max(dxt(:,1))+0.2]);

                     %%%%%%%%%%%%

    figure(8) 
    set(figure(8),'position',[3*scz(3)/4-10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
    plot(xt(:,2),dxt(:,2),'c-','Markersize',2)
    xlabel('theta2')
    ylabel('d(theta2)/dt (rad/s)')
    title('Portrait de phase en theta2')
    axis([min(xt(:,2))-0.2 max(xt(:,2))+0.2 min(dxt(:,2))-0.2 max(dxt(:,2))+0.2]);
    hold on     

    %Boucle affichage de la section de poincaré
    for i=1:np2-1
        plot(ps2(i,1),ps2(i,2),'r+','markersize', 5)
        % Possibilité de faire un affichage en temps réel
        %pause(2);
    end

                     %%%%%%%%%%%%


    figure(9)
    set(figure(9),'position',[3*scz(3)/4-10 0 scz(3)/4 scz(4)/2.2-40]);
    x2=zeros(np2-1,2);

    for i=1:np2-1
        x2(i,1)=ps2(i,1);
        x2(i,2)=ps2(i,2);
    end 
    
    plot(x2(:,1),x2(:,2),'r+','markersize', 5)
    xlabel('theta2')
    ylabel('d(theta2)/dt (rad/s)')
    title('Section de poincaré en theta2')
    axis([min(xt(:,2))-0.2 max(xt(:,2))+0.2 min(dxt(:,2))-0.2 max(dxt(:,2))+0.2]);
end 

%% Energies ODE 45
if Ener_ODE45
           %Energies cinétiques

           Ec1=0.5*m1*(l1^2)*(x(:,2).^2); %Energie cinétique pendule 1
           Ec2=0.5*m2*((l1^2)*(x(:,2).^2)+(l2^2)*(x(:,4).^2)+(2*l1*l2).*(cos(x(:,1)-x(:,3)).*x(:,2).*x(:,4))); %Energie cinétique pendule 2

           %Energies potentielles

           Ep1=(-m1*g*l1).*cos(x(:,1)); %Energie cinétique pendule 1
           Ep2=(-m2*g)*(l1.*cos(x(:,1))+l2.*cos(x(:,3))); %Energie cinétique pendule 2

    % Affichage graphique
    max1=max(Ec1); min1=min(Ec1);  %Max et min de l'énergie cinétique sur le pendule 1
    max2=max(Ec2); min2=min(Ec2);  %Max et min de l'énergie cinétique sur le pendule 2
    maxt=max(Ec1+Ec2); mint=min(Ec1+Ec2); %Max et min de l'énergie cinétique totale

    max3=max(Ep1); min3=min(Ep1);  %Max et min de l'énergie potentielle sur le pendule 1
    max4=max(Ep2); min4=min(Ep2);  %Max et min de l'énergie potentielle sur le pendule 2
    max34=max(max3,max4); %Max entre Ep1 et Ep2
    maxp=max(Ep1+Ep2); minp=min(Ep1+Ep2); %Max et min de l'énergie potentielle totale

    maxtot=max(Ep1+Ep2+Ec1+Ec2) ; mintot=min(Ep1+Ep2+Ec1+Ec2) ; 
        % détermination de la position initiale

        context_graph=1; % tracé de la position initial
        Graph_Pendule(context_graph,P1(1,1),P1(1,2),P2(1,1),P2(1,2),l1,l2,0,tf,Ec1(1),Ec2(1),Ep1(1),Ep2(1),maxt,mint,max34,minp,maxtot,mintot);

        % actualisation position
        for j = 1:Npas
            t=dt*j;
            context_graph=2; % reactualisation du tracé pour afficher la position courante
            Graph_Pendule(context_graph, P1(j,1),P1(j,2),P2(j,1),P2(j,2),l1,l2,t,tf,Ec1(j),Ec2(j),Ep1(j),Ep2(j),maxt,mint,max34,minp,maxtot,mintot);
            drawnow;
        end
end

%% Energies Newmark
if Ener_Newmark
           %Energies cinétiques

           Ec1=0.5*m1*(l1^2)*(dxt(:,1).^2); %Energie cinétique pendule 1
           Ec2=0.5*m2*((l1^2)*(dxt(:,1).^2)+(l2^2)*(dxt(:,2).^2)+(2*l1*l2).*(cos(xt(:,1)-xt(:,2)).*dxt(:,1).*dxt(:,2))); %Energie cinétique pendule 2

           %Energies potentielles

           Ep1=(-m1*g*l1).*cos(xt(:,1)); %Energie cinétique pendule 1
           Ep2=(-m2*g)*(l1.*cos(xt(:,1))+l2.*cos(xt(:,2))); %Energie cinétique pendule 2

    % Affichage graphique
    max1=max(Ec1); min1=min(Ec1);  %Max et min de l'énergie cinétique sur le pendule 1
    max2=max(Ec2); min2=min(Ec2);  %Max et min de l'énergie cinétique sur le pendule 2
    maxt=max(Ec1+Ec2); mint=min(Ec1+Ec2); %Max et min de l'énergie cinétique totale

    max3=max(Ep1); min3=min(Ep1);  %Max et min de l'énergie potentielle sur le pendule 1
    max4=max(Ep2); min4=min(Ep2);  %Max et min de l'énergie potentielle sur le pendule 2
    max34=max(max3,max4); %Max entre Ep1 et Ep2
    maxp=max(Ep1+Ep2); minp=min(Ep1+Ep2); %Max et min de l'énergie potentielle totale

    maxtot=max(Ep1+Ep2+Ec1+Ec2) ; mintot=min(Ep1+Ep2+Ec1+Ec2) ; 
        % détermination de la position initiale

        context_graph=1; % tracé de la position initial
        Graph_pendule(context_graph,P1,P2,l1,l2,0,tf,Ec1(1),Ec2(1),Ep1(1),Ep2(1),maxt,mint,max34,minp,maxtot,mintot,1);

        %%actualisation position

        for j = 2:Niter
            t=dt*j;
            context_graph=2; % reactualisation du tracé pour afficher la position courante
            Graph_pendule(context_graph, P1,P2,l1,l2,t,tf,Ec1(j),Ec2(j),Ep1(j),Ep2(j),maxt,mint,max34,minp,maxtot,mintot,j);
            drawnow;
        end
end

%% Grille erreur relative : Schéma Newmark / analytique
if grilleErr

    dtheta=0.5;                         % Pas d'angle           
    Range=10;                           % Angles extremes à atteindre
    
    t1d=-Range:dtheta:Range;
    t2d=-Range:dtheta:Range;
    t1=t1d*pi/180;
    t2=t2d*pi/180;

    W1=zeros(length(t1),length(t2));
    W2=zeros(length(t1),length(t2));

    for j=1:length(t1)                  % Boucle sur Theta1
        for i=1:length(t2)              % Boucle sur Theta2
            %ode45
            
            X0=[t1(1,j); t2(1,i)];
            dX0=[0; 0];
            [tt,xt,dxt]=newmark_Double_Pendule(X0,dX0,t0,dt,tf);
            xt=xt';
            
            %Analytique
            
            C1 = (t2(i)-A2*t1(j))/(A1-A2);
            C2 = (A1*t1(j)-(t2(i)))/(A1-A2);
            phi1 = 0;                   % Vitesse initiale nulle!
            phi2 = 0;                   % Vitesse initiale nulle!

            aTheta=zeros(Npas+1,2);
            aTheta(:,1)=C1*cos(w1*t+phi1)+C2*cos(w2*t+phi2);
            aTheta(:,2)=C1*A1*cos(w1*t+phi1)+C2*A2*cos(w2*t+phi2);

            %Erreur
            Erreur=zeros(Npas+1,2);
            Erreur(:,1)=abs((xt(:,1)-aTheta(:,1)));
            Erreur(:,2)=abs((xt(:,2)-aTheta(:,2)));

            ErreurRel1=abs(Erreur(:,1)/max(aTheta(:,1)));
            ErreurRel2=abs(Erreur(:,2)/max(aTheta(:,2)));

            erreur1=mean(ErreurRel1)*100;
            erreur2=mean(ErreurRel2)*100;
            erreur=max(erreur1,erreur2);

            W1(j,i)=erreur1;
            W2(j,i)=erreur2;


        end
            figure(25)
            pcolor(t1d,t2d,W1)
            caxis([0 5])
%             caxis('auto')
            xlabel('Theta1 en degres');
            ylabel('Theta2 en degres');
            title('Erreur relative sur theta1 en %');
            colorbar('EastOutside')
            drawnow
           
            figure(26)
            pcolor(t1d,t2d,W2)
            caxis([0 5])
%             caxis('auto')
            xlabel('Theta1 en degres');
            ylabel('Theta2 en degres');
            title('Erreur relative sur theta2 en %');
            colorbar('EastOutside')
            drawnow
    end
end