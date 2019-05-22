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
theta10 =pi/12;      % angle forme par le pendule 1 avec la verticale
theta20 = pi/6;        % angle forme par le pendule 2 avec la verticale
theta10p= 0;         %vitesse angulaire initiale du pendule 1
theta20p= 0;         % vitesse angulaire initiale du pendule 1
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
ERR_ODE_45=1; % Comparaison ODE 45
ANIM=0; %Animation
POINCARE=1; %Graphe poincare

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
if(ERR_petit_angle==1)
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
if(ERR_ODE_45==1)

    theta_NL0=[theta10 , theta10p ,theta20, theta20p];
    options = odeset('AbsTol',1e-11,'RelTol',1e-11); 
    [tt,x]=ode45(@Pendule_Double_Non_Lin, t ,theta_NL0,options);

    ERR_ODE_45=zeros(Niter+1,2);
    ERR_ODE_45(:,1)=abs(xt(:,1)-x(:,1));
    ERR_ODE_45(:,2)=abs(xt(:,2)-x(:,3));

    figure()
    plot(tt,xt(:,1))
    hold on
    plot(tt,x(:,1),'r-')
    title('Theta 1 pour ODE 45 et Newmark')
    hold off

    figure()
    plot(tt,xt(:,2))
    hold on
    plot(tt,x(:,3))
    title('Theta 2 pour ODE 45 et Newmark')
    hold off

    figure()
    subplot(1,2,1)
    plot(tt,ERR_ODE_45(:,1))
    title('Erreur ODE 45-Newmark theta 1')


    subplot(1,2,2)
    plot(tt,ERR_ODE_45(:,2))
    title('Erreur ODE 45-Newmark theta 2')

end

%% Affichage
if ANIM
    P1=zeros(Niter+1,2);
    P1(:,1)=l1.*sin(xt(:,1));
    P1(:,2)=l1.*cos(xt(:,1));

    P2=zeros(Niter+1,1);
    P2(:,1)=l2.*sin(xt(:,2));
    P2(:,2)=l2.*cos(xt(:,2));
  
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
            % detect the cros-section of the trajectory with the plane
            % y1-y2 %Penser à prendre le point le plus proche
            
            if (xt(i,1)*xt(i-1,1)<0 && dxt(i,1)>0 &&  abs(xt(i,1))<2)
                if(abs(xt(i,1))<abs(xt(i-1,1)))
                    % store detected cross-section point y1,y2 to ps1,ps2
                    ps1(np1,1)=xt(i,1);
                    ps1(np1,2)=dxt(i,1);
                else 
                    % store detected cross-section point y1,y2 to ps1,ps2
                    ps1(np1,1)=xt(i-1,1);
                    ps1(np1,2)=dxt(i-1,1);
                end 
                % increase the index of poincare point
                np1=np1+1;
            end
             if (xt(i,2)*xt(i-1,2)<0 && dxt(i,2)>0 && abs(xt(i,2))<2) 
                 if(abs(xt(i,2))<abs(xt(i-1,2)))
                    % store detected cross-section point y1,y2 to ps1,ps2
                    ps2(np2,1)=xt(i,2);
                    ps2(np2,2)=dxt(i,2);
                 else 
                    % store detected cross-section point y1,y2 to ps1,ps2
                    ps2(np2,1)=xt(i-1,2);
                    ps2(np2,2)=dxt(i-1,2);
                 end  
                % increase the index of poincare point
                np2=np2+1;
            end
    end


    figure(3) 
    set(figure(3),'position',[10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
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

    figure(4)
    set(figure(4),'position',[10 0 scz(3)/4 scz(4)/2.2-40]);
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

    figure(5) 
    set(figure(5),'position',[3*scz(3)/4-10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
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


    figure(6)
    set(figure(6),'position',[3*scz(3)/4-10 0 scz(3)/4 scz(4)/2.2-40]);
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
