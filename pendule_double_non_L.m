clear all
close all
clc


global m1 m2 l1 l2 g mu
g = 9.81;         % gravité terrestre
m1 = 2;           % masse du pendule 1
m2 = 3;           % masse du pendule 2
l1 = 3;           % longueur du pendule 1
l2 = 2;           % longueur du pendule 2

scz=get(0,'screensize'); %Taille écran
theta10 =pi/6;      % angle formé par le pendule 1 avec la verticale
theta20 = pi/6;        % angle formé par le pendule 2 avec la verticale
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

Niter= 50000; % Nombre d'itérations
dt = 0.01; % Intervalle de temps
tf = Niter * dt; %Temps de modélisation 
t =0:dt:tf ; %Matrice temps
T1=2*pi/w1; %Periode du premier pendule
T2=2*pi/w2; %Periode du second pendule

theta_NL0=[theta10 , theta10p ,theta20, theta20p];

options= odeset('RelTol',1e-10,'AbsTol',1e-12);
[tt,x]=ode45(@Pendule_Double_Non_Lin, t ,theta_NL0, options);

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
xlabel('Z(m)')
ylabel('X(m)')
title('LE pendule double');

axis([-(l1+l2) (l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]); %// freeze axes
% for j = 1:Niter
%     hold on
%     set(pendule_masse1,'XData',P1(j,1),'YData',-P1(j,2));
%     set(pendule_tige1,'XData',[0,P1(j,1)],'YData',[0,-P1(j,2)]);
%     set(pendule_masse2,'XData',P2(j,1),'YData',-P2(j,2));
%     set(pendule_tige2,'XData',[P1(j,1),P2(j,1)],'YData',[-P1(j,2),-P2(j,2)]);
%     
%     
%     if j> 100
%      set(pendule_traj,'XData',P2(j-100:j,1),'YData',-P2(j-100:j,2));
%     end
%     if j<100
%         set(pendule_traj,'XData',P2(1:j,1),'YData',-P2(1:j,2));
%     end 
%     drawnow
%     grid on 
% end






%% Section de poincaré


    n1=size(x(:,1));
    %set the index of poincare points to 1
    np1=1;
    np2=1; 
    %Création matrices ps 
    ps1=zeros(n1(1),2);
    ps2=zeros(n1(1),2);
    
    %2*x(i,1))>x(i,2)-10^(-4) && 2*x(i,1)<x(i,2)+10^(-4)





    for i=1:n1(1)
            % detect the cros-section of the trajectory with the plane y1-y2
            if((abs(x(i,1))<10^(-3))) %||(abs(x(i,3))<10^(-3))
                % store detected cross-section point y1,y2 to ps1,ps2
                ps1(np1,1)=x(i,1);
                ps1(np1,2)=x(i,2);
                % increase the index of poincare point
                np1=np1+1;
            end
             if((abs(x(i,3))<10^(-3)))
                % store detected cross-section point y1,y2 to ps1,ps2
                ps2(np2,1)=x(i,3);
                ps2(np2,2)=x(i,4);
                % increase the index of poincare point
                np2=np2+1;
            end
    end

     

    figure(3) 
    set(figure(3),'position',[10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
    plot(x(:,1),x(:,2),'c-','Markersize',2)
    xlabel('theta1')
    ylabel('d(theta1)/dt (rad/s)')
    title('Portrait de phase en theta1')
    hold on 
    axis([min(x(:,1))-0.2 max(x(:,1))+0.2 min(x(:,2))-0.2 max(x(:,2))+0.2]);
    
    %Boucle affichage de la section de poincaré
    for i=1:np1-1
        plot(ps1(i,1),ps1(i,2),'r+','markersize', 5)
        % Possibilité de faire un affichage en temps réel
        % pause(2);
    end

    
                    %%%%%%%%%%%%
                   
    figure(4)
    set(figure(4),'position',[10 0 scz(3)/4 scz(4)/2.2-40]);
    p1=zeros(np1-1,2);
    for i=1:np1-1
        p1(i,1)=ps1(i,1);
        p1(i,2)=ps1(i,2);
    end 
    plot(p1(:,1),p1(:,2),'r+','markersize', 5)
    xlabel('theta1')
    ylabel('d(theta1)/dt (rad/s)')
    title('Section de poincaré en theta1')
    axis([min(x(:,1))-0.2 max(x(:,1))+0.2 min(x(:,2))-0.2 max(x(:,2))+0.2]);
    
    
                     %%%%%%%%%%%%
    
    figure(5) 
    set(figure(5),'position',[3*scz(3)/4-10 scz(4)/2-20 scz(3)/4 scz(4)/2.2-40]);
    plot(x(:,3),x(:,4),'c-','Markersize',2)
    xlabel('theta2')
    ylabel('d(theta2)/dt (rad/s)')
    title('Portrait de phase en theta2')
    axis([min(x(:,3))-0.2 max(x(:,3))+0.2 min(x(:,4))-0.2 max(x(:,4))+0.2]);
    hold on     
    
    %Boucle affichage de la section de poincaré
    for i=1:np2-1
        plot(ps2(i,1),ps2(i,2),'r+','markersize', 5)
        % Possibilité de faire un affichage en temps réel
        % pause(2);
    end

                     %%%%%%%%%%%%
        
    
    figure(6)
    set(figure(6),'position',[3*scz(3)/4-10 0 scz(3)/4 scz(4)/2.2-40]);
    p2=zeros(np2-1,2);
    
    for i=1:np2-1
        p2(i,1)=ps2(i,1);
        p2(i,2)=ps2(i,2);
    end 
    plot(p2(:,1),p2(:,2),'r+','markersize', 5)
    xlabel('theta2')
    ylabel('d(theta2)/dt (rad/s)')
    title('Section de poincaré en theta2')
    axis([min(x(:,3))-0.2 max(x(:,3))+0.2 min(x(:,4))-0.2 max(x(:,4))+0.2]);
    
                    
 %% Diagramme de biffurcation 
 
 

for b=0.9:0.01:2
    for i=0:2
        clear x 
        options= odeset('RelTol',1e-10,'AbsTol',1e-12);
        [tt,x]=ode45(@Pendule_Double_Non_Lin, [0:l1*2*pi/w1:l1*2*50/w1],[0 1], options);
        
        if i==0
            figure(7)
            plot(c,x(30:50,2))
        end 
        if i==1
            figure(7)
            plot(c,x(30:50,2))
        end 
        if i==2
            figure(7)
            plot(c,x(30:50,2))
        end 
    end 
end 


