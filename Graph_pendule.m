function Graph_pendule(contexte,P1,P2,l1,l2,t,tf,Ec1,Ec2,Ep1,Ep2,maxt,mint,max34,minp,maxtot,mintot,j)

global pendule_masse1 pendule_tige1 pendule_masse2 pendule_tige2 Energie_cin1 Energie_cin2 Energie_cint Energie_pott Energie_t Energie_pot1 Energie_pot2 pendule_traj


Ect=Ec1+Ec2;
Ept=Ep1+Ep2;
scz=get(0,'screensize'); %Taille ecran
    
%%
        figure(12);

        set(figure(12),'position',[20+scz(3)/3 220+4*scz(4)/12 scz(3)/3 scz(4)/3]);
        
      
        xlabel('X(dm)')
        ylabel('Z(dm)')
        title('Mouvement d un double pendule');
       

          if contexte==1 %position initiale

            pendule_masse1=plot(P1(j,1),-P1(j,2),'k.','MarkerSize',35,'Color','red');
            hold on
            pendule_tige1=plot([0,P1(j,1)],[0,-P1(j,2)],'LineWidth',2,'Color','black');
            hold on
            pendule_masse2=plot(P2(j,1),-P2(j,2),'k.','MarkerSize',35,'Color','red');
            hold on
            pendule_tige2=plot([P1(j,1),P2(j,1)],[-P1(j,2),-P2(j,2)],'LineWidth',2,'Color','blue');
            hold on
            pendule_traj=plot(P2(j,1),-P2(j,2),'.b','Markersize',5);
            hold on 
    
            axis([-6 6 -6 6]);

          elseif contexte==2 %position actualise

            set(pendule_masse1,'XData',P1(j,1),'YData',-P1(j,2));
            set(pendule_tige1,'XData',[0,P1(j,1)],'YData',[0,-P1(j,2)]);
            set(pendule_masse2,'XData',P2(j,1),'YData',-P2(j,2));
            set(pendule_tige2,'XData',[P1(j,1),P2(j,1)],'YData',[-P1(j,2),-P2(j,2)]);
            grid on 
            if j> 200
            set(pendule_traj,'XData',P2(j-100:j,1),'YData',-P2(j-100:j,2));
            end
            if j<200
             set(pendule_traj,'XData',P2(1:j,1),'YData',-P2(1:j,2));
            end 
%     plot(P2(j,1),-P2(j,2),'.b','Markersize',5)
    drawnow
            
          end 
%%
    

     figure(13);
     set(figure(13),'position',[20 80 scz(3)/4 scz(4)/3]);
     
     xlabel('t(s)')
     ylabel('Energie(J)')
     title('Evolution des énergies cinétiques en fonction du temps');
     %legend('Energie cinetique pendule 1','Energie cinetique pendule 2','Energie cinetique systeme');
        
          if contexte==1 %position initiale
            Energie_cin1=plot(t,Ec1,'k.','Markersize',10,'Color','black');
            hold on;
            Energie_cin2=plot(t,Ec2,'k.','Markersize',10,'Color','blue');
            hold on;
            Energie_cint=plot(t,Ect,'k.','Markersize',10,'Color','green');
            hold on;
            axis([0 tf 0 maxt]);

          elseif contexte==2 %actualisé

            set(Energie_cin1,'XData',t,'YData',Ec1);
            set(Energie_cin2,'XData',t,'YData',Ec2);
            set(Energie_cint,'XData',t,'YData',Ect);
            
            plot(t,Ec1,'k.','Markersize',3,'Color','black'); 
            plot(t,Ec2,'k.','Markersize',3,'Color','blue'); 
            plot(t,Ect,'k.','Markersize',3,'Color','green'); 
            
            axis([0 tf 0 maxt]);
            grid on
          end
          
%%
          
     figure(14);
     set(figure(14),'position',[20 220+4*scz(4)/12 scz(3)/4 scz(4)/3]);
     
     
     xlabel('t(s)')
     ylabel('Energie(J)')
     title('Evolution des energies potentielles en fonction du temps');
     %legend('Energie potentielle pendule 1','Energie potentielle pendule 2','Energie potentielle systeme');
     
          if contexte==1 %position initiale
            Energie_pot1=plot(t,Ep1,'k.','Markersize',10,'Color','black');
            hold on;
            Energie_pot2=plot(t,Ep2,'k.','Markersize',10,'Color','blue');
            hold on;
            Energie_pott=plot(t,Ept,'k.','Markersize',10,'Color','green');
            hold on;
            axis([0 tf minp max34]);

          elseif contexte==2 %actualise

            set(Energie_pot1,'XData',t,'YData',Ep1);
            set(Energie_pot2,'XData',t,'YData',Ep2);
            set(Energie_pott,'XData',t,'YData',Ept);
            plot(t,Ept,'k.','Markersize',3,'Color','green'); 
            plot(t,Ep1,'k.','Markersize',3,'Color','black'); 
            plot(t,Ep2,'k.','Markersize',3,'Color','blue'); 
            axis([0 tf minp max34]);
            grid on
          end    
     
%%
     
     figure(15);
     set(figure(15),'position',[3*scz(3)/4-10 220+4*scz(4)/12 scz(3)/4 scz(4)/3]);
     Et=Ept+Ect;
     xlabel('t(s)')
     ylabel('Energie(J)')
     title('Evolution de l energie mecanique en fonction du temps');
     
     
      
     
         if contexte==1 %position initiale
                Energie_t=plot(t,Et,'k.','Markersize',10,'Color','blue');
                hold on;
                axis([0 tf mintot-1 maxtot+1]);
                
         elseif contexte==2 %actualise

                set(Energie_t,'XData',t,'YData',Et);
                plot(t,Et,'k.','Markersize',3,'Color','red'); 
                axis([0 tf mintot-1 maxtot+1]);
                grid on
         end    
         
end
