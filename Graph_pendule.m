function Graph_Pendule(contexte,P11,P12,P21,P22,l1,l2,t,tf,Ec1,Ec2,Ep1,Ep2,maxt,mint,max34,minp,maxtot,mintot)

global pendule_masse1 pendule_tige1 pendule_masse2 pendule_tige2 Energie_cin1 Energie_cin2 Energie_cint Energie_pott Energie_t Energie_pot1 Energie_pot2 pendule_traj


Ect=Ec1+Ec2;
Ept=Ep1+Ep2;
scz=get(0,'screensize'); %Taille écran
    
%%
        figure(1);

        set(figure(1),'position',[20+scz(3)/3 220+4*scz(4)/12 scz(3)/3 scz(4)/3]);
        
      
        xlabel('X(dm)')
        ylabel('Z(dm)')
        title('Mouvement d un double pendule');
       

          if contexte==1 %position initiale

            pendule_masse1=plot(P11,-P12,'k.','MarkerSize',35,'Color','red');
            hold on
            pendule_tige1=plot([0,P11],[0,-P12],'LineWidth',2,'Color','black');
            hold on
            pendule_masse2=plot(P21,-P22,'k.','MarkerSize',35,'Color','red');
            hold on
            pendule_tige2=plot([P11,P21],[-P12,-P22],'LineWidth',2,'Color','blue');
            hold on
            
    
            axis([-(l1+l2) (l1+l2) -1.2*(l1+l2) 6]);

          elseif contexte==2 %position actualisée

            set(pendule_masse1,'XData',P11,'YData',-P12);
            set(pendule_tige1,'XData',[0,P11],'YData',[0,-P12]);
            set(pendule_masse2,'XData',P21,'YData',-P22);
            set(pendule_tige2,'XData',[P11,P21],'YData',[-P12,-P22]);
            plot(P21,-P22,'.r','Markersize',0.5);
            grid on 
            
          end 
%%
    

     figure(2);
     set(figure(2),'position',[20 80 scz(3)/4 scz(4)/3]);
     
     xlabel('t(s)')
     ylabel('Energie(J)')
     title('Evolution des énergies cinétiques en fonction du temps');
     legend('Energie cinétique pendule 1','Energie cinétique pendule 2','Energie cinétique système');
        
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
          
     figure(3);
     set(figure(3),'position',[20 220+4*scz(4)/12 scz(3)/4 scz(4)/3]);
     
     
     xlabel('t(s)')
     ylabel('Energie(J)')
     title('Evolution des énergies potentielles en fonction du temps');
     legend('Energie potentielle pendule 1','Energie potentielle pendule 2','Energie potentielle système');
     
          if contexte==1 %position initiale
            Energie_pot1=plot(t,Ep1,'k.','Markersize',10,'Color','black');
            hold on;
            Energie_pot2=plot(t,Ep2,'k.','Markersize',10,'Color','blue');
            hold on;
            Energie_pott=plot(t,Ept,'k.','Markersize',10,'Color','green');
            hold on;
            axis([0 tf minp max34]);

          elseif contexte==2 %actualisé

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
     
     figure(4);
     set(figure(4),'position',[3*scz(3)/4-10 220+4*scz(4)/12 scz(3)/4 scz(4)/3]);
     Et=Ept+Ect;
     xlabel('t(s)')
     ylabel('Energie(J)')
     title('Evolution de l énergie mécanique en fonction du temps');
     
     
      
     
         if contexte==1 %position initiale
                Energie_t=plot(t,Et,'k.','Markersize',10,'Color','blue');
                hold on;
                axis([0 tf mintot-1 maxtot+1]);
                
         elseif contexte==2 %actualisé

                set(Energie_t,'XData',t,'YData',Et);
                plot(t,Et,'k.','Markersize',3,'Color','red'); 
                axis([0 tf mintot-1 maxtot+1]);
                grid on
         end    
         
end
