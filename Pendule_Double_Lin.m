%% Fonction utilisée dans ode45
function dxdt=Pendule_Double_Lin(t,x)

    global mu g l1 l2
    
    dxdt=zeros(4,1); 
    
    dxdt(1)=x(2); 
    dxdt(2)=(mu*g*x(3)-(1+mu)*g*x(1))/l1;
    
    dxdt(3)= x(4); 
    dxdt(4)=((1+mu)g*x(1)-(1+mu)g*x(3))/l2;
    
end