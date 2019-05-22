function dxdt=Pendule_Double_Non_Lin(t,x)

    global g l1 l2 m1 m2
    
    dxdt=zeros(4,1); 
    dxdt(1)=x(2); 
    dxdt(2)=((-m2*l1*sin(x(1)-x(3))*cos(x(1)-x(3))*x(2)^2)-(m2*l2*sin(x(1)-x(3))*x(4)^2)-(m1*g*sin(x(1)))-(m2*g*sin(x(1)-x(3))*cos(x(3))))/(m1*l1+m2*l1*sin(x(1)-x(3))^2);
    dxdt(3)= x(4); 
    dxdt(4)=(((m1+m2)*l1*sin(x(1)-x(3))*x(2)^2)+(m2*l2*sin(x(1)-x(3))*cos(x(1)-x(3))*x(4)^2)+((m1+m2)*g*sin(x(1)-x(3))*cos(x(1))))/(m1*l2+m2*l2*sin(x(1)-x(3))^2);
    
end