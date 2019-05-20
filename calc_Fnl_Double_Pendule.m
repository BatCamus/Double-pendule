function Fnl = calc_Fnl_Double_Pendule(X,dX)

global m1 m2 l1 l2 g 
ODE=0;
Fnl=zeros(size(X));
%
c=cos(X(1)-X(2));
s=sin(X(1)-X(2));

if(ODE==0)
Fnl(1)=(m2*l1*c*s*(dX(1))^2+m2*l2*s*(dX(2))^2+m1*g*sin(X(1))+m2*g*s*cos(X(2)))/(m1*l1+m2*l1*(s)^2);
Fnl(2)=-((m1+m2)*l1*s*dX(1)^2+m2*l2*c*s*dX(2)^2+(m1+m2)*g*c*sin(X(1))-(m1+m2)*g*sin(X(2)))/(m1*l2+m2*l2*s^2);
end

% Formule trouvé sur internet
if(ODE==1)
Fnl(1)=-((-m2*l1*sin(X(1)-X(2))*cos(X(1)-X(2))*dX(1)^2)-(m2*l2*sin(X(1)-X(2))*dX(2)^2)-(m1*g*sin(X(1)))-(m2*g*sin(X(1)-X(2))*cos(X(2))))/(m1*l1+m2*l1*sin(X(1)-X(2))^2);
Fnl(2)=-(((m1+m2)*l1*sin(X(1)-X(2))*dX(1)^2)+(m2*l2*sin(X(1)-X(2))*cos(X(1)-X(2))*dX(2)^2)+((m1+m2)*g*sin(X(1)-X(2))*cos(X(1))))/(m1*l2+m2*l2*sin(X(1)-X(2))^2);
end
end