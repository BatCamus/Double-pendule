function Fnl = calc_Fnl_Double_Pendule(X,dX,t)

global m1 m2 l1 l2 g mu M

Fnl=zeros(size(X));
%
Fnl(1)=(m2*l1*cos(X(1)-X(2))*sin(X(1)-X(2))*(dX(1))^2+m2*l2*sin(X(1)-X(2))*(dX(2))^2+m1*g*sin(X(1))+m2*g*cos(X(1)-X(2))*sin(X(2)))/(m1*l1+m2*l1*sin(X(1)-X(2))^2);
Fnl(2)=((m1+m2)*l1*sin(X(1)-X(2))*dX(1)^2+m2*l2*cos(X(1)-X(2))*sin(X(1)-X(2))*dX(2)^2+(m1+m2)*g*sin(X(1)-X(2))*cos(X(1)))/(m1*l2+m2*l2*sin(X(1)-X(2))^2);

end

