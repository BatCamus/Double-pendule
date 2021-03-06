function [dFX, dFdX] = calc_dFnl_Double_Pendule(X,dX)

global m1 m2 l1 l2 g 

dFX=zeros(length(X));           % D�riv�e de f par rapport � x
dFdX=zeros(length(X));          % D�riv�e de f par rapport � xpoint

%% D�claration variables
c=cos(X(1)-X(2));
s=sin(X(1)-X(2));

    % variables Fnl1
A1=m2*l1*c*s*dX(1)^2;
B1=m2*l2*s*dX(2)^2;
C1=m1*g*sin(X(1));
D1=m2*g*s*cos(X(2));
E1=m1*l1+m2*l1*(s)^2;
    
    % variables Fnl2
A2=-(m1+m2)*l1*s*dX(1)^2;
B2=-m2*l2*c*s*dX(2)^2;
C2=-(m1+m2)*g*c*sin(X(1));
D2=m1*l2+m2*l2*s^2;
E2=(m1+m2)*g*sin(X(2));

    % variables de dFnl1/X(1)
dA11=m2*l1*dX(1)^2*(-s^2+c^2);
dB11=m2*l2*dX(2)^2*c;
dC11=m1*g*cos(X(1));
dD11=m2*g*c*sin(X(2));
dE11=m2*l1*sin(2*(X(1)-X(2)));
    
    % variables de dFnl1/X(2)
dA12=m2*l1*dX(1)^2*(s^2-c^2);
dB12=-m2*l2*dX(2)^2*c;
dC12=0;
dD12=-m2*g*cos(X(1)-2*X(2));
dE12=-m2*l1*sin(2*(X(1)-X(2)));
    
    % variables dFnl2/X(1)
dA21=-(m1+m2)*l1*dX(1)^2*c;
dB21=m2*l2*dX(2)^2*(s^2-c^2);
dC21=-(m1+m2)*g*c;
dD21=m2*l2*sin(2*(X(1)-X(2)));
dE21=0;

    % variables dFnl2/X(2)
dA22=(m1+m2)*l1*dX(1)^2*cos(X(1)-X(2));
dB22=m2*l2*dX(2)^2*(cos(X(1)-X(2))^2-sin(X(1)-X(2))^2);
dC22=-(m1+m2)*g*sin(X(1))*sin(X(1)-X(2));
dD22=-m2*l2*sin(2*(X(1)-X(2)));
dE22=(m1+m2)*g*cos(X(2));



dFdX(1,1)=((m2*l1*cos(X(1)-X(2))*sin(X(1)-X(2)))/E1^2)*2*dX(1);
dFdX(1,2)=((m2*l2*sin(X(1)-X(2)))/E1)*2*dX(2);
dFdX(2,1)=(-(m1+m2)*l1*sin(X(1)-X(2))/D2)*2*dX(1);
dFdX(2,2)=((-m2*l2*cos(X(1)-X(2))*sin(X(1)-X(2)))/D2)*2*dX(2);

dFX(1,1)=(E1*(dA11+dB11+dC11+dD11)-(A1+B1+C1+D1)*dE11)/E1^2;
dFX(1,2)=(E1*(dA12+dB12+dC12+dD12)-(A1+B1+C1+D1)*dE12)/E1^2;
dFX(2,1)=(D2*(dA21+dB21+dC21)-dD21*(A2+B2+C2))/D2^2;
dFX(2,2)=(D2*(dA22+dB22+dC22+dE22)-dD22*(A2+B2+C2+E2))/D2^2;

end