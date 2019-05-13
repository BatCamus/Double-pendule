%% Jacobienne Pendule double non linéaire Newton raphson

function J=Jacobienne_double_pendule(z,delta,f1,f2)
%
global m1 m2 l1 l2 g;
%
f1theta1=(m1+m2)*l1*z(3)+m2*l2*cos(z(1)+delta-z(4))*z(6)+m2*l2*sin(z(1)+delta-z(4))*z(5)^2+(m1+m2)*g*sin(z(1)+delta);
f1dtheta1=(m1+m2)*l1*z(3)+m2*l2*cos(z(1)-z(4))*z(6)+m2*l2*sin(z(1)-z(4))*z(5)^2+(m1+m2)*g*sin(z(1));
f1ddtheta1=(m1+m2)*l1*(z(3)+delta)+m2*l2*cos(z(1)-z(4))*z(6)+m2*l2*sin(z(1)-z(4))*z(5)^2+(m1+m2)*g*sin(z(1));

f1theta2=(m1+m2)*l1*z(3)+m2*l2*cos(z(1)-z(4)-delta)*z(6)+m2*l2*sin(z(1)-z(4)-delta)*z(5)^2+(m1+m2)*g*sin(z(1));
f1dtheta2=(m1+m2)*l1*z(3)+m2*l2*cos(z(1)-z(4))*z(6)+m2*l2*sin(z(1)-z(4))*(z(5)+delta)^2+(m1+m2)*g*sin(z(1));
f1ddtheta2=(m1+m2)*l1*z(3)+m2*l2*cos(z(1)-z(4))*(z(6)+delta)+m2*l2*sin(z(1)-z(4))*z(5)^2+(m1+m2)*g*sin(z(1));

f2theta1=m2*l2*z(6)+m2*l1*cos(z(1)+delta-z(4))*z(3)-m2*l1*sin(z(1)+delta-z(4))*z(2)^2+m2*g*sin(z(4));
f2dtheta1=m2*l2*z(6)+m2*l1*cos(z(1)-z(4))*z(3)-m2*l1*sin(z(1)-z(4))*(z(2)+delta)^2+m2*g*sin(z(4));
f2ddtheta1=m2*l2*z(6)+m2*l1*cos(z(1)-z(4))*(z(3)+delta)-m2*l1*sin(z(1)-z(4))*z(2)^2+m2*g*sin(z(4));%

f2theta2=m2*l2*z(6)+m2*l1*cos(z(1)-z(4)-delta)*z(3)-m2*l1*sin(z(1)-z(4)-delta)*z(2)^2+m2*g*sin(z(4)-delta);
f2dtheta2=m2*l2*z(6)+m2*l1*cos(z(1)-z(4))*z(3)-m2*l1*sin(z(1)-z(4))*z(2)^2+m2*g*sin(z(4));
f2ddtheta2=m2*l2*(z(6)+delta)+m2*l1*cos(z(1)-z(4))*z(3)-m2*l1*sin(z(1)-z(4))*z(2)^2+m2*g*sin(z(4));

J11=(f1theta1-f1)/delta;
J12=(f1dtheta1-f1)/delta;
J13=(f1ddtheta1-f1)/delta;
J14=(f1theta2-f1)/delta;
J15=(f1dtheta2-f1)/delta;
J16=(f1ddtheta2-f1)/delta;
J21=(f2theta1-f2)/delta;
J22=(f2dtheta1-f2)/delta;
J23=(f2ddtheta1-f2)/delta;
J24=(f2theta2-f2)/delta;
J25=(f2dtheta2-f2)/delta;
J26=(f2ddtheta2-f2)/delta;
%
J=[J11 J12 J13 J14 J15 J16;J21 J22 J23 J24 J25 J26];
%
end