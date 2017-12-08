function u=RPBCpuma(z,zd)
%Robust passivity-based controller for the first 3 links of the PUMA 560

q1=z(1);
q2=z(2);
q3=z(3);
q1dot=z(4);
q2dot=z(5);
q3dot=z(6);

q1d=zd(1);
q2d=zd(2);
q3d=zd(3);
q1ddot=zd(4);
q2ddot=zd(5);
q3ddot=zd(6);
q1dddot=zd(7);
q2dddot=zd(8);
q3dddot=zd(9);

qtilde=[q1-q1d;q2-q2d;q3-q3d];
qtildedot=[q1dot-q1ddot;q2dot-q2ddot;q3dot-q3ddot];

%Gains
L=diag([2 10 2]);
K=diag([5 10 5]);

v=[q1ddot;q2ddot;q3ddot]-L*qtilde;
a=[q1dddot;q2dddot;q3dddot]-L*qtildedot;
r=qtildedot+L*qtilde;

a_1=a(1);a_2=a(2);a_3=a(3);
v_1=v(1);v_2=v(2);v_3=v(3);

%Form the regressor
Yav1=[a_1, a_1*cos(q2)^2 - q1dot*v_2*cos(q2)*sin(q2) - q2dot*v_1*cos(q2)*sin(q2), a_1*cos(q2 + q3)^2 - v_1*(q2dot*cos(q2 + q3)*sin(q2 + q3) + q3dot*cos(q2 + q3)*sin(q2 + q3)) - q1dot*v_2*cos(q2 + q3)*sin(q2 + q3) - q1dot*v_3*cos(q2 + q3)*sin(q2 + q3), 2*a_1*cos(q2 + q3)*cos(q2) - v_1*(q2dot*(cos(q2 + q3)*sin(q2) + sin(q2 + q3)*cos(q2)) + q3dot*sin(q2 + q3)*cos(q2)) - q1dot*v_2*(cos(q2 + q3)*sin(q2) + sin(q2 + q3)*cos(q2)) - q1dot*v_3*sin(q2 + q3)*cos(q2), v_2*(q2dot*cos(q2 + q3) + q3dot*cos(q2 + q3)) + v_3*(q2dot*cos(q2 + q3) + q3dot*cos(q2 + q3)) + a_2*sin(q2 + q3) + a_3*sin(q2 + q3), 0, 0, a_2*sin(q2) + q2dot*v_2*cos(q2), 0, 0];
Yav2=[ 0, q1dot*v_1*cos(q2)*sin(q2), q1dot*v_1*cos(q2 + q3)*sin(q2 + q3), 2*a_2*cos(q3) + a_3*cos(q3) - v_3*(q2dot*sin(q3) + q3dot*sin(q3)) - q3dot*v_2*sin(q3) + q1dot*v_1*(cos(q2 + q3)*sin(q2) + sin(q2 + q3)*cos(q2)), a_1*sin(q2 + q3), a_2, a_3, a_1*sin(q2), -cos(q2 + q3), -cos(q2)];
Yav3=[ 0, 0, q1dot*v_1*cos(q2 + q3)*sin(q2 + q3), a_2*cos(q3) + q2dot*v_2*sin(q3) + q1dot*v_1*sin(q2 + q3)*cos(q2), a_1*sin(q2 + q3), 0, a_2 + a_3, 0, -cos(q2 + q3), 0];

Yav=[Yav1;Yav2;Yav3];
%Nominal parameters 
TH0 =[
    2.9625
    1.5893
    0.0421
    0.3369
   -0.1171
    7.6163
    1.1870
   -0.7620
    7.6547
   39.6898];
   
deadzone=0.1;
rho=6;
s=Yav'*r;
if norm(s)>deadzone,
    dTh=-rho*s/norm(s);
else
    dTh=-rho*s/deadzone;
end
Th_hat=TH0+dTh;
%Effective input constant inverse (Volts/(Nm))
Kui=diag([0.0543 0.0806 0.1078]);
u=Kui*(Yav*Th_hat-K*r);
