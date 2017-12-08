function zdot=Puma_stateder(z,taud)

q1=z(1);
q2=z(2);
q3=z(3);
q1dot=z(4);
q2dot=z(5);
q3dot=z(6);

TH0=2.8123;
TH1=2.2623;
TH2=-0.0066;
TH3=0.3453;
TH4=-0.085;
TH5=8.8038;
TH6=1.0711;
TH7=-0.9748;
TH8=6.8067;
TH9=55.35;

D(1,1)=TH2*cos(q2 + q3)^2 + 2*TH3*cos(q2 + q3)*cos(q2) + TH1*cos(q2)^2 + TH0;
D(1,2)=TH4*sin(q2 + q3) + TH7*sin(q2);
D(1,3)=TH4*sin(q2 + q3);
D(2,1)=D(1,2);
D(2,2)=TH5 + 2*TH3*cos(q3);
D(2,3)=TH6 + TH3*cos(q3);
D(3,1)=D(1,3);
D(3,2)=D(2,3);
D(3,3)=TH6;

C(1,1)=-q2dot*(TH2*cos(q2 + q3)*sin(q2 + q3) + TH3*cos(q2 + q3)*sin(q2) + TH3*sin(q2 + q3)*cos(q2) + TH1*cos(q2)*sin(q2)) - q3dot*(TH2*cos(q2 + q3)*sin(q2 + q3) + TH3*sin(q2 + q3)*cos(q2));
C(1,2)=q2dot*(TH4*cos(q2 + q3) + TH7*cos(q2)) - q1dot*(TH2*cos(q2 + q3)*sin(q2 + q3) + TH3*cos(q2 + q3)*sin(q2) + TH3*sin(q2 + q3)*cos(q2) + TH1*cos(q2)*sin(q2)) + TH4*q3dot*cos(q2 + q3);
C(1,3)=TH4*q2dot*cos(q2 + q3) - q1dot*(TH2*cos(q2 + q3)*sin(q2 + q3) + TH3*sin(q2 + q3)*cos(q2)) + TH4*q3dot*cos(q2 + q3);
C(2,1)=q1dot*(TH2*cos(q2 + q3)*sin(q2 + q3) + TH3*cos(q2 + q3)*sin(q2) + TH3*sin(q2 + q3)*cos(q2) + TH1*cos(q2)*sin(q2));
C(2,2)=-TH3*q3dot*sin(q3);
C(2,3)=-TH3*q2dot*sin(q3) -TH3*q3dot*sin(q3);
C(3,1)=q1dot*(TH2*cos(q2 + q3)*sin(q2 + q3) + TH3*sin(q2 + q3)*cos(q2));
C(3,2)=TH3*q2dot*sin(q3);
C(3,3)=0;

gg=[0;-TH8*cos(q2 + q3)-TH9*cos(q2);-TH8*cos(q2 + q3)];

X=[q1;q2;q3;q1dot;q2dot;q3dot];

Bemf=diag([38.7871;114.0474;28.5225]);


zdot=[X(4:6);inv(D)*(taud-C*X(4:6)-Bemf*X(4:6)-gg)];
