clc
clear

syms q1 q2 q3 q1dot q2dot q3dot q1ddot q2ddot q3ddot

syms TH0 TH1 TH2 TH3 TH4 TH5 TH6 TH7 TH8 TH9

syms taud1 taud2 taud3


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

taud=[taud1;taud2;taud3];
%Bemf=[38.7871;114.0474;28.5225];
syms Bemf1 Bemf2 Bemf3
Bemf=diag([Bemf1;Bemf2;Bemf3]);

Xdot=simplify([X(4:6);inv(D)*(taud-C*X(4:6)-Bemf*X(4:6)-gg)]);

syms R1 a1 R2 a2 R3 a3
Q=zeros(6,6);
M=[zeros(3,3);-1*ones(3,3)];
W=diag([2*R1/(a1^2),2*R2/(a2^2),2*R3/(a3^2)]);

G=0.5*[X.',taud.']*[Q,M;M.',W]*[X;taud];
syms p1 p2 p3 p4 p5 p6
p=[p1;p2;p3;p4;p5;p6];
H=G+p.'*Xdot;







