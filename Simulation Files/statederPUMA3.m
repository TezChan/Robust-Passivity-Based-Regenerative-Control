function zdot=statederPUMA3(z,u,a,B_e)
%State derivatives for the first three links of the PUMA560


%These values come from perturbed nominal parameters in setupRPBC
TH0=3.1754;
TH1=1.9290;
TH2=0.0936;
TH3=0.4106;
TH4=-0.1478;
TH5=7.7600;
TH6=1.2239;
TH7=-0.9787;
TH8=8.5294;
TH9=46.0321;


%These are the effective input constants capturing 
%gear ratio, motor torque constant and amplifier gain

%Units: Nm/V
Ku=diag(a);
T=Ku*u-diag(B_e)*z(4:6);
q1=z(1);
q2=z(2);
q3=z(3);
q1dot=z(4);
q2dot=z(5);
q3dot=z(6);

%find numerical values for matrices M, C and calculate state derivatives
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
zdot=[z(4:6); inv(D)*(T-C*z(4:6)-gg)];