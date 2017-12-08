syms q1 q2 q3 q1dot q2dot q3dot q1ddot q2ddot q3ddot
Yav1=[q1ddot, q1ddot*cos(q2)^2 - q1dot*q2dot*cos(q2)*sin(q2) - q2dot*q1dot*cos(q2)*sin(q2), q1ddot*cos(q2 + q3)^2 - q1dot*(q2dot*cos(q2 + q3)*sin(q2 + q3) + q3dot*cos(q2 + q3)*sin(q2 + q3)) - q1dot*q2dot*cos(q2 + q3)*sin(q2 + q3) - q1dot*q3dot*cos(q2 + q3)*sin(q2 + q3), 2*q1ddot*cos(q2 + q3)*cos(q2) - q1dot*(q2dot*(cos(q2 + q3)*sin(q2) + sin(q2 + q3)*cos(q2)) + q3dot*sin(q2 + q3)*cos(q2)) - q1dot*q2dot*(cos(q2 + q3)*sin(q2) + sin(q2 + q3)*cos(q2)) - q1dot*q3dot*sin(q2 + q3)*cos(q2), q2dot*(q2dot*cos(q2 + q3) + q3dot*cos(q2 + q3)) + q3dot*(q2dot*cos(q2 + q3) + q3dot*cos(q2 + q3)) + q2ddot*sin(q2 + q3) + q3ddot*sin(q2 + q3), 0, 0, q2ddot*sin(q2) + q2dot*q2dot*cos(q2), 0, 0,38.7871*q1dot];
Yav2=[ 0, q1dot*q1dot*cos(q2)*sin(q2), q1dot*q1dot*cos(q2 + q3)*sin(q2 + q3), 2*q2ddot*cos(q3) + q3ddot*cos(q3) - q3dot*(q2dot*sin(q3) + q3dot*sin(q3)) - q3dot*q2dot*sin(q3) + q1dot*q1dot*(cos(q2 + q3)*sin(q2) + sin(q2 + q3)*cos(q2)), q1ddot*sin(q2 + q3), q2ddot, q3ddot, q1ddot*sin(q2), -cos(q2 + q3), -cos(q2),q2dot*114.0474];
Yav3=[ 0, 0, q1dot*q1dot*cos(q2 + q3)*sin(q2 + q3), q2ddot*cos(q3) + q2dot*q2dot*sin(q3) + q1dot*q1dot*sin(q2 + q3)*cos(q2), q1ddot*sin(q2 + q3), 0, q2ddot + q3ddot, 0, -cos(q2 + q3), 0,q3dot*28.5225];

Yav=[Yav1;Yav2;Yav3];
A=[1/8.1797,1/14.0261,1//7.0144];%alpha*n^2 of each motor
%Nominal parameters 
TH0 =[
    2.8123
    2.2623
    -0.0066
    0.3453
   -0.08499
    8.8038
    1.0711
   -0.9748
    6.8067
   55.35
   1];
%Tau=Yav*TH0