%DHkinematics3.m
%Develops forward kinematics equations for the PUMA560 (first 3 joints)
%Hanz Richter, CSU 2015

syms q1 q2 q3 d1 d2 d3 a2

%Refer to schematic. 

%Develop transformations
Rxm90=[1 0 0 0; 0 0 1 0;0 -1 0 0; 0 0 0 1]; %Rotx,-90
Rxp90=[1 0 0 0; 0 0 -1 0;0 1 0 0; 0 0 0 1]; %Rotx,+90
%Note use +/-1 instead of cos(-pi/2) and sin(-pi/2) to avoid numerical representation problems (vpa) 

%A10: (Rotz,q1)(Transz,d1)(Rotx,-pi/2)
A10=[cos(q1) -sin(q1) 0 0; sin(q1) cos(q1) 0 0; 0 0 1 d1; 0 0 0 1]*Rxm90;
%A21: (Rotz,q2)(Transz,-d2)(Transx,a2)
A21=[cos(q2) -sin(q2) 0 0; sin(q2) cos(q2) 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0;0 0 1 -d2; 0 0 0 1]*[1 0 0 a2; 0 1 0 0;0 0 1 0; 0 0 0 1];
A20=simplify(A10*A21);
%A32: (Rotz,q3)(Transz,d3)
A32=[cos(q3) -sin(q3) 0 0; sin(q3) cos(q3) 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0;0 0 1 d3; 0 0 0 1];
A30=simplify(A20*A32);
R1=simplify(A10(1:3,1:3));
R2=simplify(A20(1:3,1:3));
R3=simplify(A30(1:3,1:3));


%Verified to work on 3/2/15
