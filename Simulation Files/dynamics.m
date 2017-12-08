%dynamics.m
%Develops dynamic equations for the PUMA560 (first 3 links)
% DHkinematics.m must be in the same directory
%Hanz Richter, CSU 2015

%Find kinematic transformations by DH
DHkinematics3

%Find the Jacobians to the centers of mass
%To first link CM: J1
o1=[0;0;d1]; %We are assuming that the CM of L1 is at o1, it's irrelevant 
z0=[0;0;1];
J1v1=cross(z0,o1);
J1v2=[0;0;0]; %the velocity of link n is independent of the motion of higher links
J1v3=[0;0;0];
J1v=[J1v1 J1v2 J1v3];
J1w1=z0;
J1w2=[0;0;0];
J1w3=[0;0;0];
J1w=[J1w1 J1w2 J1w3];

%To second link CM: J2
%Assume that the CM of L2 is c22=[c2x;0;0]
%Note that c2x will be a negative parameter
syms c2x
c2=simplify(A20*[c2x;0;0;1]); 
c2=eval(c2(1:3));
o2=simplify(A20*[0;0;0;1]);
o2=eval(o2(1:3));

z1=R1*[0;0;1];
J2v1=cross(z0,c2);
J2v2=cross(z1,c2-o1);
J2v3=[0;0;0]; %the velocity of link n is independent of the motion of higher links
J2w1=z0;
J2w2=z1;
J2w3=[0;0;0];
J2v=[J2v1 J2v2 J2v3];
J2w=[J2w1 J2w2 J2w3];

%To third link CM: J3
%Assume that the CM of L3 is c33=[c3x;0;0]
%Note that c3x will be a negative parameter
syms c3x
c3=simplify(A30*[c3x;0;0;1]); 
c3=eval(c3(1:3));

z2=z1;
J3v1=cross(z0,c3);
J3v2=cross(z1,c3-o1);
J3v3=cross(z2,c3-o2);
J3w1=z0;
J3w2=z1;
J3w3=z2;
J3v=[J3v1 J3v2 J3v3];
J3w=[J3w1 J3w2 J3w3];

%Inertia matrix derivations
syms I1x I1y I1z I2x I2y I2z I3x I3y I3z m1 m2 m3
I1=diag([I1x I1y I1z]);
I2=diag([I2x I2y I2z]);
I3=diag([I3x I3y I3z]);

sum=m1*J1v.'*J1v+J1w.'*R1*I1*R1.'*J1w;
sum=simplify(sum+m2*J2v.'*J2v+J2w.'*R2*I2*R2.'*J2w);
sum=simplify(sum+m3*J3v.'*J3v+J3w.'*R3*I3*R3.'*J3w);
%c2x=0;
D_q=eval(sum);

%Coriolis/Centripetal Matrix
q=[q1;q2;q3];
for i=1:3,
    for j=1:3,
        for k=1:3,
            c(i,j,k)=(1/2)*((diff(D_q(k,j),q(i)))+(diff(D_q(k,i),q(j)))-(diff(D_q(i,j),q(k))));
        end
      end
end

syms q1dot q2dot q3dot

for k=1:3,
    for j=1:3,
        C_q(k,j)=c(1,j,k)*q1dot+c(2,j,k)*q2dot+c(3,j,k)*q3dot;
    end
end

%The above has been checked to satisfy the skew-symmetry property for
%Ddot-2C

syms g
gv=[0;0;-g];

P=-simplify(m1*gv.'*o1+m2*gv.'*c2+m3*gv.'*c3);

gq=eval([diff(P,q1);diff(P,q2);diff(P,q3)]);



