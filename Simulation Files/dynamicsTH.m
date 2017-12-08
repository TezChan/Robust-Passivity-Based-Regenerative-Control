%dynamics.m
%Develops dynamic equations for the PUMA560 (first 3 links)
% Parameterized version
%Hanz Richter, CSU 2015



syms q1 q2 q3
syms TH0 TH1 TH2 TH3 TH4 TH5 TH6 TH7 TH8 TH9
D(1,1)=TH0+TH1*cos(q2)^2+TH2*cos(q2 + q3)^2+2*TH3*cos(q2+q3)*cos(q2);
D(1,2)=TH4*sin(q2+q3)+TH7*sin(q2);
D(1,3)=TH4*sin(q2+q3);
D(2,1)=D(1,2);
D(2,2)=2*TH3*cos(q3)+TH5;
D(2,3)=TH6+TH3*cos(q3);
D(3,1)=D(1,3);
D(3,2)=D(2,3);
D(3,3)=TH6;


%Coriolis/Centripetal Matrix
q=[q1;q2;q3];
for i=1:3,
    for j=1:3,
        for k=1:3,
            c(i,j,k)=(1/2)*((diff(D(k,j),q(i)))+(diff(D(k,i),q(j)))-(diff(D(i,j),q(k))));
        end
      end
end

syms q1dot q2dot q3dot

for k=1:3,
    for j=1:3,
        C(k,j)=c(1,j,k)*q1dot+c(2,j,k)*q2dot+c(3,j,k)*q3dot;
    end
end

%The above has been checked to satisfy the skew-symmetry property for
%Ddot-2C

gTH=[0;-TH8*cos(q2+q3)-TH9*cos(q2);-TH8*cos(q2+q3)];


%Extract regressor
syms q1ddot q2ddot q3ddot
aux=D*[q1ddot;q2ddot;q3ddot]+C*[q1dot;q2dot;q3dot]+gTH;
TH=[TH0 TH1 TH2 TH3 TH4 TH5 TH6 TH7 TH8 TH9].';

for i=1:3, %number of rows (links)
    for j=1:10, %number of parameters
        Y(i,j)=diff(aux(i),TH(j));
    end
end

%Verification


syms c2x c3x d2 d3 a2 g
syms I1x I1y I1z I2x I2y I2z I3x I3y I3z m1 m2 m3

%Parameters
TH0=(m2 + m3)*d2^2 - 2*m3*d2*d3 + m3*d3^2 + I1y + I2x + I3x; %OK
TH1=m2*c2x^2 + 2*m2*c2x*a2 + (m2 + m3)*a2^2 - I2x + I2y; %OK
TH2=m3*c3x^2 - I3x + I3y; %OK
TH3=c3x*a2*m3; %OK
TH4=c3x*m3*(d3-d2); %OK
TH5=m2*c2x^2 + 2*m2*c2x*a2 + m3*c3x^2 + (m2 + m3)*a2^2 + I2z + I3z; %OK
TH6=I3z+m3*c3x^2; %OK
TH7=-d2*m2*(c2x + a2) + a2*m3*(d3 - d2); %OK
TH8=c3x*g*m3; %OK
TH9=g*m2*(c2x + a2) + a2*g*m3; %OK

aux_right=eval(Y*TH);
dynamics
aux_left=D_q*[q1ddot;q2ddot;q3ddot]+(C_q)*[q1dot;q2dot;q3dot]+gq; 

%Three zero answers should be obtained!
simplify(aux_left(1)-aux_right(1))
simplify(aux_left(2)-aux_right(2))
simplify(aux_left(3)-aux_right(3))