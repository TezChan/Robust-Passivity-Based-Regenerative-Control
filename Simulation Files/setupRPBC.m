%setupRPBC

%Calculate nominal parameters for use in controller

%These parameters match the coordinates and schematic used in MCE647,
%Spring 15

%Refer to midterm project schematic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m2=19.0157;  %mass of link 2, kg
m46=1.1720;  %combined mass of links 4-6, kg
m3=5.5340;    %mass of link 3, kg
d1=0.666;       %height of o1 from floor, m
d2=0.2435;     %link 2 offset, m
d3=0.0934;     %link 3 offset, m
a2=0.4318;     %distance between link 2 and link 3 axes, m
c2x=-0.3447;  %best estimate of CM2 location, m
c3x=0.1410;   %best estimate of CM3 location, m

%Principal moments of inertia at CM-centered frames, kg-m^2
%Refer to frame orientations in sketch!

I1y=1.3943;   
I2x=0.13;        
I2z=5.2532;
I2y=0.5432;
I3x=0.1860;
I3y=0.1181;
I3z=1.0770;
g=9.81;          %acceleration of gravity, m/s^2

%Theta parameters:
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

%Nominal (for use in simulated and realtime controller):
TH0=[TH0;TH1;TH2;TH3;TH4;TH5;TH6;TH7;TH8;TH9];


%Perturbed (for use in plant simulation):
level=0.2;  %perturbation level

m2p=m2+m2*(1-2*rand)*level;
m46p=m46+m46*(1-2*rand)*level;
m3p=m3+m3*(1-2*rand)*level;
d1p=d1+d1*(1-2*rand)*level;
d2p=d2+d2*(1-2*rand)*level;
d3p=d3+d3*(1-2*rand)*level;
a2p=a2+a2*(1-2*rand)*level;
c2xp=c2x+c2x*(1-2*rand)*level;
c3xp=c3x+c3x*(1-2*rand)*level;
I1yp=I1y+I1y*(1-2*rand)*level;
I2xp=I2x+I2x*(1-2*rand)*level;
I2zp=I2z+I2z*(1-2*rand)*level;
I2yp=I2y+I2y*(1-2*rand)*level;
I3xp=I3x+I3x*(1-2*rand)*level;
I3yp=I3y+I3y*(1-2*rand)*level;
I3zp=I3z+I3z*(1-2*rand)*level;

TH0p=(m2p + m3p)*d2p^2 - 2*m3p*d2p*d3p + m3p*d3p^2 + I1yp + I2xp + I3xp; %OK
TH1p=m2p*c2xp^2 + 2*m2p*c2xp*a2p + (m2p + m3p)*a2p^2 - I2xp + I2yp; %OK
TH2p=m3p*c3xp^2 - I3xp + I3yp; %OK
TH3p=c3xp*a2p*m3p; %OK
TH4p=c3xp*m3p*(d3p-d2p); %OK
TH5p=m2p*c2xp^2 + 2*m2p*c2xp*a2p + m3p*c3xp^2 + (m2p + m3p)*a2p^2 + I2zp + I3zp; %OK
TH6p=I3zp+m3p*c3xp^2; %OK
TH7p=-d2p*m2p*(c2xp + a2p) + a2p*m3p*(d3p - d2p); %OK
TH8p=c3xp*g*m3p; %OK
TH9p=g*m2p*(c2xp + a2p) + a2p*g*m3p; %OK

THpert=[TH0p;TH1p;TH2p;TH3p;TH4p;TH5p;TH6p;TH7p;TH8p;TH9p];

%Calculate parameter error bound
rhob=norm(THpert-TH0)
%Percentage:
rhob/norm(TH0)

%These are the effective input constants capturing 
%gear ratio, motor torque constant and amplifier gain
systemident
%Units: Nm/V
Ku=diag(A_Matrix);

%Simulation settings
%Initial condition for robot
Z0=[0;0;0;0;-pi/2;0];  %reference position after manual control
%Trajectories: (amplitude, frequency and bias)
A1=1;w1=2;b1=0;
A2=0.25;w2=2;b2=0;
A3=0.5;w3=2;b3=-pi/2;

