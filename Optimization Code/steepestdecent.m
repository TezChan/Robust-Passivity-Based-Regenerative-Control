%Steepest Decent for Puma

clc
close all
clear



t0=0; %initial time
tf=2; %final time
dt=0.001; %time step
time=t0:dt:tf; %time vector
N=length(time); %number of time steps

%intial u
u0=zeros(3,N);

u=u0;
%initial conditions
x0=[0 pi/3 0 0 0 0]';
xf=[pi -pi/3 pi/3 0 0 0];
H=diag([10 10 10 10 10 10]);
J=inf;
iter=1;

%system parameters
R1=1.725;
R2=1.725;
R3=1.725;
a1=0.2254*62.61;
a2=0.2254*107.36;
a3=0.2254*53.69;

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

Q=zeros(6,6);
M=[zeros(3,3);-1*ones(3,3)];
W=diag([2*R1/(a1^2),2*R2/(a2^2),2*R3/(a3^2)]);
eps=[];
n=0
while 1
    
% system simulation
xx=x0;
% Initialize arrays for plotting at the end of the program
x= zeros(6, N);
% Begin simulation loop
for i=1:N-1
% Save data for plotting
x(:,i) = xx;
% Update x
xxdot = Puma_stateder(xx, u(:,i));
xxdot1 = Puma_stateder(xx+xxdot*dt,u(:,i+1));
xx = xx + (xxdot + xxdot1) * dt / 2;
end

p0=(-(xf'-xx)'*H)';

%costate simulation
pp=p0;
% Initialize arrays for plotting at the end of the program
p= zeros(6, N);
% Begin simulation loop
for i=1:N-1
% Save data for plotting
p(:,N+1-i) = pp;
% Update x
ppdot = Puma_costateder(pp,x(:,N+1-i), u(:,N+1-i));
ppdot1 = Puma_costateder(pp+ppdot*dt,x(:,N-i),u(:,N-i));
pp = pp + (ppdot + ppdot1) * dt / 2;
end



%Compute the partial of the Hamiltonian with respect to the control
parHparu(1,:)=(p(5,:).*(TH6.*TH7.*sin(x(2,:)) - TH3.*TH4.*sin(x(2,:) + x(3,:)).*cos(x(3,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) - x(5,:) - x(6,:) - x(4,:) - (p(6,:).*(TH4.*TH6.*sin(x(2,:) + x(3,:)) - TH4.*TH5.*sin(x(2,:) + x(3,:)) + TH6.*TH7.*sin(x(2,:)) - TH3.*TH4.*sin(x(2,:) + x(3,:)).*cos(x(3,:)) + TH3.*TH7.*cos(x(3,:)).*sin(x(2,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) + (p(4,:).*(TH3.^2.*cos(x(3,:)).^2 - TH5.*TH6 + TH6.^2))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) + (2.*R1.*u(1,:))./a1.^2;
parHparu(2,:)=(p(4,:).*(TH6.*TH7.*sin(x(2,:)) - TH3.*TH4.*sin(x(2,:) + x(3,:)).*cos(x(3,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) - x(5,:) - x(6,:) - x(4,:) + (p(6,:).*(TH0.*TH6 - TH4.^2.*sin(x(2,:) + x(3,:)).^2 + TH2.*TH6.*cos(x(2,:) + x(3,:)).^2 + TH1.*TH6.*cos(x(2,:)).^2 + TH0.*TH3.*cos(x(3,:)) + TH1.*TH3.*cos(x(2,:)).^2.*cos(x(3,:)) + 2.*TH3.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)) + 2.*TH3.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - TH4.*TH7.*sin(x(2,:) + x(3,:)).*sin(x(2,:)) + TH2.*TH3.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) - (p(5,:).*(TH0.*TH6 - TH4.^2.*sin(x(2,:) + x(3,:)).^2 + TH2.*TH6.*cos(x(2,:) + x(3,:)).^2 + TH1.*TH6.*cos(x(2,:)).^2 + 2.*TH3.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) + (2.*R2.*u(2,:))./a2.^2;
parHparu(3,:)=(p(5,:).*(TH0.*TH6 - TH4.^2.*sin(x(2,:) + x(3,:)).^2 + TH2.*TH6.*cos(x(2,:) + x(3,:)).^2 + TH1.*TH6.*cos(x(2,:)).^2 + TH0.*TH3.*cos(x(3,:)) + TH1.*TH3.*cos(x(2,:)).^2.*cos(x(3,:)) + 2.*TH3.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)) + 2.*TH3.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - TH4.*TH7.*sin(x(2,:) + x(3,:)).*sin(x(2,:)) + TH2.*TH3.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) - x(5,:) - x(6,:) - (p(4,:).*(TH4.*TH6.*sin(x(2,:) + x(3,:)) - TH4.*TH5.*sin(x(2,:) + x(3,:)) + TH6.*TH7.*sin(x(2,:)) - TH3.*TH4.*sin(x(2,:) + x(3,:)).*cos(x(3,:)) + TH3.*TH7.*cos(x(3,:)).*sin(x(2,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) - x(4,:) - (p(6,:).*(TH0.*TH5 - TH4.^2.*sin(x(2,:) + x(3,:)).^2 - TH7.^2.*sin(x(2,:)).^2 + TH2.*TH5.*cos(x(2,:) + x(3,:)).^2 + TH1.*TH5.*cos(x(2,:)).^2 + 2.*TH0.*TH3.*cos(x(3,:)) + 2.*TH1.*TH3.*cos(x(2,:)).^2.*cos(x(3,:)) + 4.*TH3.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)) + 2.*TH3.*TH5.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*sin(x(2,:)) + 2.*TH2.*TH3.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:))))./(TH0.*TH6.^2 - TH0.*TH5.*TH6 + TH2.*TH6.^2.*cos(x(2,:) + x(3,:)).^2 + TH4.^2.*TH5.*sin(x(2,:) + x(3,:)).^2 - TH4.^2.*TH6.*sin(x(2,:) + x(3,:)).^2 + TH0.*TH3.^2.*cos(x(3,:)).^2 + TH1.*TH6.^2.*cos(x(2,:)).^2 + TH6.*TH7.^2.*sin(x(2,:)).^2 - TH2.*TH5.*TH6.*cos(x(2,:) + x(3,:)).^2 - TH1.*TH5.*TH6.*cos(x(2,:)).^2 + TH2.*TH3.^2.*cos(x(2,:) + x(3,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.^3.*cos(x(2,:) + x(3,:)).*cos(x(2,:)).*cos(x(3,:)).^2 + TH1.*TH3.^2.*cos(x(2,:)).^2.*cos(x(3,:)).^2 + 2.*TH3.*TH6.^2.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH5.*TH6.*cos(x(2,:) + x(3,:)).*cos(x(2,:)) - 2.*TH3.*TH4.*TH7.*sin(x(2,:) + x(3,:)).*cos(x(3,:)).*sin(x(2,:))) + (2.*R3.*u(3,:))./a3.^2;

% Norm of the partial of the Hamiltonian with respect to the control

NormparHparu = (trapz(time,parHparu(1,:).^2+parHparu(2,:).^2+parHparu(3,:).^2))^0.5;

%Compute cost
Jold=J;
n=n+1
J=0.5*trapz(time,diag([x.',u.']*[Q,M;M.',W]*[x;u]));
JJ(iter)=J;
disp(['Iteration # ',num2str(iter),', Hu = ',num2str(NormparHparu),', J = ',num2str(J)]);

%Stopping criteria
if (NormparHparu < 1e-9) 
    break
end

if (abs(J-Jold) < 1e-9)
    break
end

% Update the control using gradient descent.
eps=[eps 10^-round(log10(NormparHparu))];
eps=max(eps)
% if (NormparHparu < 1e-3) 
%     eps=0.01;
% elseif (NormparHparu < 1e-5)
%     eps=0.001;
% end

u = u - eps * parHparu;

iter=iter+1;
end


