function [displacement1,displacement2,displacement3]= DCopt(alpha,R,r)

close all
% Direct collocation solution to Puma optimal control problem

global Problem Puma_model 

Problem.N = 100; % number of collocation points

Problem.t_ini=0; % initial time

Problem.t_f=2; % final time

Problem.h = Problem.t_f/(Problem.N-1);    % time step for direct collocation

Problem.t = Problem.t_ini:Problem.h:Problem.t_f;    % time points for direct collocation

Problem.discretization='BE';
%%system parameters
Puma_model.R1=R;
Puma_model.R2=R;
Puma_model.R3=R;
Puma_model.a1=alpha*62.61;
Puma_model.a2=alpha*107.36;
Puma_model.a3=alpha*53.69;
rate=(1+r);
Puma_model.TH0=2.8123*rate;
Puma_model.TH1=2.2623*rate;
Puma_model.TH2=-0.0066*rate;
Puma_model.TH3=0.3453*rate;
Puma_model.TH4=-0.085*rate;
Puma_model.TH5=8.8038*rate;
Puma_model.TH6=1.0711*rate;
Puma_model.TH7=-0.9748*rate;
Puma_model.TH8=6.8067*rate;
Puma_model.TH9=55.35*rate;

% collocation grid and unknowns
    Puma_model.Nstates=6;               %number of states
    Puma_model.Nvarpernode = 9;			% number of unknowns per node: 6 states and 3 controls q1,q2,q3,q1dot,q2dot,q3dot,u1,u2,u3
	Puma_model.Nconpernode = 6;         % number of constraint equations per node(number of state equations)
    Puma_model.Nvar = Problem.N * Puma_model.Nvarpernode;   % number of unknowns 
	Puma_model.Ncon = (Problem.N-1)* Puma_model.Nconpernode + Puma_model.Nstates*2 ;		% number of state equation constraints + initial and final cconditions
    
    
% variable inexes
    Problem.iq1= 1:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;     % index of q1 within X
    Problem.iq2= 2:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;     % index of q2 within X
    Problem.iq3= 3:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;     % index of q3 within X
    Problem.iq1dot= 4:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;  % index of q1dot within X
    Problem.iq2dot= 5:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;  % index of q2dot within X
    Problem.iq3dot= 6:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;  % index of q3dot within X
    Problem.iu1 = 7:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;    % index of u1 within X
    Problem.iu2 = 8:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;    % index of u2 within X
    Problem.iu3 = 9:Puma_model.Nvarpernode:(Problem.N)*Puma_model.Nvarpernode;    % index of u3 within X
    
% Initial and final conditions x=[q1;q2;q3;q1dot;q2dot;q3dot]

Problem.x_ini=[0;-pi/3;pi/3;0;0;0];
Problem.x_f=[pi;pi/3;-pi/3;0;0;0];




% Determine Jacobian structure and number of non-zeros
Xr = rand(Puma_model.Nvar,1);
Puma_model.J_nnz = 1;
J = jac_confun(Xr);
Puma_model.structure_jac_confun = double(J ~= 0);
Puma_model.J_nnz = nnz(J);
    
     % Check the derivations
     %checkderiv
     
    % Bounds on the variables (q1,q2,q3,q1dot,q2dot,q3dot,u1,u2,u3)   
    Vcap=25;
%     xlb = [-inf ; -inf ; -inf ; -inf ; -inf ; -inf ; -Puma_model.a1/Puma_model.R1*Vcap; -Puma_model.a2/Puma_model.R2*Vcap;-Puma_model.a3/Puma_model.R3*Vcap];
%     xub = [ inf ;  inf ; inf;  inf ;  inf ;  inf ; Puma_model.a1/Puma_model.R1*Vcap; Puma_model.a2/Puma_model.R2*Vcap;Puma_model.a3/Puma_model.R3*Vcap];

    xlb = [-inf ; -inf ; -inf ; -inf ; -inf ; -inf ; -inf; -inf;-inf];
    xub = [ inf ;  inf ; inf;  inf ;  inf ;  inf ; inf; inf;inf];
    
    
    %initial guess
    a=exist('X0');
    if a==0
%     X0 = zeros(Puma_model.Nvar,1);       %Included the trajectories all the 9 variable
     X0 = rand(Puma_model.Nvar,1);
     X0(Problem.iu1)=400*X0(Problem.iu1);
     X0(Problem.iu2)=800*X0(Problem.iu2);
     X0(Problem.iu3)=400*X0(Problem.iu3);
    end
    
    
    % solve the NLP with IPOPT
    funcs.objective = @objfun;
    funcs.gradient  = @grad_objfun;
    funcs.constraints = @confun;
    funcs.jacobian    = @jac_confun;
    funcs.jacobianstructure = @structure_jac_confun;

    options.lb = repmat(xlb, Problem.N, 1);
    options.ub = repmat(xub, Problem.N, 1);
    options.cl = zeros(Puma_model.Ncon,1);
    options.cu = zeros(Puma_model.Ncon,1);	
    options.ipopt.max_iter = 50000;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.tol = 1e-4;
% 	options.ipopt.print_level = 0;       
    [Xres, info] = ipopt(X0,funcs,options);
    
	if info.status ~= 0
		error('IPOPT did not solve')  
    end
    
    disp(['Energy Regenerated: ',num2str(objfun(Xres))])
    
% figure(1)
% subplot(2,3,1)
% hold on
% plot(Problem.t,Xres(Problem.iq1))
% ylabel('q1')
% subplot(2,3,2)
% hold on
% plot(Problem.t,Xres(Problem.iq2))
% ylabel('q2')
% subplot(2,3,3)
% hold on
% plot(Problem.t,Xres(Problem.iq3))
% ylabel('q3')
% subplot(2,3,4)
% hold on
% plot(Problem.t,Xres(Problem.iu1))
% ylabel('u1')
% subplot(2,3,5)
% hold on
% plot(Problem.t,Xres(Problem.iu2))
% ylabel('u2')
% subplot(2,3,6)
% hold on
% plot(Problem.t,Xres(Problem.iu3))
% ylabel('u3')
 displacement1=Xres(Problem.iq1);
 displacement2=Xres(Problem.iq2);
 displacement3=Xres(Problem.iq3);
% p1=-Xres(Problem.iq1dot).*Xres(Problem.iu1)+Puma_model.R1/(Puma_model.a1^2)*Xres(Problem.iu1).^2;
% p2=-Xres(Problem.iq2dot).*Xres(Problem.iu2)+Puma_model.R2/(Puma_model.a2^2)*Xres(Problem.iu2).^2;
% p3=-Xres(Problem.iq3dot).*Xres(Problem.iu3)+Puma_model.R3/(Puma_model.a3^2)*Xres(Problem.iu3).^2;
% p=p1+p2+p3;
% 
% figure
% hold on
% plot(Problem.t,p1,Problem.t,p2,Problem.t,p3,Problem.t,p)
% legend('p1','p2','p3','p')
% xlabel('time (s)')
% ylabel('Power (w)')
% 
function [f] = objfun(X)
    
% energy stored in common capacitor
    f=Problem.h*sum(-X(Problem.iq1dot).*X(Problem.iu1)+Puma_model.R1/(Puma_model.a1^2)*X(Problem.iu1).^2+...
        -X(Problem.iq2dot).*X(Problem.iu2)+Puma_model.R2/(Puma_model.a2^2)*X(Problem.iu2).^2+...
    -X(Problem.iq3dot).*X(Problem.iu3)+Puma_model.R3/(Puma_model.a3^2)*X(Problem.iu3).^2);

    % a short pause to make sure that IPOPT screen output appears
    % continuously (if print_level was set accordingly)
    pause(1e-6);	
end



function g = grad_objfun(X)
   
    % initialize gradient
    g = zeros(size(X));
	
    g(Problem.iq1dot)=-Problem.h*X(Problem.iu1);
    g(Problem.iq2dot)=-Problem.h*X(Problem.iu2);
    g(Problem.iq3dot)=-Problem.h*X(Problem.iu3);
    g(Problem.iu1)=Problem.h*(-X(Problem.iq1dot)+2*Puma_model.R1/(Puma_model.a1^2)*X(Problem.iu1));
    g(Problem.iu2)=Problem.h*(-X(Problem.iq2dot)+2*Puma_model.R2/(Puma_model.a2^2)*X(Problem.iu2));
    g(Problem.iu3)=Problem.h*(-X(Problem.iq3dot)+2*Puma_model.R3/(Puma_model.a3^2)*X(Problem.iu3));

end

function f = dynfun(z,zdot)
    % the five dynamics equations of the system, in the form
    %   f(x, dx/dt) = 0
    % where x contains the state variables and the controls
    
q1=z(1);
q2=z(2);
q3=z(3);
q1dot=z(4);
q2dot=z(5);
q3dot=z(6);

taud=z(7:9);


D(1,1)=Puma_model.TH2*cos(q2 + q3)^2 + 2*Puma_model.TH3*cos(q2 + q3)*cos(q2) + Puma_model.TH1*cos(q2)^2 + Puma_model.TH0;
D(1,2)=Puma_model.TH4*sin(q2 + q3) + Puma_model.TH7*sin(q2);
D(1,3)=Puma_model.TH4*sin(q2 + q3);
D(2,1)=D(1,2);
D(2,2)=Puma_model.TH5 + 2*Puma_model.TH3*cos(q3);
D(2,3)=Puma_model.TH6 + Puma_model.TH3*cos(q3);
D(3,1)=D(1,3);
D(3,2)=D(2,3);
D(3,3)=Puma_model.TH6;

C(1,1)=-q2dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*cos(q2 + q3)*sin(q2) + Puma_model.TH3*sin(q2 + q3)*cos(q2) + Puma_model.TH1*cos(q2)*sin(q2)) - q3dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*sin(q2 + q3)*cos(q2));
C(1,2)=q2dot*(Puma_model.TH4*cos(q2 + q3) + Puma_model.TH7*cos(q2)) - q1dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*cos(q2 + q3)*sin(q2) + Puma_model.TH3*sin(q2 + q3)*cos(q2) + Puma_model.TH1*cos(q2)*sin(q2)) + Puma_model.TH4*q3dot*cos(q2 + q3);
C(1,3)=Puma_model.TH4*q2dot*cos(q2 + q3) - q1dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*sin(q2 + q3)*cos(q2)) + Puma_model.TH4*q3dot*cos(q2 + q3);
C(2,1)=q1dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*cos(q2 + q3)*sin(q2) + Puma_model.TH3*sin(q2 + q3)*cos(q2) + Puma_model.TH1*cos(q2)*sin(q2));
C(2,2)=-Puma_model.TH3*q3dot*sin(q3);
C(2,3)=-Puma_model.TH3*q2dot*sin(q3) -Puma_model.TH3*q3dot*sin(q3);
C(3,1)=q1dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*sin(q2 + q3)*cos(q2));
C(3,2)=Puma_model.TH3*q2dot*sin(q3);
C(3,3)=0;

gg=[0;-Puma_model.TH8*cos(q2 + q3)-Puma_model.TH9*cos(q2);-Puma_model.TH8*cos(q2 + q3)];

Bemf=diag([Puma_model.a1^2/Puma_model.R1;Puma_model.a2^2/Puma_model.R2;Puma_model.a3^2/Puma_model.R3]);

f=[zdot(1:3)-z(4:6);D*zdot(4:6)-(taud-C*z(4:6)-Bemf*z(4:6)-gg)];

% Scaling of the elements 
S = diag([1 1 1 1 1 1]);
f = S*f;
    
end

function [dfdx, dfdxdot]=jac_dynfun(z,zdot)
    
q1=z(1);
q2=z(2);
q3=z(3);
q1dot=z(4);
q2dot=z(5);
q3dot=z(6);

taud1=z(7);
taud2=z(8);
taud3=z(9);

q1ddot=zdot(4);
q2ddot=zdot(5);
q3ddot=zdot(6);

% initialize the Jacobian matrices df/dx, df/dxdot
dfdx = spalloc(Puma_model.Nconpernode ,Puma_model.Nvarpernode,21);
dfdxdot = spalloc(Puma_model.Nconpernode ,Puma_model.Nvarpernode,12);
    
% equation 1:3 zdot(1:3)-z(4:6)= 0:
dfdx(1:3,4:6) = -eye(3);
dfdxdot(1:3,1:3) = eye(3);

 
% equation 4 
dfdx(4,2) = Puma_model.TH4*q2ddot*cos(q2 + q3) - Puma_model.TH4*q2dot^2*sin(q2 + q3) - Puma_model.TH4*q3dot^2*sin(q2 + q3) ...
- Puma_model.TH1*q1ddot*sin(2*q2) - Puma_model.TH7*q2dot^2*sin(q2) - Puma_model.TH2*q1ddot*sin(2*q2 + 2*q3) - 2*Puma_model.TH3*q1ddot*sin(2*q2 + q3)...
+ Puma_model.TH4*q3ddot*cos(q2 + q3) + Puma_model.TH7*q2ddot*cos(q2) - 2*Puma_model.TH4*q2dot*q3dot*sin(q2 + q3) - 4*Puma_model.TH3*q1dot*q2dot*cos(2*q2 + q3)...
- 2*Puma_model.TH3*q1dot*q3dot*cos(2*q2 + q3) - 2*Puma_model.TH1*q1dot*q2dot*cos(2*q2) - 2*Puma_model.TH2*q1dot*q2dot*cos(2*q2 + 2*q3) ...
- 2*Puma_model.TH2*q1dot*q3dot*cos(2*q2 + 2*q3);

dfdx(4,3)= Puma_model.TH4*q2ddot*cos(q2 + q3) - Puma_model.TH3*q1ddot*sin(2*q2 + q3) - Puma_model.TH4*q2dot^2*sin(q2 + q3) - Puma_model.TH4*q3dot^2*sin(q2 + q3)...
- Puma_model.TH2*q1ddot*sin(2*q2 + 2*q3) - Puma_model.TH3*q1ddot*sin(q3) + Puma_model.TH4*q3ddot*cos(q2 + q3) - 2*Puma_model.TH4*q2dot*q3dot*sin(q2 + q3) ...
- Puma_model.TH3*q1dot*q3dot*cos(q3) - 2*Puma_model.TH3*q1dot*q2dot*cos(2*q2 + q3) - Puma_model.TH3*q1dot*q3dot*cos(2*q2 + q3) - 2*Puma_model.TH2*q1dot*q2dot*cos(2*q2 + 2*q3)...
- 2*Puma_model.TH2*q1dot*q3dot*cos(2*q2 + 2*q3); 

dfdx(4,4)=Puma_model.a1^2/Puma_model.R1 - Puma_model.TH3*q3dot*sin(q3) - 2*Puma_model.TH3*q2dot*sin(2*q2 + q3) - Puma_model.TH3*q3dot*sin(2*q2 + q3)...
    - Puma_model.TH1*q2dot*sin(2*q2) - Puma_model.TH2*q2dot*sin(2*q2 + 2*q3) - Puma_model.TH2*q3dot*sin(2*q2 + 2*q3);

dfdx(4,5)=2*Puma_model.TH4*q2dot*cos(q2 + q3) - Puma_model.TH1*q1dot*sin(2*q2) - Puma_model.TH2*q1dot*sin(2*q2 + 2*q3) - 2*Puma_model.TH3*q1dot*sin(2*q2 + q3)...
    + 2*Puma_model.TH4*q3dot*cos(q2 + q3) + 2*Puma_model.TH7*q2dot*cos(q2);

dfdx(4,6)=2*Puma_model.TH4*q2dot*cos(q2 + q3) - 2*q1dot*(Puma_model.TH2*cos(q2 + q3)*sin(q2 + q3) + Puma_model.TH3*sin(q2 + q3)*cos(q2))...
    + 2*Puma_model.TH4*q3dot*cos(q2 + q3);

dfdx(4,7)=-1;


dfdxdot(4,4) = Puma_model.TH0 + Puma_model.TH2*cos(q2 + q3)^2 + Puma_model.TH1*cos(q2)^2 + 2*Puma_model.TH3*cos(q2 + q3)*cos(q2);

dfdxdot(4,5) =Puma_model.TH4*sin(q2 + q3) + Puma_model.TH7*sin(q2);

dfdxdot(4,6) =Puma_model.TH4*sin(q2 + q3);


% equation 5

dfdx(5,2)= Puma_model.TH8*sin(q2 + q3) + Puma_model.TH9*sin(q2) + 2*Puma_model.TH3*q1dot^2*cos(2*q2 + q3) + Puma_model.TH1*q1dot^2*cos(2*q2) ...
    + Puma_model.TH4*q1ddot*cos(q2 + q3) + Puma_model.TH2*q1dot^2*cos(2*q2 + 2*q3) + Puma_model.TH7*q1ddot*cos(q2);

dfdx(5,3)= Puma_model.TH8*sin(q2 + q3) - 2*Puma_model.TH3*q2ddot*sin(q3) - Puma_model.TH3*q3ddot*sin(q3) - Puma_model.TH3*q3dot^2*cos(q3)...
    + Puma_model.TH3*q1dot^2*cos(2*q2 + q3) + Puma_model.TH4*q1ddot*cos(q2 + q3) + Puma_model.TH2*q1dot^2*cos(2*q2 + 2*q3) -...
    2*Puma_model.TH3*q2dot*q3dot*cos(q3);

dfdx(5,4)=q1dot*(Puma_model.TH2*sin(2*q2 + 2*q3) + 2*Puma_model.TH3*sin(2*q2 + q3) + Puma_model.TH1*sin(2*q2));

dfdx(5,5)=Puma_model.a2^2/Puma_model.R2 - 2*Puma_model.TH3*q3dot*sin(q3);

dfdx(5,6)=-2*Puma_model.TH3*sin(q3)*(q2dot + q3dot);

dfdx(5,8)=-1;

dfdxdot(5,4) =Puma_model.TH4*sin(q2 + q3) + Puma_model.TH7*sin(q2);

dfdxdot(5,5)=Puma_model.TH5 + 2*Puma_model.TH3*cos(q3);

dfdxdot(5,6)=Puma_model.TH6 + Puma_model.TH3*cos(q3);

% equation 6

dfdx(6,2)=Puma_model.TH8*sin(q2 + q3) + Puma_model.TH3*q1dot^2*cos(2*q2 + q3) + Puma_model.TH4*q1ddot*cos(q2 + q3) + Puma_model.TH2*q1dot^2*cos(2*q2 + 2*q3);

dfdx(6,3)=Puma_model.TH8*sin(q2 + q3) - Puma_model.TH3*q2ddot*sin(q3) + Puma_model.TH3*q2dot^2*cos(q3) + Puma_model.TH2*q1dot^2*cos(q2 + q3)^2- ...
 Puma_model.TH2*q1dot^2*sin(q2 + q3)^2 + Puma_model.TH4*q1ddot*cos(q2 + q3) + Puma_model.TH3*q1dot^2*cos(q2 + q3)*cos(q2);

dfdx(6,4)=2*q1dot*sin(q2 + q3)*(Puma_model.TH2*cos(q2 + q3) + Puma_model.TH3*cos(q2));

dfdx(6,5)=2*Puma_model.TH3*q2dot*sin(q3);

dfdx(6,6)=Puma_model.a3^2/Puma_model.R3;

dfdx(6,9)=-1;

dfdxdot(6,4) =Puma_model.TH4*sin(q2 + q3);

dfdxdot(6,5) =Puma_model.TH6 + Puma_model.TH3*cos(q3);

dfdxdot(6,6) =Puma_model.TH6;
     
   
        
% Scaling of the elements 
S = diag([1 1 1 1 1 1]);
dfdx = S*dfdx;
dfdxdot = S*dfdxdot;
     
end

function c = confun(X)
    
    %function containing the constraints of the problem
   
	c = zeros(Puma_model.Ncon,1);     % initialize the constraints
        
    ix1 = 1:Puma_model.Nvarpernode;   % index for the variables of node 1
    ic  = 1:Puma_model.Nconpernode;   % index for constraints from node 1
    
    for i = 1:Problem.N-1
        % extract variables from successive nodes
        x1 = X(ix1);
        x2 = X(ix1 + Puma_model.Nvarpernode);  
        
      
    	% dynamics constraints are calculated by dynfun function
        
        % use Backward Euler formula or Midpoint Euler formula as dynamics constraint
        if strcmp(Problem.discretization, 'BE')
            c(ic) = dynfun(x2 , (x2-x1)/Problem.h);
        else
            c(ic) = dynfun((x1+x2)/2 , (x2-x1)/Problem.h);
        end
        
        ix1 = ix1 + Puma_model.Nvarpernode;
        ic  = ic  + Puma_model.Nconpernode;
    end

% initial and final conditions constraints
c(Puma_model.Nconpernode*(Problem.N-1)+(1:2*Puma_model.Nstates)) = ...
    [X((1:Puma_model.Nstates))-Problem.x_ini;X((Problem.N-1)*Puma_model.Nvarpernode + (1:Puma_model.Nstates))-Problem.x_f];


end

function J = jac_confun(X)

% initialize the sparse Jacobian matrix
J= spalloc(Puma_model.Ncon,Puma_model.Nvar, Puma_model.J_nnz);		
    
 	ix1 = 1:Puma_model.Nvarpernode;   % index for the variables of node 1
    ic  = 1:Puma_model.Nconpernode;   % index for constraints from node 1
    
    for i=1:Problem.N-1
        
		% extract variables from two successive nodes
		x1 = X(ix1);
		x2 = X(ix1 + Puma_model.Nvarpernode);
        
        if strcmp(Problem.discretization, 'BE')
            [dfdx, dfdxdot] = jac_dynfun(x2, (x2-x1)/Problem.h);
            J(ic,ix1) = -dfdxdot/Problem.h;
            J(ic,ix1 + Puma_model.Nvarpernode) = dfdx + dfdxdot/Problem.h;
        else
            [dfdx, dfdxdot] = jac_dynfun((x1+x2)/2, (x2-x1)/Problem.h);
            J(ic,ix1) = dfdx/2 - dfdxdot/h;
            J(ic,ix1 + Puma_model.Nvarpernode) = dfdx/2 + dfdxdot/Problem.h;
        end
        
		%  advance ix1 and ic to next node
		ix1 = ix1 + Puma_model.Nvarpernode;
		ic  = ic  + Puma_model.Nconpernode;
        
    end
    % Jacobian elements from the intitla and final conditions constraints
    J(Puma_model.Nconpernode*(Problem.N-1)+(1:Puma_model.Nstates),(1:Puma_model.Nstates))=eye(Puma_model.Nstates);
    J(Puma_model.Nconpernode*(Problem.N-1)+(Puma_model.Nstates+1:2*Puma_model.Nstates),Puma_model.Nvarpernode*(Problem.N-1)+(1:Puma_model.Nstates))=eye(Puma_model.Nstates);
    

end

function J = structure_jac_confun(X)
    
J = Puma_model.structure_jac_confun;

end

 function checkderiv
	% using finite differences to check that the code in grad_objfun and jac_confun is correct
    N_Check=10;
    NX=N_Check * Puma_model.Nvarpernode;
	hh= 1e-6;
    X_Check = randn(NX,1);
	f        = objfun(X_Check);
	grad     = grad_objfun(X_Check);
	c        = confun(X_Check);
    Ncon     = size(c,1);
	cjac     = jac_confun(X_Check);
	cjac_num = zeros(Ncon,NX);
	grad_num = zeros(NX,1);
	
    for i=1:NX
        fprintf('checking derivatives for unknown %4d of %4d\n',i,NX);
        Xisave        = X_Check(i);
        X_Check(i)    = X_Check(i) + hh;
        cjac_num(:,i) = (confun(X_Check) - c)/hh;
        grad_num(i)   = (objfun(X_Check) - f)/hh;
        X_Check(i)          = Xisave;
    end
	
	% report maximal differences between analytical derivatives and numerical results
	%fprintf('Max. error in constraint jacobian: ');
	%matcompare(cjac, cjac_num);
	%fprintf('Max. error in objective gradient: ');
	%matcompare(grad, grad_num);
	%disp('Type dbcont to continue');

    %keyboard	
    end
 %====================================================================
    function matcompare(a,b)
	% compares two matrices and prints element that has greatest difference
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('%9.6f at %d %d (%9.6f vs. %9.6f)\n', full(maxerr), irow, icol, full(a(irow,icol)), full(b(irow,icol)));
    end


end

