alpha=[0.189 0.219 0.202 0.223 0.226 0.240 0.255 0.220 0.202 0.258 0.245];
R=[1.6 1.6 1.6 2.1];
alpha=mean(alpha);
R=mean(R);
load('Values.mat')
load('Gear Ratios.mat')
L=fieldnames(Val);
A_Matrix=[];
B_Matrix=[];
J_Matrix=[];
%L=[L(2) L(1) L(3)];
for i=1:numel(L)
A=eval(strcat('Val.',L{i}));
n=eval(strcat('N.',L{i}));
[ a,J,B,B_e ] = Sys_Id( n,alpha,R,A(1),A(2) );
A_Matrix=[A_Matrix a];
B_Matrix=[B_Matrix B_e];
J_Matrix=[J_Matrix J];
end
function [ a,J,B,B_e ] = Sys_Id( n,alpha,R,A_J,B_J )
a=n*alpha/R;
J=a/A_J;
B_e=a^2/R;
B=J*(B_J)-B_e;

end


