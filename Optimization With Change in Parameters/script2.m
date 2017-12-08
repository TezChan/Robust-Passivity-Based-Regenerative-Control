clear all
load('861.mat')
Q1=rem(Q1,pi);Q2=rem(Q2,pi);Q3=rem(Q3,pi);
figure(1)
for i =1:861
D(i,:)=[sum((Q1(i,:)-q1opt').^2)/100;sum(((Q2(i,:))-q2opt').^2)/100;sum((Q3(i,:)-q3opt').^2)/100];%Meansquare change compared to values
subplot(2,3,1)
hold on
plot(1:100,Q1(i,:))
ylabel('q1')
subplot(2,3,2)
hold on
plot(1:100,Q2(i,:))
ylabel('q2')
subplot(2,3,3)
hold on
plot(1:100,Q3(i,:))
ylabel('q3')
end
R=reshape(trails(1,:),c,861/c);
rate=reshape(trails(2,:),c,861/c);
% q_1=reshape(D(:,1),c,861/c);
% q_2=reshape(D(:,2),c,861/c);
% q_3=reshape(D(:,3),c,861/c);
% 
% subplot(2,3,4)
% mesh(R,rate,q_1)
% xlabel('R')
% ylabel('Change in parameters (ratio)')
% zlabel('Change in q1')
% 
% subplot(2,3,5)
% mesh(R,rate,q_2)
% xlabel('R')
% ylabel('Change in parameters (ratio)')
% zlabel('Change in q2')
% 
% subplot(2,3,6)
% mesh(R,rate,q_3)
% xlabel('R')
% ylabel('Change in parameters (ratio)')
% zlabel('Change in q3')
