setupRPBC;
C=10;
sim('RPBC560.mdl');
figure
plot(t,u);legend('u1','u2','u3');
figure
plot(t,q1d,t,q1);
legend('q1d','q1')
figure
plot(t,q2d,t,q2);legend('q2d','q2')
figure
plot(t,q3d,t,q3);
legend('q3d','q3')
figure
plot(t,V_c1,t,V_c2,t,V_c3)
% legend('Link one semi-active','Link two  semi-active','Link three semi active');
% figure
% plot(t,V_c12,t,V_c13,t,V_c23)
% legend('Link one and two semi-active','Link one and three semi active','Link two and three semi-active');
% figure
% plot(t,V_c)
% title('All three links semi-active')
% figure
% plot(t,V_b12,t,V_c);
% Regen_po=(E_po(end)^2-E_po(1)^2)/(2*C)
% Regen_RaCo=E_all(end)
% Regen_RaCap=E_b1(end)-E_b1(1)
% Regen_Ra=-Regen_RaCo+Regen_RaCap