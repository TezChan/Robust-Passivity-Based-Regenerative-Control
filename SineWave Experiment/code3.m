close all
figure
plot(Power1.fq.t,Power1.fq.M,Power1.fq.t,Power1.fq.C)
title('Power first quadrent')
Enegry1.fq.M=trapz(Power1.fq.t,Power1.fq.M);
Enegry1.fq.C=trapz(Power1.fq.t,Power1.fq.C);
figure
plot(Power1.sq.t,Power1.sq.M,Power1.sq.t,Power1.sq.C)
Enegry1.sq.M=trapz(Power1.sq.t,Power1.sq.M);
Enegry1.sq.C=trapz(Power1.sq.t,Power1.sq.C);
figure
plot(Power1.tq.t,Power1.tq.M,Power1.tq.t,Power1.tq.C)
Enegry1.tq.M=trapz(Power1.tq.t,Power1.tq.M);
Enegry1.tq.C=trapz(Power1.tq.t,Power1.tq.C);
figure
plot(Power1.ftq.t,Power1.ftq.M,Power1.ftq.t,Power1.ftq.C)
Enegry1.ftq.M=trapz(Power1.ftq.t,Power1.ftq.M);
Enegry1.ftq.C=trapz(Power1.ftq.t,Power1.ftq.C);
close all
figure
plot(Power2.fq.t,Power2.fq.M,Power2.fq.t,Power2.fq.C)
title('Power first quadrent')
Enegry2.fq.M=trapz(Power2.fq.t,Power2.fq.M);
Enegry2.fq.C=trapz(Power2.fq.t,Power2.fq.C);
figure
plot(Power2.sq.t,Power2.sq.M,Power2.sq.t,Power2.sq.C)
Enegry2.sq.M=trapz(Power2.sq.t,Power2.sq.M);
Enegry2.sq.C=trapz(Power2.sq.t,Power2.sq.C);
figure
plot(Power2.tq.t,Power2.tq.M,Power2.tq.t,Power2.tq.C)
Enegry2.tq.M=trapz(Power2.tq.t,Power2.tq.M);
Enegry2.tq.C=trapz(Power2.tq.t,Power2.tq.C);
figure
plot(Power2.ftq.t,Power2.ftq.M,Power2.ftq.t,Power2.ftq.C)
Enegry2.ftq.M=trapz(Power2.ftq.t,Power2.ftq.M);
Enegry2.ftq.C=trapz(Power2.ftq.t,Power2.ftq.C);