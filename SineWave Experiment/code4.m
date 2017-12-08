len1='Actual';
len2='Control Signal';
ylab='Displacement (rad)';
figure('units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
plot(TEST.t,TEST.q1,'-r',TEST.t,TEST.q1d,'-b')
legend(len1, len2)
xlabel('Time (Sec)')
ylabel(ylab)
suptitle('')
title('Motor 1')
subplot(3,1,2)
plot(TEST.t,TEST.q2,'-r',TEST.t,TEST.q2d,'-b')
legend(len1, len2)
xlabel('Time (Sec)')
ylabel(ylab)
suptitle('')
title('Motor 2')
subplot(3,1,3)
plot(TEST.t,TEST.q3,'-r',TEST.t,TEST.q3d,'-b')
legend(len1, len2)
xlabel('Time (Sec)')
ylabel(ylab)
suptitle('')
title('Motor 3')
print -depsc Voltages.eps