script;
script2;
for k=1:numel(I)
subplot(3,1)
hold on
plot(1:100,Q1(I(k),:))
ylabel('q1')
subplot(3,2)
hold on
plot(1:100,Q2(I(k),:))
ylabel('q2')
subplot(3,3)
hold on
plot(1:100,Q3(I(k),:))
ylabel('q3')
end