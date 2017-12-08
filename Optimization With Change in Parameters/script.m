clear all
[q1opt,q2opt,q3opt]=DCopt(0.2254,1.725,0);
Q1=[];Q2=[];Q3=[];
zee=[];
trails=[];
n=0;
c=0;
for R=1.725:0.05*1.725:1.725*2
    n=n+1;
    c=0;
    for rate=-0.5:0.05:1.5
        zee=[zee [R;rate]];
        c=c+1;
    end
end
n=1;
while ~isempty(zee)
    [Q1(n,:),Q2(n,:),Q3(n,:)]=DCopt(zee(1,1),0.2254,zee(2,1));
    trails=[trails zee(:,1)];
    zee(:,1)=[];
    n=n+1;
end