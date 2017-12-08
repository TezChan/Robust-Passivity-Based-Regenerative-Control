%controlregr
%Forms and verifies the control regressor  Yav

dynamicsTH
syms a_1 a_2 a_3 v_1 v_2 v_3

%Replace every instance of joint acceleration with virtual 
%acceleration

%Virtual velocity substitution is more nuanced
%Only qdot=[q1dot;q2dot;q3dot] appearing in C*qdot is substituted,
%not the qdots "inside" C

aux=D*[a_1;a_2;a_3]+C*[v_1;v_2;v_3]+gTH;

%Use partial differentiation to extract control regressor
for i=1:3, %number of rows (links)
    for j=1:10, %number of parameters
        Yav(i,j)=diff(aux(i),TH(j));
    end
end

%Verify
aux1=D_q*[a_1;a_2;a_3]+C_q*[v_1;v_2;v_3]+gq;

simplify(eval(aux-aux1))