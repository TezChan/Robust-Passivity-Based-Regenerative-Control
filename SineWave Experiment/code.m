clear all
Files={'with_regen2'};%Change this into with_regen2 for the other case
Folders={'ALL0002'};%change this into ALL0002 for the other case 
Tests={'Test'};
pul=1;
for i=1:length(Files)
eval(strcat(upper(Tests{i}),'=run(Files{i},pul);'));
end
for i=1:length(Folders)
eval(strcat(Tests{i},'=loadfiles(','Folders{i},pul',');'));
end
% Test.m1=-[Test.m1 ;0];
% Test.m2=-[Test.m2 ;0];
% Test.m3=-[Test.m3 ;0];
% TEST.C_Current_2=-2*TEST.C_Current_2;
power=[Test.m1(1:end-1)'.*TEST.M_Current_1;Test.m2(1:end-1)'.*TEST.M_Current_2;Test.m3(1:end-1)'.*TEST.M_Current_3];
power_cap=[TEST.CapVolt.*TEST.C_Current_1;TEST.CapVolt.*TEST.C_Current_2;TEST.CapVolt.*TEST.C_Current_3];
I=find(power<=0);
Regen_Motors=trapz(power(I))*(TEST.t(2)-TEST.t(1));
p_c=sum(power_cap);
I2=find(p_c<=0);
Regen_Cap=trapz(p_c(I2))*(TEST.t(2)-TEST.t(1));
TEST.Ct=TEST.C_Current_1+TEST.C_Current_2+TEST.C_Current_3;
Test.m1=Test.m1(1:end-1);
Test.m2=Test.m2(1:end-1);
Test.m3=Test.m3(1:end-1);
k=4;
figure()
subplot(3,k,1)
plot(TEST.t,TEST.C_Current_1,'-r',TEST.t,TEST.M_Current_1,'-b')
title('Motor 1 Currents')
legend('Capactior current','Motor current')
subplot(3,k,2)
plot(TEST.t,TEST.V1,'-r',TEST.t,Test.m1,'-b')
title('Motor 1 Voltages')
legend('Signal','Motor Voltage')
subplot(3,k,1+k)
plot(TEST.t,TEST.C_Current_2,'-r',TEST.t,TEST.M_Current_2,'-b')
title('Motor 2 Currents')
legend('Capactior current','Motor current')
subplot(3,k,1+2*k)
plot(TEST.t,TEST.C_Current_3,'-r',TEST.t,TEST.M_Current_3,'-b')
title('Motor 3 Currents')
legend('Capactior current','Motor current')
subplot(3,k,2+k)
plot(TEST.t,-TEST.V2,'-r',TEST.t,Test.m2,'-b')
title('Motor 2 signal and Motor Voltage')
legend('Signal','Motor Voltage')
subplot(3,k,2+2*k)
plot(TEST.t,TEST.V3,'-r',TEST.t,Test.m3,'-b')
title('Motor 3 signal and Motor Voltage')
legend('Signal','Motor Voltage')
titles={'Motor 1 ','Motor 2 ','Motor 3 '};
for i=1:3
subplot(3,k,3+k*(i-1))
plot(TEST.t,power_cap(i,:),TEST.t,power(i,:))
title(strcat(titles{i},' power'));
legend('Power sent into Capacitor','Motor Power');
subplot(3,k,4+k*(i-1))
eval(['plot(TEST.t,TEST.q',int2str(i),'d,TEST.t,TEST.q',int2str(i),')'])
title(strcat(titles{i},' displacement'));
legend('qd','q')
end
figure()
subplot(2,1,1)
plot(TEST.t,TEST.CapVolt)
title('Capacitor Voltage VS Time')
subplot(2,1,2)
plot(TEST.t,sum(power_cap))
title('Capacitor power VS Time')
%Enegries
enegry1=trapz(TEST.t,Test.m1.*TEST.M_Current_1');
enegry2=trapz(TEST.t,Test.m2.*TEST.M_Current_2');
enegry3=trapz(TEST.t,Test.m3.*TEST.M_Current_3');
enegry=trapz(TEST.t,TEST.Ct.*TEST.CapVolt);
eff=(enegry1+enegry2+enegry3)/enegry;
clearvars Folders Files a b i titles Tests
function Data=run(file,pul)
load(file)
Data.t=eval(strcat(file,'.X.Data'));
for i=1:17
     a=strcat(file,'.Y(',int2str(i),')');
     b=eval(strcat(a,'.Name(11:end-1)'));
     b=regexprep(b,'[^\w'']','');
     eval(strcat('Data.',b,'=',a,'.Data;'));
end
I=find(diff(Data.pulse)>=4);
fields=fieldnames(Data);
for i=1:numel(fields)
eval(strcat('Data.',fields{i},'=Data.',fields{i},'(',int2str(I(end-pul)),':',int2str(I(end)),');'));
end
Data.t=Data.t(:)-Data.t(1);
for i=1:numel(fields)
eval(strcat('Data.',fields{i},'=Data.',fields{i},'(1:4:end);'));
end
end
function Data=loadfiles(Folder,pul)
files=dir(strcat(Folder,'/*.CSV'));
files={files.name};
C={'m1','m2','m3','pulse'};
for i=1:4
a=[Folder,'/',files{i}];
eval(strcat('[Data.t Data.',C{i},']=importfile(a);'));
end
I=find(((diff(Data.pulse))>4) &((diff(Data.pulse))<6));
Data.m1=Data.m1(I(end-pul):I(end));
Data.m2=Data.m2(I(end-pul):I(end));
Data.m3=Data.m3(I(end-pul):I(end));
Data.pulse=Data.pulse(I(end-pul):I(end));
Data.t=Data.t(I(end-pul):I(end));
Data.t=Data.t(:)-Data.t(1);
end
function [VarName4,VarName5]=importfile(filename)
%% Initialize variables.
delimiter = ',';

%% Format for each line of text:
%   column4: double (%f)
%	column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*q%*q%*q%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
VarName4 = dataArray{:, 1};
VarName5 = dataArray{:, 2};


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
end