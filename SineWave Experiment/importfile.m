unction [VarName4,VarName5]=importfile(filename)
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