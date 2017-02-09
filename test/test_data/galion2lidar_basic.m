function data = galion2lidar_basic(allow_stacking, force_monospace_time, output_file, varargin) %#codegen
%GALION2LIDAR_BASIC 

data = read_galion_day_file(varargin{1});
for iFile = 2:numel(varargin)
    fprintf('Processing file %i of %i, %s\n', iFile, numel(varargin), varargin{iFile})
    data_chunk = read_galion_day_file(varargin{iFile});
    
    % Append the data chunk. Horrid memory acces pattern
    data.Datetime       = vertcat(data.Datetime, data_chunk.Datetime);
    data.Heightrange    = vertcat(data.Heightrange, data_chunk.Heightrange);
    data.Hsp            = vertcat(data.Hsp, data_chunk.Hsp);
    data.Vsp            = vertcat(data.Vsp, data_chunk.Vsp);
    data.Wdir           = vertcat(data.Wdir, data_chunk.Wdir);
    data.Turb           = vertcat(data.Turb, data_chunk.Turb);
    data.Minintensity   = vertcat(data.Minintensity, data_chunk.Minintensity);
    data.Meanintensity  = vertcat(data.Meanintensity, data_chunk.Meanintensity);
    data.NE             = vertcat(data.NE, data_chunk.NE);
    data.ES             = vertcat(data.ES, data_chunk.ES);
    data.SW             = vertcat(data.SW, data_chunk.SW);
    data.WN             = vertcat(data.WN, data_chunk.WN);
    data.Vector         = vertcat(data.Vector, data_chunk.Vector);

end

end

function data = read_galion_day_file(filename)

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of tab-delimited data as a string (for date) then floating point
% values, ignoring the header line
delimiter = '\t';
startRow = 2;
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Convert the contents of columns with dates to MATLAB datetimes using date format string
my_dates = nan(size(dataArray{1},1),1);
for i = 1:size(dataArray{1},1)
    my_dates(i) = datenum(dataArray{1}(i,1), 'yyyy-mm-dd HH:MM:SS.FFF');
end

% Allocate imported array to column variable names
data.Datetime = my_dates;
data.Heightrange = dataArray{2};
data.Hsp = dataArray{3};
data.Vsp = dataArray{4};
data.Wdir = dataArray{5};
data.Turb = dataArray{6};
data.Minintensity = dataArray{7};
data.Meanintensity = dataArray{8};
data.NE = dataArray{9};
data.ES = dataArray{10};
data.SW = dataArray{11};
data.WN = dataArray{12};
data.Vector = dataArray{13};

% Concatenate the matrices

end