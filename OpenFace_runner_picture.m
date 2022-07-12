% A demo script that demonstrates how to process a single video file using
% OpenFace and extract and visualize all of the features

% The location executable will depend on the OS
if(isunix)
    executable = '"C:\Users\Itaru Minami\Downloads\OpenFace(編集したやつ)\FeatureExtraction"';
else
    executable = '"C:\Users\Itaru Minami\Downloads\OpenFace(編集したやつ)\FeatureExtraction.exe"';
end

% Input file
[file,path] = uigetfile('*.jpg');
picture_file = fullfile(path,file);
in_file = picture_file;

% Where to store the output
output_dir = pwd;

% This will take file after -f and output all the features to directory
% after -out_dir
command = sprintf('%s -f "%s" -out_dir "%s" -verbose', executable, in_file, output_dir);
                 
if(isunix)
    unix(command);
else
    dos(command);
end

%% Demonstrating reading the output files

% Most of the features will be in the csv file in the output directory with
% the same name as the input file
[~,name2,~] = fileparts(in_file);
output_csv = sprintf('%s/%s.csv', output_dir, name2);

% First read in the column names, to know which columns to read for
% particular features
tab = readtable(output_csv);
column_names = tab.Properties.VariableNames;

% Read all of the data
all_params  = dlmread(output_csv, ',', 1, 0);

% This indicates which frames were succesfully tracked

% Find which column contains success of tracking data and timestamp data
valid_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'success'));
time_stamp_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'timestamp'));

% Extract tracking success data and only read those frame
valid_frames = logical(all_params(:,valid_ind));

% Get the timestamp data
time_stamps = all_params(valid_frames, time_stamp_ind);

disp('CSV出力完了！');
