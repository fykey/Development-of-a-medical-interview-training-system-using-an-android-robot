[file,path] = uigetfile('*.csv');   %choose movie
csv_file = fullfile(path,file);
in_file = csv_file;
confidences_threshold = 0;
noise_param = 5;
i = 1;
j = 1;
% Where to store the output
output_dir = pwd;
% Demonstrating reading the output files

% Most of the features will be in the csv file in the output directory with
% the same name as the input file
[~,name,~] = fileparts(in_file);
output_csv = sprintf('%s/%s.csv', output_dir, name);

% First read in teatures
tab = readtable(output_csv);
% particular f(output_csv);
column_names = tab.Properties.VariableNames;

% Read all of the data
all_params  = dlmread(output_csv, ',', 1, 0);

% This indicates which frames were succesfully tracked

% Find which column contains success of tracking data and timestamp data
valid_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'success'));
time_stamp_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'timestamp'));
%confidences
confidences_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'confidence'));
for j = 1:height(all_params)
    if(all_params(j, confidences_inds) >= confidences_threshold)
        all_params(j, confidences_inds) = true;
    else
        all_params(j, confidences_inds) = false;
    end
end


% Extract tracking success data and only read those frame
valid_cells = logical(all_params(:,valid_ind));
confidences_cells = logical(all_params(:, confidences_inds));

for i = 1:height(all_params)
    if (valid_cells(i, :) == true) && (confidences_cells(i, :) == true)
       all_params(i, valid_ind) = true;
    else
        all_params(i, valid_ind) = false;
    end
end

%confidences = all_params(valid_frames, confidences_inds);
valid_frames = logical(all_params(:,valid_ind));


% Get the timestamp data
time_stamps = all_params(valid_frames, time_stamp_ind);
data_len = length(time_stamps);
% Finding which header line starts with p_ (basically model params)
shape_inds = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'p_'));

% Output rigid (first 6) and non-rigid shape parameters
shape_params  = all_params(valid_frames, shape_inds);

%{
figure
plot(time_stamps, shape_params);
title('Shape parameters');
xlabel('Time (s)');
%}

% Demonstrate 2D landmarks
landmark_inds_x = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'x_'));
landmark_inds_y = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'y_'));

xs = all_params(valid_frames, landmark_inds_x);
ys = all_params(valid_frames, landmark_inds_y);

%{
figure
for j = 1:size(xs,1)
    plot(xs(j,:), -ys(j,:), '.');
    hold on;
    plot(eye_xs(j,:), -eye_ys(j,:), '.r');
    hold off;
    
    xlim([min(xs(1,:)) * 0.5, max(xs(2,:))*1.4]);
    ylim([min(-ys(1,:)) * 1.4, max(-ys(2,:))*0.5]);
    xlabel('x (px)');
    ylabel('y (px)');
    drawnow
end
%}

% Demonstrate 3D landmarks
landmark_inds_x = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'X_'));
landmark_inds_y = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'Y_'));
landmark_inds_z = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'Z_'));

Xs = all_params(valid_frames, landmark_inds_x);
Ys = all_params(valid_frames, landmark_inds_y);
Zs = all_params(valid_frames, landmark_inds_z);



% Demonstrate AUs
au_reg_inds = cellfun(@(x) ~isempty(x) && x==5, strfind(column_names, '_r'));

au_class_inds = cellfun(@(x) ~isempty(x) && x==5, strfind(column_names, '_c'));

aus_c = all_params(valid_frames, au_class_inds);
figure
plot(time_stamps, aus_c);
title('Facial Action Units (presense)');
xlabel('Time (s)');
ylim([0,2]);

% Demo pose
%pose_inds = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'pose_'));

%pose = all_params(valid_frames, pose_inds);
%figure
%plot(time_stamps, pose);
%title('Pose (rotation and translation)');
%xlabel('Time (s)');

% Demo gaze
gaze_inds = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'gaze_angle'));

% Read gaze (x,y,z) for one eye and (x,y,z) for another
gaze = all_params(valid_frames, gaze_inds);

plot(time_stamps, gaze(:,1), 'DisplayName', 'Left - right');
hold on;
plot(time_stamps, gaze(:,2), 'DisplayName', 'Up - down');
xlabel('Time(s)') % x-axis label
ylabel('Angle radians') % y-axis label
legend('show');
hold off;
oformat2='_Left right_up down.fig';
ofigure2=strcat(name,oformat2);
saveas(gcf,ofigure2);

%face_angle(turn, up_down)

fangleys_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'up_down'));
fangley = all_params(valid_frames, fangleys_inds);

fabglexs_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'turn'));
fanglex = all_params(valid_frames, fabglexs_inds);

a = 1;
b = 1;

for n = 1:data_len
    if fanglex(n, 1) > 0
        pfanglex(a, 1) = fanglex(n, 1);
        a = a + 1;
    else
        mfanglex(b, 1) = fanglex(n, 1);
        b = b + 1;
    end
end

case1 = zeros(data_len + 1, 1);
case2 = zeros(data_len + 1, 1);
m = median(fanglex);    %median→中央値
s = std(fanglex);       %std→標準偏差

pave = mean(pfanglex);
mave = mean(mfanglex);
mean1 = mean(fanglex);

ave = (pave + mave)/2;

pthreshold = ave + s / 2;
mthreshold = ave - s / 2;


%turn_threshold = 15;
a = 1;
b = 1;
for k = 1 : data_len

    if(fanglex(k, 1) > pthreshold)
        case1(k, 1) = 50;

    else if (fanglex(k, 1) < mthreshold)
            case2(k, 1) = 50;

        end
    end

end






