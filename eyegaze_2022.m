clear

% The location executable will depend on the OS
if(isunix)
    executable = '"C:\Users\fuyuk\Documents\研究_大野\OpenFace(編集したやつ)\FeatureExtraction.exe"';
else
    executable = '"C:\Users\fuyuk\Documents\研究_大野\OpenFace(編集したやつ)\FeatureExtraction"';
end

% Input file
[file,path] = uigetfile('*.mp4');   %choose movie
movie_file = fullfile(path,file);
in_file = movie_file;

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

confidences_threshold = 0.85;
% Where to store the output
output_dir = pwd;
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
windows環境では出力されない？(R4.4.18 大野)
eye_landmark_inds_x = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'eye_lmk_x_'));
eye_landmark_inds_y = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'eye_lmk_y_'));

eye_xs = all_params(valid_frames, eye_landmark_inds_x);
eye_ys = all_params(valid_frames, eye_landmark_inds_y);
%}

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
%{
windows環境では出力されない？(R4.4.18 大野)
eye_landmark_inds_x = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'eye_lmk_X_'));
eye_landmark_inds_y = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'eye_lmk_Y_'));
eye_landmark_inds_z = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'eye_lmk_Z_'));

eye_xs = all_params(valid_frames, eye_landmark_inds_x);
eye_ys = all_params(valid_frames, eye_landmark_inds_y);
eye_zs = all_params(valid_frames, eye_landmark_inds_z);
%}

%{
figure
for j = 1:size(xs,1)
    plot3(xs(j,:), ys(j,:), zs(j,:), '.');axis equal;
    hold on;
    plot3(eye_xs(j,:), eye_ys(j,:), eye_zs(j,:), '.r');
    hold off;
    xlabel('X (mm)');
    ylabel('Y (mm)');    
    zlabel('Z (mm)');    
    drawnow
end
%}

% Demonstrate AUs
au_reg_inds = cellfun(@(x) ~isempty(x) && x==5, strfind(column_names, '_r'));

aus_r = all_params(valid_frames, au_reg_inds);
figure
plot(time_stamps, aus_r);
title('Facial Action Units (intensity)');
xlabel('Time (s)');
ylabel('Intensity');
ylim([0,6]);
oformat1='_Action units.fig';
ofigure1=strcat(name,oformat1);
saveas(gcf,ofigure1);

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
ganglex = gaze(:, 1);
gangley = gaze(:, 2);
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
jformat = '_Moment.jpg';
jpg_file = strcat(name,jformat);
in_file = jpg_file;

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




%face_angle(turn, up_down)

fangleys_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'up_down'));
fangley = all_params(valid_frames, fangleys_inds);

fabglexs_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'turn'));
fanglex = all_params(valid_frames, fabglexs_inds);


%移動平均フィルタ

%FPS数/2
fpsper2 = 15;

coeff15flame = ones(1, fpsper2)/fpsper2;

MAFfil_ganglex = filter(coeff15flame, 1, ganglex);
MAFfil_gangley = filter(coeff15flame, 1, gangley);


MAFfil_ganglex = filter(coeff15flame, 1, MAFfil_ganglex);
MAFfil_gangley = filter(coeff15flame, 1, MAFfil_gangley);


fil_time = time_stamps - 2*((fpsper2 - 1) / 2 / (fpsper2 * 2));
time = fil_time;

renum = length(MAFfil_gangley);
pformat='_Moment.csv';
picture_file = strcat(name,pformat);
N = csvread(picture_file, 1, 0, [1, 0, 1, 427]);

left_right = N(:,6);
up_down = N(:,7);

x = left_right(1,1);
y = up_down(1,1);
theta = atan2(y, x);

%アイコンタクト閾値
con_threshold = 8;
theta_threshold = pi()/12;

max_con_x = x + con_threshold;
min_con_x = x - con_threshold;
max_con_y = y + con_threshold;
min_con_y = y - con_threshold;

min_theta = theta - theta_threshold;
max_theta = theta + theta_threshold;


%秒数算出(アイコンタクト)
gazesum = 0;
gazesum2 = 0;
for i = 1:(renum) - 1
    MAFfil_theta(i,1) = atan2(MAFfil_gangley(i,1), MAFfil_ganglex(i,1));
    if min_con_x < MAFfil_ganglex(i,1) && MAFfil_ganglex(i,1) < max_con_x
        if min_con_y < MAFfil_gangley(i,1) && MAFfil_gangley(i,1) < max_con_y
            gazesum = (time(i+1,1) - time(i,1)) + gazesum;
        end
    end

    if min_theta < MAFfil_theta(i, 1) && MAFfil_theta(i, 1) < max_theta;
    gazesum2 = gazesum + (time(i+1,1) - time(i,1));
    
    end
end

%%%点数算出%%%

extime = time(renum) - time(1);
eyeper = gazesum * 100 / extime;
eyepoint = eyeper * 100 / 65;  %65％を100点とした

    if (eyepoint >= 60.0)   %60点以上で合格=1それ未満で不合格=0
        eyedet = 1;
    else
        eyedet = 0;
    end

%%%グラフ&結果出力%%%

figure('Name','Eye contact')
plot(time(1:renum,1),MAFfil_ganglex(1:renum,1));
hold on
plot(time(1:renum,1),MAFfil_gangley(1:renum,1));
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('yaw','pitch');
hold off;
format4='_Gaze angle.fig';
figure4=strcat(name,format4);
saveas(gcf,figure4);


figure("Name", "Eye angle")
plot(time(1:renum - 1, 1), MAFfil_theta(1:renum - 1, 1));
hold on
hline7 = refline([0, theta]);
hline7.Color = "b";
hline7.LineStyle = '--';
hline8 = refline([0, max_theta]);
hline8.Color = "b";
hline9 = refline([0, min_theta]);
hline9.Color = "b";
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angle [deg]');

hold off;

figure('Name','Eye contact')
plot(time(1:renum,1),MAFfil_ganglex(1:renum,1));
hold on
plot(time(1:renum,1),MAFfil_gangley(1:renum,1));
hline5 = refline([0 x]);
hline5.Color = 'b';
hline5.LineStyle = '--';
hline = refline([0 max_con_x]);
hline.Color = 'b';
hline6 = refline([0 y]);
hline6.Color = 'r';
hline6.LineStyle = '--';
hline2 = refline([0 max_con_y]);
hline2.Color = 'r';
hline3 = refline([0 min_con_x]);
hline3.Color = 'b';
hline4 = refline([0 min_con_y]);
hline4.Color = 'r';
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('yaw','pitch','Eye contact(yaw)','Threshold(yaw)','Eye contact(pitch)','Threshold(pitch)');
hold off;
format5='_Eye contact.fig';
figure5=strcat(name,format5);
saveas(gcf,figure5);

X1 = sprintf('患者を診た時間　　 %.1f秒\n', gazesum);
X2 = sprintf('Left-lightアイコンタクト　%.3f°\nUp-downアイコンタクト　　 %.3f°\n',x,y);
X3 = sprintf('実験時間　%0.1f秒\nアイコンタクト　%0.1f％ 点数　%0.1f点　判定　%d\n'...
    ,extime,eyeper,eyepoint,eyedet);

%f = msgbox(X1);
%f = msgbox(X2);
%f = msgbox(X3);

rformat='_Result.txt';
result_file=strcat(name,rformat);
writematrix(X1,result_file,'WriteMode','append')
writematrix(X2,result_file,'WriteMode','append')
writematrix(X3,result_file,'WriteMode','append')

disp('動作解析終了！');
