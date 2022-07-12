%すでに解析したデータをもちいてつかう
% csv, _Moment.csvを用意
% 
% 

[file,path] = uigetfile('*.csv');   %choose movie
csv_file = fullfile(path,file);
in_file = csv_file;

confidences_threshold = 0;
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





% Output HOG files
%output_hog_file = sprintf('%s/%s.hog', output_dir, name);
%[hog_data, valid_inds] = Read_HOG_file(output_hog_file);

% Output aligned images
%{
output_aligned_dir = sprintf('%s/%s_aligned/', output_dir, name);
img_files = dir([output_aligned_dir, '/*.bmp']);
imgs = cell(numel(img_files, 1));
%for i=1:numel(img_files)
%   imgs{i} = imread([ output_aligned_dir, '/', img_files(i).name]);
%   imshow(imgs{i})
%   drawnow
%end
%}

%% Voice_analysis2021.m
%音声の解析
%wavファイル読み込み
%{
vformat='.mp3';
voice_file=strcat(name,vformat);
[x,fs] = audioread(voice_file);
t = (0:length(x)-1)/fs;

%j=1612;
j=1612;
volumeMatrix = zeros(round(length(x)/j),2);

%音量の平均の算出
volume_sum = 0;
v_count = 0;

for i=1:length(x)
    x(i,1) = abs(x(i,1));
    if rem(i,j)==0
        volumeMatrix(round(i/j),1) = t(i);
        volumeMatrix(round(i/j),2) = x(i,1);
    end
    if(x(i,1) > 0.1)
        volume_sum = volume_sum + x(i,1);
        v_count = v_count + 1;
    end
end

volume_ave = volume_sum / v_count;          %音量の平均値

fig = figure;
plot(t,x);
%hline = refline([0 volume_ave]);
hline = refline([0 0.05]);
hline.Color = 'r';
xlabel('Time[s]');
ylabel('Amplitude[V]');
ylim([0,1.0]);
xlim([0,10]);
vformat='_Volume.fig';
vfigure=strcat(name,vformat);
saveas(gcf,vfigure);

volume_file='volume.csv';
writematrix(volumeMatrix,volume_file)
d1=readtable(volume_file);
d2=readtable('add.csv');
D=[d1;d2];
writetable(D,volume_file)

vX1=sprintf('声の大きさ              %.2fV\n', volume_ave);
%f = msgbox(vX1);

rformat='_Result.txt';
result_file=strcat(name,rformat);
writematrix(vX1,result_file,'WriteMode','append')

%% Csv_copypaste.m
%OpenFaceのcsv429,430行目にvolume.csvの2行をコピー&ペースト、上書き保存
cformat='.csv';
csv_file=strcat(name,cformat);
t1=readtable(csv_file);
t2=readtable('volume.csv');
a=readtable('zero.csv');
h1=height(t1);
h2=height(t2);

if h1>h2
    while h1-h2>0
        t2=[t2;a];
        writetable(t2,volume_file)
        h2=h2+1;
    end
elseif h1<h2
    while h2-h1>0
        t2([h2],:)=[];
        h2=h2-1;
    end
else
end

T=[t1 t2];
writetable(T,csv_file)

%}
%% OpenFace_runner_picture.m
% A demo script that demonstrates how to process a single video file using
% OpenFace and extract and visualize all of the features

% The location executable will depend on the OS
%{
if(isunix)
    executable = "C:\Program Files\OpenFace\FeatureExtraction";
else
    executable = "C:\Program Files\OpenFace\FeatureExtraction.exe";
end
%}
% Input file
jformat = '_Moment.jpg';
jpg_file = strcat(name,jformat);
in_file = jpg_file;

% Where to store the output
output_dir = pwd;

% This will take file after -f and output all the features to directory
% after -out_dir
%command = sprintf('%s -f "%s" -out_dir "%s" -verbose', executable, in_file, output_dir);
                 


% Demonstrating reading the output files

% Most of the features will be in the csv file in the output directory with
% the same name as the input file
[~,name2,~] = fileparts(in_file);
output_csv_Moment = sprintf('%s/%s.csv', output_dir, name2);

% First read in the column names, to know which columns to read for
% particular features
tab_Moment = readtable(output_csv_Moment);
column_names_Moment = tab_Moment.Properties.VariableNames;

% Read all of the data
all_params_Moment  = dlmread(output_csv_Moment, ',', 1, 0);

% This indicates which frames were succesfully tracked

% Find which column contains success of tracking data and timestamp data
valid_ind_Moment = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names_Moment, 'success'));
time_stamp_ind_Moment = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names_Moment, 'timestamp'));

% Extract tracking success data and only read those frame
valid_frames_Moment = logical(all_params_Moment(:,valid_ind_Moment));

% Get the timestamp data
time_stamps_Moment = all_params_Moment(valid_frames_Moment, time_stamp_ind_Moment);

%% Nod_analysis22021.m
%目線の解析
%ファイル読みこみ
dformat='.csv';%%
data_file=strcat(name,dformat);
t1=readtable(data_file);
h1=height(time_stamps);
fstart = 1;
fend = h1;

[file,path] = uigetfile('*.csv');   %choose movie
movie_file = fullfile(path,file);
in_file = movie_file;

% Where to store the output
output_dir = pwd;


% Most of the features will be in the csv file in the output directory with
% the same name as the input file
[~,name,~] = fileparts(in_file);
output_csv = sprintf('%s/%s.csv', output_dir, name);

% First read in the column names, to know which columns to read for
% particular features
tab = readtable(output_csv);
column_names = tab.Properties.VariableNames;

% Read all of the data
all_params  = dlmread(output_csv, ',', 1, 0);


%{

M = csvread(data_file, fstart, 0, [fstart, 0, fend, 429]);
%M = csvread(filename, 1, 0);
frame = M(:,1);
time = M(:,3);
confidence = M(:,4);
ganglex = M(:,6);
gangley = M(:,7);
fangley = M(:,11);
fanglex = M(:,12);
%}

frame = valid_frames;
time = time_stamps;

ganglex = gaze(:, 1);
gangley = gaze(:, 2);

volume_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'Var'));
volume_params = all_params(valid_frames, volume_inds);


%特徴点

%眉
%{
X_21=M(:,171);X_22=M(:,172);
Y_21=M(:,239);Y_22=M(:,240);
Z_21=M(:,307);Z_22=M(:,308);
%}

X_21 = Xs(:, 22);X_22 = Xs(:, 23);
Y_21 = Ys(:, 22);Y_22 = Ys(:, 23);
Z_21 = Zs(:, 22);Z_22 = Zs(:, 23);

%鼻
%{
X_27=M(:,177);X_30=M(:,180);
Y_27=M(:,245);Y_30=M(:,248);
Z_27=M(:,313);Z_30=M(:,316);
%}

X_27 = Xs(:, 28);X_30 = Xs(:, 31);
Y_27 = Ys(:, 28);Y_30 = Ys(:, 31);
Z_27 = Zs(:, 28);Z_30 = Zs(:, 31);

%口
%{
X_61=M(:,211);X_63=M(:,213);
Y_61=M(:,279);Y_63=M(:,281);
Z_61=M(:,347);Z_63=M(:,349);
%}

X_61 = Xs(:, 62);X_63 = Xs(:, 64);
Y_61 = Ys(:, 62);Y_63 = Ys(:, 64);
Z_61 = Zs(:, 62);Z_63 = Zs(:, 64);



volume  = volume_params(:,2);      %15行目にvolumeの数値があるか確認する

renum = fend-fstart;

%%%頷き%%%

fvangleX = zeros(renum,1);
fvangleY = zeros(renum,1);
fvangle = zeros(renum,1);
speaking = zeros(renum,1);

%角速度算出
for i = 1:(renum)
    fvangleX(i,1) = (fanglex(i+1,1) - fanglex(i,1)) / (time(i+1,1) - time(i,1));
    fvangleY(i,1) = (fangley(i+1,1) - fangley(i,1)) / (time(i+1,1) - time(i,1));
end

for i = 1:(renum)
    if fvangleY(i,1) >= 0
        fvangle(i,1) = fvangleY(i,1) - abs(fvangleX(i,1));
        if fvangle(i,1) < 0
            fvangle(i,1) = 0;
        end
    else
        fvangle(i,1) = fvangleY(i,1) + abs(fvangleX(i,1));
        if fvangle(i,1) > 0
            fvangle(i,1) = 0;
        end
    end
end

%大きすぎる角速度を0に
for i = 1:renum
    if fvangle(i,1) > 1000 || fvangle(i,1) < -1000
        fvangle(i,1) = 0;
    end
end

%閾値算出 1.15から2へ        %数人の動画撮って変えてみる
m = median(fvangle);    %median→中央値
s = std(fvangle);       %std→標準偏差
pthreshold = m + 2 * s;
mthreshold = m - 2 * s;

%連続したうなずきは一回の相槌とする(重複する閾値超えを減らす)
for i = 1:renum-30
    if fvangle(i,1) > pthreshold
        for j = i+1:i+30
            if fvangle(j,1) > pthreshold
                fvangle(j,1) = pthreshold;
            end
        end
    elseif fvangle(i,1) < mthreshold
        for j = i+1:i+30
            if fvangle(j,1) < mthreshold
                fvangle(j,1) = mthreshold;
            end
        end
    end
end

%医師が話しているときは角速度を0にする("はい"や"なるほど"等の短い返事はそのまま)
for i = 1:renum-30
    if volume(i,1) > 0.05
        for j = i+1:i+30
            if volume(j,1) > 0.05
                for k = i:j
                    fvangle(k,1) = 0;
                    speaking(k,1) = 100;
                end
            end
        end
    end
end

%回数算出
nodsum = 0;
nodtimes = zeros(1, 1);

i = 1;
k = 1;
while i < (renum-1)
    if ((fvangle(i,1) <= pthreshold) && (fvangle(i+1,1) > pthreshold))
        for j = 1:(renum-1)
            if ((i < j) && (j < i + 30))
                if ((fvangle(j,1) >= mthreshold) && (fvangle(j+1,1) < mthreshold))
                    nodtimes(1, k) = time_stamps(j, 1);
                    k = k + 1;


                    nodsum = nodsum + 1;
                    i = j;
                    break
                end
            end
        end
    end
    i = i + 1;
end

%%%目線%%%

%アイコンタクト角度(x,y)入力(写真撮っておく)
pformat='_Moment.csv';
picture_file = strcat(name,pformat);
N = csvread(picture_file, 1, 0, [1, 0, 1, 427]);

left_right = N(:,6);
up_down = N(:,7);

x = left_right(1,1);
y = up_down(1,1);

%アイコンタクト閾値
con_threshold = 8;
max_con_x = x + con_threshold;
min_con_x = x - con_threshold;
max_con_y = y + con_threshold;
min_con_y = y - con_threshold;

%秒数算出(アイコンタクト)
gazesum = 0;
for i = 1:(renum)
    if min_con_x < ganglex(i,1) && ganglex(i,1) < max_con_x
        if min_con_y < gangley(i,1) && gangley(i,1) < max_con_y
            gazesum = (time(i+1,1) - time(i,1)) + gazesum;
        end
    end
end

%%%表情%%%

%特徴点距離(眉は鼻頭(30)から、口は鼻尻(27)からの距離)
d_21 = zeros(renum,1);
d_22 = zeros(renum,1);
d_27 = zeros(renum,1);
d_30 = zeros(renum,1);
d_61 = zeros(renum,1);
d_63 = zeros(renum,1);

for i = 1:(renum)
    d_21(i,1) = sqrt((X_21(i,1) - X_30(i,1)).^2 + (Y_21(i,1) - Y_30(i,1)).^2 + (Z_21(i,1) - Z_30(i,1)).^2);
    d_22(i,1) = sqrt((X_22(i,1) - X_30(i,1)).^2 + (Y_22(i,1) - Y_30(i,1)).^2 + (Z_22(i,1) - Z_30(i,1)).^2);
    d_27(i,1) = sqrt((X_27(i,1) - X_30(i,1)).^2 + (Y_27(i,1) - Y_30(i,1)).^2 + (Z_27(i,1) - Z_30(i,1)).^2);
    d_30(i,1) = sqrt((X_30(i,1) - X_27(i,1)).^2 + (Y_30(i,1) - Y_27(i,1)).^2 + (Z_30(i,1) - Z_27(i,1)).^2);
    d_61(i,1) = sqrt((X_61(i,1) - X_27(i,1)).^2 + (Y_61(i,1) - Y_27(i,1)).^2 + (Z_61(i,1) - Z_27(i,1)).^2);
    d_63(i,1) = sqrt((X_63(i,1) - X_27(i,1)).^2 + (Y_63(i,1) - Y_27(i,1)).^2 + (Z_63(i,1) - Z_27(i,1)).^2);
end

%特徴点移動距離(中央値からの)
m_21 = zeros(renum,1);
m_22 = zeros(renum,1);
m_27 = zeros(renum,1);
m_30 = zeros(renum,1);
m_61 = zeros(renum,1);
m_63 = zeros(renum,1);

for i = 1:(renum)
    m_21(i,1) = d_21(i,1) - median(d_21);
    m_22(i,1) = d_22(i,1) - median(d_22);
    m_27(i,1) = d_27(i,1) - median(d_27);
    m_30(i,1) = d_30(i,1) - median(d_30);
    m_61(i,1) = d_61(i,1) - median(d_61);
    m_63(i,1) = d_63(i,1) - median(d_63);
end

%特徴点表情(特徴点距離比)
smile = zeros(renum+1,1);
surprise = zeros(renum+1,1);

for i = 1:(renum)
    smile(i,1) = (m_61(i,1) + m_63(i,1)) / abs(2 * m_30(i,1));
    surprise(i,1) = (m_21(i,1) + m_22(i,1)) / abs(2 * m_27(i,1));
end

%大きすぎる特徴点距離比を0に
for i = 1:renum
    if smile(i,1) > 50 || smile(i,1) < -50
        smile(i,1) = 0;
    end
    if surprise(i,1) > 50 || surprise(i,1) < -50
        surprise(i,1) = 0;
    end
end

%秒数算出
smilesum = 0;
surprisesum = 0;

for i = 1:(renum)
    if smile(i,1) < -2
        smilesum = (time(i+1,1) - time(i,1)) + smilesum;
    end
    
    if surprise(i,1) > 2
        surprisesum = (time(i+1,1) - time(i,1)) + surprisesum;
    end
end

%%%点数算出%%%

extime = (renum) / 30;
nodper = nodsum / extime;
nodpoint = nodper * 100 / 0.226;  %0.226回/1秒を100点とした
eyeper = gazesum * 100 / extime;
eyepoint = eyeper * 100 / 65;  %65％を100点とした
faceper = (smilesum + surprisesum) * 100 / extime;
facepoint = faceper * 100 / 50;  %50％を100点とした

    if (nodpoint >= 60.0)   %60点以上で合格=1それ未満で不合格=0   
        noddet = 1;
    else
        noddet = 0;
    end
    
    if (eyepoint >= 60.0)
        eyedet = 1;
    else
        eyedet = 0;
    end
    
    if (facepoint >= 60.0)   
        facedet = 1;
    else
        facedet = 0;
    end

%%%グラフ&結果出力%%%

figure('Name','Nod counts')         %うなずき
plot(time(1:renum,1),fvangle)
hline = refline([0 pthreshold]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
hline2 = refline([0 mthreshold]);     %閾値のライン出力
hline2.Color = 'b';                  %その色青
xlim([0,30]);
ylim([-50,50]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
format1='_Nod counts.fig';
figure1=strcat(name,format1);
saveas(gcf,figure1);

figure('Name','Angular velocity of pitch') %up-down
plot(time(1:renum,1),fvangleY)
xlim([0,30]);
ylim([-50,50]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
format2='_Angular velocity of pitch.fig';
figure2=strcat(name,format2);
saveas(gcf,figure2);

figure('Name','Angular velocity of yaw') %left-right
plot(time(1:renum,1),fvangleX)
xlim([0,30]);
ylim([-50,50]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
format3='_Angular velocity of yaw.fig';
figure3=strcat(name,format3);
saveas(gcf,figure3);

figure('Name','Speaking or not')
plot(time(1:renum,1),speaking)
bar(time(1:renum,1),speaking)
xlim([0,30]);
ylim([0,100]);
xlabel('Time [s]');
ylabel('Speaking or not');
format4='_Speaking or not.fig';
figure4=strcat(name,format4);
saveas(gcf,figure4);

figure('Name','Eye contact')
plot(time(1:renum,1),ganglex(1:renum,1));
hold on
plot(time(1:renum,1),gangley(1:renum,1));
hline = refline([0 max_con_x]);
hline.Color = 'b';
hline2 = refline([0 max_con_y]);
hline2.Color = 'r';
hline3 = refline([0 min_con_x]);
hline3.Color = 'b';
hline4 = refline([0 min_con_y]);
hline4.Color = 'r';
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('Left - right','Up - down');
hold off;
format5='_Eye contact.fig';
figure5=strcat(name,format5);
saveas(gcf,figure5);

figure('Name','Smile or not')       %微笑み
plot(time(1:renum,1),smile(1:renum,1))
hline = refline([0 -2]);             %閾値のライン出力
hline.Color = 'b';                  %その色赤
xlim([0,30]);
ylim([-1,1]);
xlabel('Time');
ylabel('Facial Expression');
format6='_Smile or not.fig';
figure6=strcat(name,format6);
saveas(gcf,figure6);

figure('Name','Surprise or not')    %驚き
plot(time(1:renum,1),surprise(1:renum,1))
hline = refline([0 2]);             %閾値のライン出力
hline.Color = 'r';                  %その色赤
xlim([0,30]);
ylim([-1,1]);
xlabel('Time');
ylabel('Facial Expression');
format7='_Surprise or not.fig';
figure7=strcat(name,format7);
saveas(gcf,figure7);

X1 = sprintf('頷いた回数　　　%d回\n患者を診た時間　%.1f秒\n表情あり　　　　%.1f秒\n', nodsum,gazesum,smilesum + surprisesum);
X2 = sprintf('中央値　　%f\n標準偏差　%f\n+閾値 　　%f\n-閾値　　%f\n',m,s,pthreshold,mthreshold);
X3 = sprintf('Left-lightアイコンタクト　%.3f°\nUp-downアイコンタクト　　 %.3f°\n',x,y);
X4 = sprintf('微笑み　%.1f秒\n驚き　　%.1f秒\n',smilesum,surprisesum);
X5 = sprintf('実験時間　　　　　%0.1f秒\n１秒毎の頷き回数　%0.3f回　点数　%0.1f点　判定　%d\nアイコンタクト　　%0.1f％　 点数　%0.1f点　判定　%d\n表情　　　　　　　%0.1f％　 点数　%0.1f点　判定　%d\n'...
    ,extime,nodper,nodpoint,noddet,eyeper,eyepoint,eyedet,faceper,facepoint,facedet);

%f = msgbox(X1);
%f = msgbox(X2);
%f = msgbox(X3);
%f = msgbox(X4);
%f = msgbox(X5);

rformat='_Result.txt';
result_file=strcat(name,rformat);
writematrix(X1,result_file,'WriteMode','append')
writematrix(X2,result_file,'WriteMode','append')
writematrix(X3,result_file,'WriteMode','append')
writematrix(X4,result_file,'WriteMode','append')
writematrix(X5,result_file,'WriteMode','append')



%% Move_file.m
%新しくフォルダを作成し，撮影した動画，解析結果，グラフを分かりやすいようにまとめる
folder_name = strcat('confidences_threshold_',string(confidences_threshold), '_', name);
mkdir(folder_name);
movefile(strcat(name,'*'),folder_name);

disp('動作解析終了！');
