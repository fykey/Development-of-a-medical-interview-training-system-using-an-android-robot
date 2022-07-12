%目線の解析
%ファイル読みこみ
[file,path] = uigetfile('*.csv');
data_file = fullfile(path,file);
in_file=data_file;
[~,name,~] = fileparts(in_file);
% Where to store the output

t1=readtable(data_file);
h1=height(t1);
fstart = 1;
fend = h1;

% First read in the column names, to know which columns to read for
% particular features
tab = readtable(in_file);
column_names = tab.Properties.VariableNames;

% Read all of the data
all_params  = dlmread(in_file, ',', 1, 0);

% This indicates which frames were succesfully tracked

% Find which column contains success of tracking data and timestamp data
valid_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'success'));
time_stamp_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'timestamp'));

% Extract tracking success data and only read those frame
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

%M = csvread(filename, 1, 0);
frame = valid_frames;
time = time_stamps;

% Demo gaze
gaze_inds = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'gaze_angle'));

% Read gaze (x,y,z) for one eye and (x,y,z) for another
gaze = all_params(valid_frames, gaze_inds);

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

% Demonstrate AUs
au_reg_inds = cellfun(@(x) ~isempty(x) && x==5, strfind(column_names, '_r'));

aus_r = all_params(valid_frames, au_reg_inds);

au_class_inds = cellfun(@(x) ~isempty(x) && x==5, strfind(column_names, '_c'));

aus_c = all_params(valid_frames, au_class_inds);

% Demo gaze
gaze_inds = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'gaze_angle'));

% Read gaze (x,y,z) for one eye and (x,y,z) for another
gaze = all_params(valid_frames, gaze_inds);

%confidences
confidences_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'confidence'));
confidences = all_params(valid_frames, confidences_inds);

%face_angle(turn, up_down)

fangleys_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'up_down'));
fangley = all_params(valid_frames, fangleys_inds);

fabglexs_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'turn'));
fanglex = all_params(valid_frames, fabglexs_inds);


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
    if volume(i,1) > 0.1
        for j = i+1:i+30
            if volume(j,1) > 0.1
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
i = 1;

while i < (renum-1)
    if ((fvangle(i,1) <= pthreshold) && (fvangle(i+1,1) > pthreshold))
        for j = 1:(renum-1)
            if ((i < j) && (j < i + 30))
                if ((fvangle(j,1) >= mthreshold) && (fvangle(j+1,1) < mthreshold))
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

%{
figure('Name','Angular velocity of up-down') %up-down
plot(time(1:renum,1),fvangleY)
hline = refline([0 pthreshold]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
hline2 = refline([0 mthreshold]);     %閾値のライン出力
hline2.Color = 'b';                  %その色青
xlim([0,30]);
ylim([-50,50]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
format2='_Angular velocity of up-down.fig';
figure2=strcat(name,format2);
saveas(gcf,figure2);

figure('Name','Angular velocity of left-right') %left-right
plot(time(1:renum,1),fvangleX)
hline = refline([0 pthreshold]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
hline2 = refline([0 mthreshold]);     %閾値のライン出力
hline2.Color = 'b';                  %その色青
xlim([0,30]);
ylim([-50,50]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
format3='_Angular velocity of left-right.fig';
figure3=strcat(name,format3);
saveas(gcf,figure3);
%}

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

%%ディレクトリ移動
move_dir = strcat('nod_', name);
mkdir(move_dir);
movefile(strcat(name,'*'),move_dir);



disp('動作解析終了！');
