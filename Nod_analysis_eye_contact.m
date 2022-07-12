%目線の解析
%ファイル読みこみ
[file,path] = uigetfile('*.csv');
data_file = fullfile(path,file);
in_file=data_file;
[~,name,~] = fileparts(in_file);
t1=readtable(data_file);
h1=height(t1);
fstart = 1;
fend = h1;

M = csvread(data_file, fstart, 0, [fstart, 0, fend, 429]);
%M = csvread(filename, 1, 0);
frame = M(:,1);
time = M(:,3);
confidence = M(:,4);
ganglex = M(:,6);
gangley = M(:,7);
fangley = M(:,11);
fanglex = M(:,12);

renum = fend-fstart;

%%%目線%%%

%アイコンタクト角度(x,y)入力(写真撮ってその角度を記録しておく)
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

%%%点数算出%%%

extime = (renum) / 30;
eyeper = gazesum * 100 / extime;
eyepoint = eyeper * 100 / 65;  %65％を100点とした

    if (eyepoint >= 60.0)   %60点以上で合格=1それ未満で不合格=0
        eyedet = 1;
    else
        eyedet = 0;
    end

%%%グラフ&結果出力%%%

figure('Name','Eye contact')
plot(time(1:renum,1),ganglex(1:renum,1));
hold on
plot(time(1:renum,1),gangley(1:renum,1));
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('yaw','pitch');
hold off;
format4='_Gaze angle.fig';
figure4=strcat(name,format4);
saveas(gcf,figure4);

figure('Name','Eye contact')
plot(time(1:renum,1),ganglex(1:renum,1));
hold on
plot(time(1:renum,1),gangley(1:renum,1));
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
