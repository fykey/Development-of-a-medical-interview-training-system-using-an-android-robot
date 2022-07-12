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

volume  = M(:,430);      %15行目にvolumeの数値があるか確認する

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

figure('Name','Angular velocity') %(up-down)-(left-right)
plot(time(1:renum,1),fvangle)
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
legend('ω');
format0='_Angular velocity(Y-X).fig';
figure0=strcat(name,format0);
saveas(gcf,figure0);

%閾値算出 1.15から2へ        %数人の動画撮って変えてみる
f = mode(fvangle);
a = mean(fvangle);
m = median(fvangle);    %median→中央値
s = std(fvangle);       %std→標準偏差
pthreshold = m + 2 * s;
mthreshold = m - 2 * s;

figure('Name','Angular velocity') %連続除外
plot(time(1:renum,1),fvangle)
hline = refline([0 pthreshold]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
hline2 = refline([0 mthreshold]);     %閾値のライン出力
hline2.Color = 'b';                  %その色青
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
legend('ω','Threshold(+)','Threshold(-)');
format9='_Angular velocity(before).fig';
figure9=strcat(name,format9);
saveas(gcf,figure9);

%連続したうなずきは一回の相槌とする(重複する閾値超えを減らす)
for i = 1:renum-30
    if fvangle(i,1) > pthreshold
        for j = i+1:i+30
            if fvangle(j,1) > pthreshold
                fvangle(j,1) = 0;
            end
        end
    elseif fvangle(i,1) < mthreshold
        for j = i+1:i+30
            if fvangle(j,1) < mthreshold
                fvangle(j,1) = 0;
            end
        end
    end
end

figure('Name','Angular velocity') %連続除外
plot(time(1:renum,1),fvangle)
hline = refline([0 pthreshold]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
hline2 = refline([0 mthreshold]);     %閾値のライン出力
hline2.Color = 'b';                  %その色青
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
legend('ω','Threshold(+)','Threshold(-)');
format5='_Angular velocity(continuous).fig';
figure5=strcat(name,format5);
saveas(gcf,figure5);

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

figure('Name','Voice volume')
plot(time(1:renum+1,1),volume)
hline = refline([0 0.05]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
xlim([0,20]);
xlabel('Time [s]');
ylabel('Amplitude [V]');
legend('Voice amplitude','Threshold');
format6='_Voice volume.fig';
figure6=strcat(name,format6);
saveas(gcf,figure6);

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

%%%点数算出%%%

extime = (renum) / 30;
nodper = nodsum / extime;
nodpoint = nodper * 100 / 0.226;  %0.226回/1秒を100点とした
    
    if (nodpoint >= 60.0)   %60点以上で合格=1それ未満で不合格=0
        noddet = 1;
    else
        noddet = 0;
    end

%%%グラフ&結果出力%%%

figure('Name','Nod counts')         %うなずき
plot(time(1:renum,1),fvangle)
hline = refline([0 pthreshold]);     %閾値のライン出力
hline.Color = 'r';                  %その色赤
hline2 = refline([0 mthreshold]);     %閾値のライン出力
hline2.Color = 'b';                  %その色青
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
legend('ω','Threshold(+)','Threshold(-)');
format1='_Nod counts.fig';
figure1=strcat(name,format1);
saveas(gcf,figure1);

figure('Name','Angular velocity of up-down') %up-down
plot(time(1:renum,1),fvangleY)
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angular velocity(pitch) [deg/s]');
legend('ωp');
format2='_Angular velocity of up-down.fig';
figure2=strcat(name,format2);
saveas(gcf,figure2);

figure('Name','Angular velocity of left-right') %left-right
plot(time(1:renum,1),fvangleX)
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angular velocity(yaw) [deg/s]');
legend('ωy');
format3='_Angular velocity of left-right.fig';
figure3=strcat(name,format3);
saveas(gcf,figure3);

figure('Name','Speaking or not')
plot(time(1:renum,1),speaking)
bar(time(1:renum,1),speaking)
xlim([0,20]);
xlabel('Time [s]');
ylabel('Speaking or not');
format4='_Speaking or not.fig';
figure4=strcat(name,format4);
saveas(gcf,figure4);

figure('Name','Angle(pitch)')
plot(time(1:renum+1,1),fangley)
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angle(pitch) [deg]');
format7='_Angle(pitch).fig';
figure7=strcat(name,format7);
saveas(gcf,figure7);

figure('Name','Angle(pitch)')
plot(time(1:renum+1,1),fanglex)
xlim([0,20]);
xlabel('Time [s]');
ylabel('Angle(yaw) [deg]');
format8='_Angle(yaw).fig';
figure8=strcat(name,format8);
saveas(gcf,figure8);

X1 = sprintf('頷いた回数　　　    %d回\n', nodsum);
X2 = sprintf('最頻値　　%f\n平均値　　%f\n中央値　　%f\n標準偏差　%f\n+閾値 　　%f\n-閾値　　%f\n',f,a,m,s,pthreshold,mthreshold);
X3 = sprintf('実験時間　%0.1f秒\n１秒毎の頷き回数　%0.1f回　点数　%0.1f点　判定　%d\n'...
    ,extime,nodper,nodpoint,noddet);

%f = msgbox(X1);
%f = msgbox(X2);
%f = msgbox(X3);

rformat='_Result.txt';
result_file=strcat(name,rformat);
writematrix(X1,result_file,'WriteMode','append')
writematrix(X2,result_file,'WriteMode','append')
writematrix(X3,result_file,'WriteMode','append')

disp('動作解析終了！');
