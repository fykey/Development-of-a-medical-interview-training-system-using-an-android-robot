%音声の解析
%wavファイル読み込み
[file,path] = uigetfile('*.mp3');
voice_file = fullfile(path,file);
in_file=voice_file;
[~,name,~] = fileparts(in_file);
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
