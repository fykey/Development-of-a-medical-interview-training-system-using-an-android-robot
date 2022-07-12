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

%特徴点

%輪郭
X_0=M(:,150);X_1=M(:,151);X_2=M(:,152);X_3=M(:,153);X_4=M(:,154);X_5=M(:,155);X_6=M(:,156);X_7=M(:,157);X_8=M(:,158);
X_9=M(:,159);X_10=M(:,160);X_11=M(:,161);X_12=M(:,162);X_13=M(:,163);X_14=M(:,164);X_15=M(:,165);X_16=M(:,166);
Y_0=M(:,218);Y_1=M(:,219);Y_2=M(:,220);Y_3=M(:,221);Y_4=M(:,222);Y_5=M(:,223);Y_6=M(:,224);Y_7=M(:,225);Y_8=M(:,226);
Y_9=M(:,227);Y_10=M(:,228);Y_11=M(:,229);Y_12=M(:,230);Y_13=M(:,231);Y_14=M(:,232);Y_15=M(:,233);Y_16=M(:,234);
Z_0=M(:,286);Z_1=M(:,287);Z_2=M(:,288);Z_3=M(:,289);Z_4=M(:,290);Z_5=M(:,291);Z_6=M(:,292);Z_7=M(:,293);Z_8=M(:,294);
Z_9=M(:,295);Z_10=M(:,296);Z_11=M(:,297);Z_12=M(:,298);Z_13=M(:,299);Z_14=M(:,300);Z_15=M(:,301);Z_16=M(:,302);

%眉
X_17=M(:,167);X_18=M(:,168);X_19=M(:,169);X_20=M(:,170);X_21=M(:,171);X_22=M(:,172);X_23=M(:,173);X_24=M(:,174);X_25=M(:,175);X_26=M(:,176);
Y_17=M(:,235);Y_18=M(:,236);Y_19=M(:,237);Y_20=M(:,238);Y_21=M(:,239);Y_22=M(:,240);Y_23=M(:,241);Y_24=M(:,242);Y_25=M(:,243);Y_26=M(:,244);
Z_17=M(:,303);Z_18=M(:,304);Z_19=M(:,305);Z_20=M(:,306);Z_21=M(:,307);Z_22=M(:,308);Z_23=M(:,309);Z_24=M(:,310);Z_25=M(:,311);Z_26=M(:,312);

%鼻
X_27=M(:,177);X_28=M(:,178);X_29=M(:,179);X_30=M(:,180);X_31=M(:,181);X_32=M(:,182);X_33=M(:,183);X_34=M(:,184);X_35=M(:,185);
Y_27=M(:,245);Y_28=M(:,246);Y_29=M(:,247);Y_30=M(:,248);Y_31=M(:,249);Y_32=M(:,250);Y_33=M(:,251);Y_34=M(:,252);Y_35=M(:,253);
Z_27=M(:,313);Z_28=M(:,314);Z_29=M(:,315);Z_30=M(:,316);Z_31=M(:,317);Z_32=M(:,318);Z_33=M(:,319);Z_34=M(:,320);Z_35=M(:,321);

%目
X_36=M(:,186);X_37=M(:,187);X_38=M(:,188);X_39=M(:,189);X_40=M(:,190);X_41=M(:,191);
X_42=M(:,192);X_43=M(:,193);X_44=M(:,194);X_45=M(:,195);X_46=M(:,196);X_47=M(:,197);
Y_36=M(:,254);Y_37=M(:,255);Y_38=M(:,256);Y_39=M(:,257);Y_40=M(:,258);Y_41=M(:,259);
Y_42=M(:,260);Y_43=M(:,261);Y_44=M(:,262);Y_45=M(:,263);Y_46=M(:,264);Y_47=M(:,265);
Z_36=M(:,322);Z_37=M(:,323);Z_38=M(:,324);Z_39=M(:,325);Z_40=M(:,326);Z_41=M(:,327);
Z_42=M(:,328);Z_43=M(:,329);Z_44=M(:,330);Z_45=M(:,331);Z_46=M(:,332);Z_47=M(:,333);

%口
X_48=M(:,198);X_49=M(:,199);X_50=M(:,200);X_51=M(:,201);X_52=M(:,202);X_53=M(:,203);X_54=M(:,204);X_55=M(:,205);X_56=M(:,206);X_57=M(:,207);
X_58=M(:,208);X_59=M(:,209);X_60=M(:,210);X_61=M(:,211);X_62=M(:,212);X_63=M(:,213);X_64=M(:,214);X_65=M(:,215);X_66=M(:,216);X_67=M(:,217);
Y_48=M(:,266);Y_49=M(:,267);Y_50=M(:,268);Y_51=M(:,269);Y_52=M(:,270);Y_53=M(:,271);Y_54=M(:,272);Y_55=M(:,273);Y_56=M(:,274);Y_57=M(:,275);
Y_58=M(:,276);Y_59=M(:,277);Y_60=M(:,278);Y_61=M(:,279);Y_62=M(:,280);Y_63=M(:,281);Y_64=M(:,282);Y_65=M(:,283);Y_66=M(:,284);Y_67=M(:,285);
Z_48=M(:,334);Z_49=M(:,335);Z_50=M(:,336);Z_51=M(:,337);Z_52=M(:,338);Z_53=M(:,339);Z_54=M(:,340);Z_55=M(:,341);Z_56=M(:,342);Z_57=M(:,343);
Z_58=M(:,344);Z_59=M(:,345);Z_60=M(:,346);Z_61=M(:,347);Z_62=M(:,348);Z_63=M(:,349);Z_64=M(:,350);Z_65=M(:,351);Z_66=M(:,352);Z_67=M(:,353);

volume  = M(:,430);      %430行目にvolumeの数値があるか確認する

renum = fend-fstart;

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
faceper = (smilesum + surprisesum) * 100 / extime;
facepoint = faceper * 100 / 50;  %50％を100点とした
    
    if (facepoint >= 60.0)   %60点以上で合格=1それ未満で不合格=0
        facedet = 1;
    else
        facedet = 0;
    end

%%%グラフ&結果出力%%%

figure('Name','Feature point distance')
plot(time(1:renum,1),d_21(1:renum,1),time(1:renum,1),d_22(1:renum,1),time(1:renum,1),d_27(1:renum,1),time(1:renum,1),d_61(1:renum,1),time(1:renum,1),d_63(1:renum,1))
xlim([0,20]);
xlabel('Time [s]');
ylabel('3D distance [mm]');
legend('21-30','22-30','27-30','27-61','27-63');
format1='_Feature point distance.fig';
figure1=strcat(name,format1);
saveas(gcf,figure1);

figure('Name','Amount of movement')
plot(time(1:renum,1),(m_21(1:renum,1)+m_22(1:renum,1))/2,time(1:renum,1),(m_61(1:renum,1)+m_63(1:renum,1))/2,time(1:renum,1),m_27(1:renum,1))
xlim([0,20]);
xlabel('Time [s]');
ylabel('Movement distance [mm]');
legend('surprise','smile','nose');
format2='_Amount of movement.fig';
figure2=strcat(name,format2);
saveas(gcf,figure2);

figure('Name','Distance ratio')
plot(time(1:renum,1),surprise(1:renum,1),time(1:renum,1),smile(1:renum,1))
xlim([0,20]);
xlabel('Time [s]');
ylabel('Distance ratio');
legend('surprise','smile');
format3='_Distance ratio.fig';
figure3=strcat(name,format3);
saveas(gcf,figure3);

figure('Name','Smile or not')       %微笑み
colororder([0.8500 0.3250 0.0980])
plot(time(1:renum,1),smile(1:renum,1))
hline = refline([0 -2]);             %閾値のライン出力
hline.Color = 'b';                  %その色赤
xlim([0,20]);
xlabel('Time [s]');
ylabel('Distance ratio');
legend('smile','Threshold');
format6='_Smile or not.fig';
figure6=strcat(name,format6);
saveas(gcf,figure6);

figure('Name','Surprise or not')    %驚き
plot(time(1:renum,1),surprise(1:renum,1))
hline = refline([0 2]);             %閾値のライン出力
hline.Color = 'r';                  %その色赤
xlim([0,20]);
xlabel('Time [s]');
ylabel('Distance ratio');
legend('surprise','Threshold');
format7='_Surprise or not.fig';
figure7=strcat(name,format7);
saveas(gcf,figure7);

X1 = sprintf('表情あり　%.1f秒\n', smilesum + surprisesum);
X2 = sprintf('微笑み　%.1f秒\n驚き　　%.1f秒\n',smilesum,surprisesum);
X3 = sprintf('実験時間　%0.1f秒\n表情　%0.1f％ 点数　%0.1f点　判定　%d\n'...
    ,extime,faceper,facepoint,facedet);

%f = msgbox(X1);
%f = msgbox(X2);
%f = msgbox(X3);

rformat='_Result.txt';
result_file=strcat(name,rformat);
writematrix(X1,result_file,'WriteMode','append')
writematrix(X2,result_file,'WriteMode','append')
writematrix(X3,result_file,'WriteMode','append')

disp('動作解析終了！');
