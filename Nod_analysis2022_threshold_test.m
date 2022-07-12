%{
R4.4.19
閾値をマイナス方向に平行移動した場合，個人でどのくらいの差がでるのかを試験するプログラム．
qの値を0~10に0.2ずつ増加させたとき，うなずき量の推移を試験する．

tracking成功フレームのみ抽出するのは，confidencesなどのパラメータを用いたうえで効果を評価できて
からにする．
%}

decide_value = 1;
loop_num = 0;
nodcount = cell(2, 20);
filename = cell(1, 1);
result_file = 'nod_alldata.csv';



while decide_value == 1
    loop_num = loop_num + 1;
    
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
    filename{loop_num + 1, 1} = data_file;
    
    
    
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
    
    %眉
    X_21=M(:,171);X_22=M(:,172);
    Y_21=M(:,239);Y_22=M(:,240);
    Z_21=M(:,307);Z_22=M(:,308);
    
    %鼻
    X_27=M(:,177);X_30=M(:,180);
    Y_27=M(:,245);Y_30=M(:,248);
    Z_27=M(:,313);Z_30=M(:,316);
    
    %口
    X_61=M(:,211);X_63=M(:,213);
    Y_61=M(:,279);Y_63=M(:,281);
    Z_61=M(:,347);Z_63=M(:,349);
    
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
    
    
    k = 1;
    q = 0;
   
   
    for k = 1:20
   nodsum = 0;
   i = 1;

    while i < (renum-1)
        if ((fvangle(i,1) <= pthreshold - q) && (fvangle(i+1,1) > pthreshold - q))
            for j = 1:(renum-1)
                if ((i < j) && (j < i + 30))
                    if ((fvangle(j,1) >= mthreshold - q) && (fvangle(j+1,1) < mthreshold - q))
                        nodsum = nodsum + 1;
                        i = j;
                        break
                    end
                end
            end
        end
        i = i + 1;
    end
    if loop_num == 1
        nodcount{1, k} = q;
    end
    q = q + 0.5;
    nodcount{loop_num + 1, k} = nodsum;
    
    end
    %csvwrite(result_file, name, loop_num, 0);


 


    prompt = "続ける場合は1を，やめる場合はそれ以外を入力";
    decide_value = input(prompt);


end

alldata = [filename, nodcount];

writecell(alldata, result_file)

%{

csvwrite(result_file, nodcount(1,:), 0, 1)

csvnum = 0;

while csvnum < loop_num
    csvwrite(result_file, nodcount(csvnum + 1, :), loop_num, 1);

    csvnum = csvnum + 1;
end

%}


%{
%%ディレクトリ移動
move_dir = strcat('nod_', name);
mkdir(move_dir);
movefile(strcat(name,'*'),move_dir);

%}

disp('動作解析終了！');
