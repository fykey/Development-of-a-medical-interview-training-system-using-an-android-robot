%OpenFaceのcsv429,430行目にvolume.csvの2行をコピー&ペースト、上書き保存
[file,path] = uigetfile('*.csv');
csv_file = fullfile(path,file);
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
