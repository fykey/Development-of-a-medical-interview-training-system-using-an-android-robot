%新しくフォルダを作成し，撮影した動画，解析結果，グラフを分かりやすいようにまとめる
[file,path] = uigetfile('*.csv');
data_file = fullfile(path,file);
in_file=data_file;
[~,name,~] = fileparts(in_file);

mkdir(name);
movefile(strcat(name,'*'),name);

disp('移動完了！');
