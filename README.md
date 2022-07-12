# Development-of-a-medical-interview-training-system-using-an-android-robot
研究内容．人工筋肉で実際に表情が動くロボットを使った，医学部生のための医療面接訓練の開発を行っている．

![S__49709099](https://user-images.githubusercontent.com/73433285/178558213-45778046-a477-41ac-883f-d15cf3587e5d.jpg)

診察動画を取って，OpenFaceのFeatureExtractionのアプリで，動画上の特徴点の変位や，そこから算出されるパラメータを取得．
そのパラメータをつかって，うなずき回数，アイコンタクト，発話区間，表情などを自動で計測するプログラムの開発を行っています．
また，ロボットと医学部生が対話をするためのプログラムも作成している．

![image](https://user-images.githubusercontent.com/73433285/178556497-783a2883-3dbf-48cc-93f8-dca851bf0bef.png)

以下，プログラムの説明．
・MATLAB動作解析プログラム
Motion_analysis2021.m
撮影した動画を解析し，非言語動作の評価を出すプログラム．
以下のプログラムから構成されている．

①	OpenFace_runner.m
MATLABでOpenFaceを実行し，動画解析を行うプログラム．

②	Voice_analysis2021.m
動画の音声解析を行うプログラム．

③	Csv_copypaste.m
OpenFaceで出力されたCSVファイルと音声解析で出力されたCSVファイルを結合するプログラム．

④	OpenFace_runner_picture.m
医師が患者に対してアイコンタクトを取っている視線を撮影した画像からOpenFaceにてその視線角度を調べるプログラム．

⑤	Nod_analysis22021.m
頷き回数・アイコンタクト量・表情変化の算出を行うプログラム．
以下のプログラムは頷き・アイコンタクト・表情を別々に検出できるプログラムになる．

A)	Nod_analysis_nod_detection.m
頷き回数を算出するプログラム．

B)	Nod_analysis_eye_contact.m
アイコンタクト量を測るプログラム．

C)	Nod_analysis_facial_expression.m
表情変化を調べるプログラム．

⑥	Move_file.m
撮影した動画とその解析結果を1つのフォルダにまとめるプログラム．
