# Athena++ Project(修士研究用)

## 1. プロジェクト概要
このプロジェクトは、Athena++ を用いた自己重力流体の数値シミュレーションを管理・再現するためのリポジトリです。
主な目的は以下の通りです：
- 各実験の再現性確保
- コード変更履歴管理 (git commit)
- 実験条件、入力ファイル、結果の体系的管理
- Python/Jupyter 解析コードとの連携

対象例：
- AMR（Adaptive Mesh Refinement）の実装検証
- ジーンズ条件下のガス分布解析
- 中心星ポテンシャル下の降着解析

## 2. プロジェクト構成
athena-project/
├── README.md # このファイル (git管理内)
├── docs/ # 個別実験の詳細記録 (git管理内)
├── src/ # Athena++ のソースコード (git管理内)
├── inputs/ # 各実験の入力ファイル (git管理内)
├── run/ # 実行用スクリプト (run.sh (git管理内))
├── results/ # 計算結果ディレクトリ（git管理外）
├── analysis/ # Python/Jupyter解析用コード (git管理内)
├── Makefile # ビルド用
├── configure.py # Athena++ ビルド設定
└── .gitignore # 管理外ファイルの指定

## 3. 必要環境
- Ubuntu 22.04 以上
- Python 3.10 以上
- C++17 対応コンパイラ（g++ / clang++）
- MPI (必要に応じて)
- Athena++ ビルド済み環境
- git

## 4.実験管理方法
1. 各実験ごとに入力ファイルを inputs/ に作成
   例: athinput.Toyouchi
2. 実行用スクリプトを作成 run/run.sh
   ・ 使用バイナリ、入力ファイル、出力ディレクトリを指定
   ・ コードの commit ID をコメントに記載
3. 結果は results/ に保存
4. 実験記録は docs/ にまとめる
   ・ 使用 commit ID
   ・ 入力ファイル名
   ・ 実験内容・注意点

## 5.git管理ツール
   ・ src/, inputs/, run/, analysis/, docs/ を追跡
   ・ results/, bin/, obj/ は .gitignore に追加	
   ・ 各実験やコード変更ごとに commit 作成
   ・ commit メッセージは簡潔かつ意味がわかる内容にする
   　 例: AMR: add log-gradient condition

## 6.実験再現フロー
1. docs/ にある対象実験ファイルを確認
2. 使用 commit ID を確認
   Bash: git checkout <commit-id>
3. 対応する run/run.sh を実行
   Bash:./run/run.sh
4. results/ で計算結果を確認

## 7.参考文献・先行研究
・先行研究(Toyouchi+23): https://arxiv.org/pdf/2206.14459
・Athena++ Documentation: https://github.com/PrincetonUniversity/Athena

