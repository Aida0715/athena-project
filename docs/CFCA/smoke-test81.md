# CFCA初回動作確認：test81ベース

## コードとビルド
- ソースのコミット：83d167c
- ビルドスクリプト：scripts/build/cfca.sh
- 環境：PrgEnv-oneapi/2024.2
- コンパイラ：Intel icpx 2024.2.1
- MPI：Intel MPI 2021.13
- ビルドノード：cna100（システムM）
- ビルドジョブID：335540
- ビルド終了コード：0
- 構成：Toyouchi、cartesian、isothermal、hlld、MHD、自己重力Multigrid
- 倍精度、NGHOST=4、MPI有効、OpenMP無効
- 最適化に-xhostを使用

## 計算条件
- 元input：inputs/hydro/Toyouchi_test/athinput.Toyouchi_71
- 動作確認input：inputs/hydro/Toyouchi_CFCA/athinput.cfca-smoke-test81
- 元inputからの変更：nlim=-1 → 10、SMR level=5 → 1
- MHD有効のビルドで初期磁場を0に設定
- ジョブスクリプト：run/CFCA/smoke-test81.sh
- パーティション：M-test-md
- 1ノード、MPI 4プロセス、1プロセス1コア、メモリ16GB
- 制限時間：10分

## 実行結果
- 実行日：2026-11-05
- ジョブID：335545
- Slurm状態：COMPLETED
- 終了コード：0:0
- ジョブ経過時間：11秒
- cycle=10、time=5.7393339835057784
- Terminating on cycle limitで正常終了
- Slurm標準エラーログは空
- 結果保存先：/work/aidarj/results/cfca-smoke81-335545/
- 実行ファイル・inputのコピーとSHA256は結果保存先に記録

## 注意・未確認事項
- Multigridの最粗格子が5×5×5となる効率上の警告あり
- 小規模な起動・並列実行確認であり、科学的な結果の検証ではない
- Novaとの数値結果の比較、長時間安定性、再スタートは未確認
- 複数ノードでのMPI通信とシステムPでの動作は未確認
