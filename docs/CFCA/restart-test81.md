# CFCA再スタート動作確認

- 実施日：2026-11-06
- 元ジョブ：335545
- 元データ：/work/aidarj/results/cfca-smoke81-335545/
- 元ジョブで保存した実行ファイルとToyouchi.final.rstを使用
- ジョブスクリプト：run/CFCA/restart-test81.sh
- 実行：srun ./athena -r restart-source.rst time/nlim=20
- M-test-md、1ノード、MPI 4プロセス、OpenMP無効

## 結果

- ジョブID：336403
- Slurm：COMPLETED、終了コード0:0
- 経過時間：9秒
- 再開：cycle=10、time=5.7393339835057784
- 終了：cycle=20、time=9.4736242947174354
- Terminating on cycle limitで正常終了
- Slurm標準エラーは空
- 保存先：/work/aidarj/results/cfca-restart81-336403/
- 使用した実行ファイル・再スタートファイルのSHA256はrun-info.txtに記録

## 検証範囲

- 同一実行ファイル・同一MPIプロセス数での再開動作を確認
- 連続20ステップ計算との数値一致は未確認
- 本計算の格子規模・非ゼロ磁場・複数ノードでの再開は未確認
