# model71：本計算前の運用・設計確認用テスト

## 計算条件
- 目標：1 Myr（tlim=1000、nlim=-1）
- 元input：inputs/hydro/Toyouchi_test/athinput.Toyouchi_71
- 使用input：inputs/hydro/Toyouchi_CFCA/athinput.cfca-model71
- SMR level=5、538 MeshBlocks、初期磁場0、fM=1.0
- VTK：2000年間隔、HST：200年間隔
- 定期restart：約10.032 kyr間隔
- 実行ファイル：ジョブ335545で保存したもの（ソース83d167c）
- スクリプト：run/CFCA/model71.sh
- 長時間実行設定：M-bulk-md、1ノード、MPI 32プロセス、32 GiB
- Slurm制限12時間、Athena++制限11時間45分
- 各ジョブの結果は別ディレクトリに保存。再投入は手動

## 時間制限・再開確認（2026-11-09）
M-test-md、Slurm制限5分、Athena++制限2分で確認。

- 初回ジョブ338526：
  cycle=311、time=1.4784262313612877で時間制限終了
- 再開ジョブ338530：
  上記状態から再開し、
  cycle=630、time=2.0085491554103423で時間制限終了
- 両ジョブともCOMPLETED、終了コード0:0、標準エラーなし
- 両ジョブでToyouchi.final.rst（約121 MiB）を保存
- 結果：/work/aidarj/results/cfca-model71-338526/
- 結果：/work/aidarj/results/cfca-model71-338530/

## 状況・注意
- 1 Myrには未到達。次回は338530の保存状態から継続可能
- 再開時は前回の実行ファイルとrestartを使用
- athinputは条件の記録用。再開時の設定はrestart内から読み込む
- 解析には各ジョブの出力を時系列で合わせて使用する
- 連続計算との数値一致と長時間安定性は未確認
- 保存領域の個人容量・ファイル数上限は未確認
