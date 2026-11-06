# Nova・CFCAの小規模計算比較

- 比較日：2026-11-06
- 両環境のソース：83d167c
- Nova：g++ / Open MPI
- CFCA：icpx 2024.2.1 / Intel MPI 2021.13
- 両環境：MPI 4プロセス、OpenMP無効、倍精度
- input：test81で使用するathinput.Toyouchi_71を基に、
  nlim=10、SMR level=1へ変更
- 実行に使用した両inputの差は末尾空白と空行のみ
- 両環境ともcycle=10で正常終了

## 比較結果

- Toyouchi.hst：ファイル全体が一致（出力精度内）
- cycle=0〜10のログを比較
- time最大絶対差：1.651124e-12（コード時間単位）
- time最大相対差：4.623547e-13
- dt最大絶対差：5.397904e-13（コード時間単位）
- dt最大相対差：8.721719e-13
- この小規模テストの全体量・時刻・時間刻みは良好に一致

## データ保存先

- Nova：/home/aida/athena-cfca-compare-83d167c/results-smoke81-oIlN4A/
- CFCA：/work/aidarj/results/cfca-smoke81-335545/
- CFCAから取得した比較資料：
  /home/aida/athena-cfca-compare-83d167c/cfca-compare-335545/

## 検証の範囲

- 格子ごとの物理量、長時間安定性、再スタートは未確認
- 初期磁場は0であり、非ゼロ磁場での検証ではない
- 複数ノードでのMPI通信は未確認
- Novaは別計算と同時に低優先度で実行したため、性能比較はしない
