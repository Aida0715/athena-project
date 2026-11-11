# vcap77：速度上限100・300・1000 km/sの20 kyr比較

- 基準input：inputs/hydro/Toyouchi_test/athinput.Toyouchi_77_cp
- 3つのinputは空白整形を除き、velocity_ceiling_kmsだけを変更。
- tlim=20、nlim=-1、SMR level=5、fM=1、初期磁場0。
- VTK/HST間隔は基準input通り0.2 kyr（旧model71の2 kyrではない）。
- velocity_ceiling_diagnostics=true。各実行のToyouchi.velocity_ceiling.csvを比較する。
- 1モデルにつきM-bulk-md、1ノード、MPI 32、1プロセス1コア、32 GiB。
- Slurm 12時間、Athena++ 11時間45分。20 kyr未到達なら手動再開が必要。

## 旧テストの廃止

model71用input・runスクリプト・説明書を削除。
CFCAで `bash run/CFCA/cleanup-model71.sh` を実行すると、ユーザー自身の
cfca-model71ジョブを停止し、停止確認後に次だけを削除する。

- /work/$USER/results/cfca-model71-<数字>/
- リポジトリ直下のcfca-model71-<数字>.out/.err

smoke81・restart81・benchmark71は検証記録および既存ビルドの由来として保持。
Novaや手元PCへコピーした解析用データは、このスクリプトの削除対象外。

## CFCAでの投入

リポジトリ直下でoneAPIを有効化して実行する。

```bash
source /work/opt/local/bin/enable-oneapi.sh
bash run/CFCA/cleanup-model71.sh
bash run/CFCA/submit-vcap77.sh
```

submit-vcap77.shはコミット済みソースをgit archiveで独立した作業領域に固定し、
inputと投入スクリプトもコピーする。計算ノードで新たにIntel MPI版をビルドし、
成功を条件として3要素のジョブ配列を投入する。Novaの実行ファイルは使わない。

- 配列0：100 km/s
- 配列1：300 km/s
- 配列2：1000 km/s

3モデルを一度に投入するが、同時実行数はQOSと空き資源に従う。
保存先はsubmitスクリプトが表示するBundle以下。submission.txtにジョブID、
source-commit.txtにソース版、build-info.txtにビルド環境、各run-info.txtに
実行条件・終了状態・SHA256を記録する。現時点で実際の投入IDは未記録。

## 再開

時間制限で止まったモデルだけを、同じBundleのスクリプトで再投入する。
例えば300 km/s（配列1）は以下。新しい結果ディレクトリへ保存される。

```bash
sbatch --array=1 --output="$BUNDLE/logs/run-%A_%a.out" \
  --error="$BUNDLE/logs/run-%A_%a.err" \
  "$BUNDLE/vcap77.sh" "$BUNDLE" "$PREVIOUS_RESULT"
```

BUNDLEとPREVIOUS_RESULTには実際の絶対パスを設定する。
再開時は実行ファイルとinputが前回と同一であることを検証する。

## 次回の解析

全モデルが20 kyrに到達したかを確認し、同じ物理時刻で比較する。
降着率・密度構造・速度分布に加え、CSVの補正質量、除去エネルギー、
除去角運動量を確認する。補正質量の累積は同じガスを重複して数える。
速度制限は各ステップ末に適用されるため、CSVのmax_before_kmsは上限を超え得る。
