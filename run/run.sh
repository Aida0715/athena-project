#!/bin/bash

# commit: 6ad162c //本シミュレーションに対応するjeans.cppをgitの履歴から追跡可
# AMR: add log-gradient condition

cd "$(dirname "$0")/.."

OUTDIR=results/Toyouchi-test6
mkdir -p $OUTDIR
cd $OUTDIR

# 実行（ファイル名はinputに任せる）
# 環境作成中のテスト計算用入力ファイルはToyouchi_testディレクトリに格納
../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi
