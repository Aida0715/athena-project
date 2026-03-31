#!/bin/bash

# commit: 6ad162c //本シミュレーションに対応するjeans.cppをgitの履歴から追跡可
# AMR: add log-gradient condition

cd "$(dirname "$0")/.."

OUTDIR=results/Toyouchi-test6
mkdir -p $OUTDIR
cd $OUTDIR

# 実行（ファイル名はinputに任せる）
../../bin/athena -i ../../inputs/athinput.Toyouchi
