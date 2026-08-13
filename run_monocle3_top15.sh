#!/usr/bin/env bash
set -e

cd "$(dirname "$0")"

if [ ! -f modified_banjo/build/edu/duke/cs/banjo/application/Banjo.class ]; then
    echo "Banjo has not been compiled. Run ./compile_banjo.sh first."
    exit 1
fi

rm -rf results/reproduced/monocle3/top15
mkdir -p results/reproduced/monocle3/top15

java \
  -cp modified_banjo/build \
  edu.duke.cs.banjo.application.Banjo \
  settingsFile=configs/banjo_monocle3_top15.txt
