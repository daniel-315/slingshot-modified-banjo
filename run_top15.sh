#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "========================================"
echo " Slingshot + Modified Banjo: Top 15"
echo "========================================"

if ! command -v java >/dev/null 2>&1; then
    echo "ERROR: Java was not found."
    echo "Install it with:"
    echo "  sudo apt install -y openjdk-17-jdk"
    exit 1
fi

if [ ! -x /usr/bin/dot ]; then
    echo "ERROR: Graphviz was not found."
    echo "Install it with:"
    echo "  sudo apt install -y graphviz"
    exit 1
fi

if [ ! -f modified_banjo/build/edu/duke/cs/banjo/application/Banjo.class ]; then
    echo "Compiled Banjo was not found."
    echo "Compiling it now..."
    ./compile_banjo.sh
fi

mkdir -p results/reproduced/top15

echo
echo "Input: data/banjo_inputs/top15_slingshot_input.txt"
echo "Genes: 15"
echo "Observations: 959"
echo "Maximum networks requested: 5"
echo "Search time: 5 minutes"
echo
echo "Starting Banjo..."
echo

java \
    -cp modified_banjo/build \
    edu.duke.cs.banjo.application.Banjo \
    settingsFile=configs/banjo_top15.txt

echo
echo "Top-15 analysis complete."
echo "Results are in:"
echo "  results/reproduced/top15/"
