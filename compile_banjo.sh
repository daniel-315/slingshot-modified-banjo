#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "========================================"
echo " Compiling Modified Banjo"
echo "========================================"

if ! command -v javac >/dev/null 2>&1; then
    echo "ERROR: Java compiler (javac) was not found."
    echo
    echo "On Ubuntu, install Java with:"
    echo "  sudo apt update"
    echo "  sudo apt install -y openjdk-17-jdk"
    exit 1
fi

rm -rf modified_banjo/build
mkdir -p modified_banjo/build

find modified_banjo/source \
    -name '*.java' \
    -print0 |
    xargs -0 javac \
        --release 8 \
        -d modified_banjo/build

echo
echo "Compilation successful."
echo "Compiled class files:"
find modified_banjo/build -name '*.class' | wc -l
