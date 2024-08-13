#!/usr/bin/env bash

terminal="$(tty)"
columns=$(stty -a <"$terminal" | grep -Po '(?<=columns )\d+')
rows=$(stty -a <"$terminal" | grep -Po '(?<=rows )\d+')

echo "cols:$columns"
echo "rows:$rows"
