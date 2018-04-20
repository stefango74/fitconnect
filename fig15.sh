#!/bin/sh

./src/TestReconstruct2D fig15a
./src/TestReconstruct2D fig15b
./crop.sh fig_limit_a.png
./crop.sh fig_limit_b.png
