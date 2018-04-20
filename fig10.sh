#!/bin/sh

./src/TestReconstruct2D fig10a
./src/TestReconstruct2D fig10b
./src/TestReconstruct2D fig10c
./crop.sh fig_comparison2_a.png
./crop.sh fig_comparison2_b.png
./crop.sh fig_highnoise.png
