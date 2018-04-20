#!/bin/sh

./src/TestReconstruct2D fig9a
./src/TestReconstruct2D fig9b
./src/TestReconstruct2D fig9c
./src/TestReconstruct2D fig9d
./crop.sh fig_comparison_a.png
./crop.sh fig_comparison_b.png
./crop.sh fig_comparison_c.png
./crop.sh fig_comparison_d.png
