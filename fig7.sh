#!/bin/sh

./src/TestReconstruct2D fig7a
./src/TestReconstruct2D fig7b
./src/TestReconstruct2D fig7c
./src/TestReconstruct2D fig7d
./crop.sh fig_basic_a.png
./crop.sh fig_basic_b.png
./crop.sh fig_basic_c.png
./crop.sh fig_basic_d.png
