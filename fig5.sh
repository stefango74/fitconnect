#!/bin/sh

./src/TestReconstruct2D fig5a
./src/TestReconstruct2D fig5b
./src/TestReconstruct2D fig5c
./crop.sh fig_blend_a.png
./crop.sh fig_blend_b.png
./crop.sh fig_blend_c.png
