#!/bin/sh

./src/TestReconstruct2D fig6a
./src/TestReconstruct2D fig6b
./src/TestReconstruct2D fig6c
./src/TestReconstruct2D fig6d
./src/TestReconstruct2D fig6e
./src/TestReconstruct2D fig6f
./src/TestReconstruct2D fig6g
./src/TestReconstruct2D fig6h
./crop.sh fig_realdata_a.png
./crop.sh fig_realdata_b.png
./crop.sh fig_realdata_c.png
./crop.sh fig_realdata_d.png
./crop.sh fig_realdata_e.png
./crop.sh fig_realdata_f.png
./crop.sh fig_mouse_a.png
./crop.sh fig_mouse_b.png
