#!/bin/sh

./src/TestReconstruct2D fig1a
./src/TestReconstruct2D fig1b
./src/TestReconstruct2D fig1c
./src/TestReconstruct2D fig1d
./crop.sh fig_teaser_a.png
./crop.sh fig_teaser_b.png
./crop.sh fig_teaser_c.png
./crop.sh fig_teaser_d.png
