#!/bin/sh

./src/TestReconstruct2D fig14a
./src/TestReconstruct2D fig14b
./crop.sh fig_feature_a.png
./crop.sh fig_feature_b.png
