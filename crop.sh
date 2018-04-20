#!/bin/sh

convert $1 -trim +repage $1.trimmed
mv $1.trimmed $1
