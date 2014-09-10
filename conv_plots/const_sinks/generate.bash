#!/bin/bash

cd mean_conv
gnuplot mean_conv.plt
gnuplot mean_conv_fit.plt
cd ..
cd cov_conv
gnuplot cov_conv.plt
gnuplot cov_conv_fit.plt
cd ..
latexmk -g const_sinks
