#!/bin/bash

mkdir -p bin
cd src/vdw
g++ -O3 hbs_main.cpp -o ../../bin/hbonds
