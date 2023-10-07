#!/bin/bash

echo "Compile 'main.c' ..."
gcc main.c utils.c -o GaussMethod -Wall -pedantic 
echo "Compile 'main_OpenMP.c' ..."
gcc main_OpenMP.c utils.c -o GaussMethodOpenMP -Wall -pedantic -fopenmp
echo "Compilation is complete!"
