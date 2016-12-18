icc matrixmul_lin.c -o matrixmul_lin_N$2_$1_icc -g -Wall -$1 -D N=$2 -openmp
icc matrixmul_red.c -o matrixmul_red_N$2_$1_icc -g -Wall -$1 -D N=$2 -openmp
icc matrixmul_in.c -o matrixmul_in_N$2_$1_icc -g -Wall -$1 -D N=$2 -openmp
icc matrixmul_out.c -o matrixmul_out_N$2_$1_icc -g -Wall -$1 -D N=$2 -openmp
gcc matrixmul_lin.c -o matrixmul_lin_N$2_$1_icc -g -Wall -$1 -D N=$2 -fopenmp
gcc matrixmul_red.c -o matrixmul_red_N$2_$1_icc -g -Wall -$1 -D N=$2 -fopenmp
gcc matrixmul_in.c -o matrixmul_in_N$2_$1_icc -g -Wall -$1 -D N=$2 -fopenmp
gcc matrixmul_out.c -o matrixmul_out_N$2_$1_icc -g -Wall -$1 -D N=$2 -fopenmp
