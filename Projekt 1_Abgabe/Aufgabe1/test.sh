echo static
export OMP_SCHEDULE=static,1
./mandelbrot
export OMP_SCHEDULE=static,10
./mandelbrot
export OMP_SCHEDULE=static,100
./mandelbrot
export OMP_SCHEDULE=static,1000
./mandelbrot
echo dynamic
export OMP_SCHEDULE=dynamic,1
./mandelbrot
export OMP_SCHEDULE=dynamic,10
./mandelbrot
export OMP_SCHEDULE=dynamic,100
./mandelbrot
export OMP_SCHEDULE=dynamic,1000
./mandelbrot

