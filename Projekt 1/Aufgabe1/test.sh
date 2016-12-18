echo static
export OMP_SCHEDULE=static,1
./mandelbrot_ompp
export OMP_SCHEDULE=static,10
./mandelbrot_ompp
export OMP_SCHEDULE=static,100
./mandelbrot_ompp
export OMP_SCHEDULE=static,1000
./mandelbrot_ompp
echo dynamic
export OMP_SCHEDULE=dynamic,1
./mandelbrot_ompp
export OMP_SCHEDULE=dynamic,10
./mandelbrot_ompp
export OMP_SCHEDULE=dynamic,100
./mandelbrot_ompp
export OMP_SCHEDULE=dynamic,1000
./mandelbrot_ompp
echo guided
export OMP_SCHEDULE=guided,1
./mandelbrot_ompp
export OMP_SCHEDULE=guided,10
./mandelbrot_ompp
export OMP_SCHEDULE=guided,100
./mandelbrot_ompp
export OMP_SCHEDULE=guided,1000
./mandelbrot_ompp
echo auto
export OMP_SCHEDULE=auto,1
./mandelbrot_ompp
export OMP_SCHEDULE=auto,10
./mandelbrot_ompp
export OMP_SCHEDULE=auto,100
./mandelbrot_ompp
export OMP_SCHEDULE=auto,1000
./mandelbrot_ompp

