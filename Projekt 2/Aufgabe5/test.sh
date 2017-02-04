sh compileAll.sh
printf "\nGPU only:\n"
time ./cg
printf "\n\nMixed:\n"
time ./cg_mix
