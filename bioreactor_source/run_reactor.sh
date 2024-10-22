file="reactor.c"

# For the moment his has to be set to 1 otw the code breaks
THREADS=1

# Loop over paramter to allow error study more easily (feel free to change this)

for VISCOSITY in 0.0015 0.003; do

  # Create filename based on parameter of intrest
  cp -r reactor_master/ Reactor-Visc$VISCOSITY
	cd Reactor-Visc$VISCOSITY

  # Compile Basilisk including MPI
  CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 "$file" -o "out_reactor" -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11

  echo "Running code"
	
  # Parameters in SI units: liquid_density gas_density liquid_viscosiy gas_viscosity surface_tension shaking_frequency reactor_width acceleration_gravity level max_time 
  mpirun -np $THREADS ./out_reactor 1000.0 1.0 $VISCOSITY 0.00001 0.07 0.5 0.05 9.81 7 80.0

  cd ..

done
