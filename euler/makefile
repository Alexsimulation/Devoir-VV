INCLUDES := -I .
LIBS := -lGL -lGLU -lGLEW
OPTIM := -O3 -fopenmp -fopenmp-simd -lgomp
STD := -std=c++17
GTKMM := `pkg-config --cflags --libs gtkmm-3.0`

build:
	@echo "Building gui app"
	@g++ src/gui.cpp -o bin/euler -w ${STD} ${GTKMM} ${LIBS} ${INCLUDES} ${OPTIM}

profile:
	@echo "Building gui profiled app"
	@g++ src/gui.cpp -o bin/euler_prof -w ${STD} ${GTKMM} ${LIBS} ${INCLUDES} ${OPTIM} -g
