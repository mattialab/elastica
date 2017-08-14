hn ?= $(shell hostname)
username ?= $(shell whoami)

ifeq "$(hn)" "monolith.mechse.illinois.edu"
include make.monolithmechse
endif

ifeq "$(hn)" "Mattias-MacBook-Pro.local"
include make.mattia
endif

ifeq "$(hn)" "Wyss11139.local"
include make.monolith
endif

ifeq "$(hn)" "monolith.seas.harvard.edu"
include make.monolith
endif

ifeq "$(hn)" "Raghunaths-MacBook-Pro.local"
include make.ragu
endif

ifeq "$(hn)" "Levis-MacBook-Pro.local"
include make.levi
endif

ifeq "$(hn)" "levis-mbp"
include make.levi
endif

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
endif 

ifeq "$(viz)" "true"
CPPFLAGS += -DSNAKE_VIZ
endif 

ifeq "$(andrew)" "true"
CPPFLAGS += -DANDREW
endif  

SNAKE_OBJS = \
	main.o \
	InterfaceCma.o \
	GeometryFunctions.o \
	Polymer.o \
	PolymerIntegrator.o \
	VelocityVerlet2nd.o \
	PositionVerlet2nd.o \
	PositionVerletGeneral.o \
	AdaptiveVelocityVerlet.o \
	Rod.o \
	SplineProfileZeroEnds.o \
	RodInitialConfigurations.o \
	MRAGProfiler.o \
	Vector3.o \
	Matrix3.o \
	SpeedFunctions.o \
	MathFunctions.o \
	EulerBuckling.o \
	EulerBeam.o \
	SolenoidsJCP.o \
	Solenoids.o \
	Snake.o \
	StigmaticStart.o \
	MitchellBuckling.o \
	InstabilityHelical.o \
	MassSpringSystem.o \
	TorsionalWavesCouple.o \
	TorsionalWavesCoupleStretch.o \
	QuasistaticTimoshenkoBeam.o \
	RollingFrictionInitialVelocity.o \
	RollingFrictionInclinedPlane.o \
	RollingFrictionTorque.o \
	AxialFriction.o \
	SlenderBodyStokes.o

VPATH := ../source/
.DEFAULT_GOAL := snake

#Linking stage
snake: $(SNAKE_OBJS)
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Compiling stage	
%.o: %.cpp
	$(CC) $(CPPSETTINGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o
	rm -f *.s
	rm -f snake
