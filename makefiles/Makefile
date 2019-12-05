hn = $(shell hostname)
username = $(shell whoami)

#CC=$(shell which g++)
#CC=/usr/local/Cellar/gcc/8.2.0/bin/g++-8
CPPSETTINGS+= -std=c++11
CPPFLAGS+= -O3 -Wno-all -Wno-deprecated -w -Wno-reorder #-march=native

COMMON_OBJS = \
	main.o \
	GeometryFunctions.o \
	Polymer.o \
	PolymerIntegrator.o \
	PositionVerlet2nd.o \
	Rod.o \
	SplineProfileZeroEnds.o \
	RodInitialConfigurations.o \
	MRAGProfiler.o \
	Vector3.o \
	Matrix3.o \
	SpeedFunctions.o \
	MathFunctions.o

VPATH := ../source/
.DEFAULT_GOAL := elbow

#Linking stage
elbow: CPPSETTINGS += -DFLAGELBOW
elbow: $(COMMON_OBJS) Elbow.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
flagella: CPPSETTINGS += -DFLAGFLAGELLA
flagella: $(COMMON_OBJS) Flagella.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
walker: CPPSETTINGS += -DFLAGWALKER
walker: $(COMMON_OBJS) Walker.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
muscularsnake: CPPSETTINGS += -DFLAGMUSCULARSNAKE
muscularsnake: $(COMMON_OBJS) MuscularSnake.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
wing: CPPSETTINGS += -DFLAGWING
wing: $(COMMON_OBJS) Wing.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
sphericaljoint: CPPSETTINGS += -DFLAGSPHERICALJOINT
sphericaljoint: $(COMMON_OBJS) SphericalJoint.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
hingejoint: CPPSETTINGS += -DFLAGHINGEJOINT
hingejoint: $(COMMON_OBJS) HingeJoint.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
fixedjoint: CPPSETTINGS += -DFLAGFIXEDJOINT
fixedjoint: $(COMMON_OBJS) FixedJoint.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
pullingmuscle: CPPSETTINGS += -DFLAGPULLINGMUSCLE
pullingmuscle: $(COMMON_OBJS) PullingMuscle.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
timoshenkobeam: CPPSETTINGS += -DFLAGQUASISTATICTIMOSHENKOBEAM
timoshenkobeam: $(COMMON_OBJS) QuasistaticTimoshenkoBeam.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
helicalbuckling: CPPSETTINGS += -DFLAGHELICALBUCKLING
helicalbuckling: $(COMMON_OBJS) InstabilityHelical.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Linking stage
snake: CPPSETTINGS += -DFLAGSNAKE
snake: $(COMMON_OBJS) Snake.o
	$(CC) $(CPPSETTINGS) $(LIB) $^ -o $@
	@echo done

#Compiling stage
%.o: %.cpp
	$(CC) $(CPPSETTINGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o
	rm -f *.s
	rm -f timoshenkobeam
	rm -f helicalbuckling
	rm -f snake
	rm -f sphericaljoint
	rm -f hingejoint
	rm -f fixedjoint
	rm -f pullingmuscle
	rm -f elbow
	rm -f flagella
	rm -f walker
	rm -f muscularsnake
	rm -f wing
