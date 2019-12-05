#!/usr/bin/env bash
# Not posix compliant

read -rd '' globalhelp <<-EOF
	usage
	-----
	./run.sh <case_name>
	
	case_names and explanations
	---------------------------
		help : Print this help message
	
		timoshenkobeam : Runs a cantilevered slender beam (clamped at the wall at one
		end, free at another) bending under the influence of a downward point force.
		We compare against analytical results using Timoshenko Beam theory.
	
		helicalbuckling : Runs the localized helical buckling benchmark which tests
		the mixed bending modes (two bending modes + twist mode). When an unstretchable
		rod is twisted at two ends, a characteristic instability arises which forms a
		localized helical shape at the center of the beam. There is an analytical solution
		(see Gazzola et. al, RSOS, 2018), which we compare against.
	
		snake : Runs an example of an snake, actuated by a continuum analytical bend-
		ing torque profile slithering on a plane-ground with anisotropic friction.
		This snake has been optimized for maximal forward speed.
	
		sphericaljoint : Runs an example of two rods connected by a spherical joint
		that allows for arbitrary rotations (see SphericalJoint.cpp)
	
		hingejoint : Runs an example of two rods connected by a hinge joint that
		allows for motion only along a constrained plane (see HingeJoint.cpp)
	
		fixedjoint : Runs an example of two rods connected by a fixed joint that
		does not allow any motion between the rods (see FixedJoint.cpp)
	
		pullingmuscle : Runs an example of two hanging vertical rods (connected to
		ground by a hinge) connected to each another's centers by a horizontal
		muscle. This muscle is actuated and pulls the rods together.
		(see PullingMuscle.cpp)
	
		elbow : Runs the injured elbow (50% Strength) with artificial muscles,
		as seen in Fig.7 in Supplementary Material.
	
		flagella : Runs the original biohybrid flagella design seen in Fig.2(a,b).
	
		walker : Runs the original biohybrid walker design at a frequency of 2 Hz,
		as seen in Fig.2(e,f). SLOW without parallelization.
	
		muscularsnake : Runs the muscular snake seen in Fig.3. SLOW without parallelization.
	
		wing : Runs the wing seen in Fig.4. EXTREMELY SLOW without parallelization.
	
	case_names
	----------
	timoshenkobeam, helicalbuckling, sphericaljoint, hingejoint, fixedjoint, pullingmuscle, snake
	flagella, muscularsnake, walker, elbow, wing
	
	All results are stored in folder 'run_<case_name>' in the main directory
EOF

if [[ $1 =~ ^([hH][eE][lL][pP]|[hH])$ || $# -eq 0 ]]; then
	echo "${globalhelp}"
	exit 0
fi

rm -fr "../run_$1"
mkdir -p "../run_$1"
if command -v make >/dev/null 2>&1; then

	# check gcc version starting from 9 on to 4
	CC_ver_arr=($(seq 9 -1 4))
	CC_ver_arr=("${CC_ver_arr[@]/#/g++-}")
	# Try and detect GNU g++ from the shell, if not use default CC
	for cc_ver in "${CC_ver_arr[@]}"; do
		if command -v "${cc_ver}" >/dev/null 2>&1 && "${cc_ver}" --version | grep -q '[Gg][Cc][Cc]'; then
			CC="${cc_ver}"
			break
		fi
	done
	# Check if not set, else set it
	if [ -z "${CC}" ] && g++ --version | grep -q '[Gg][Cc][Cc]'; then
		CC="g++"
	fi
	if [ -z "${CC}" ]; then
		# We have no hope, fall back on XZ's g++
		CC="/usr/local/Cellar/gcc/8.2.0/bin/g++-8"
	fi
	make clean
	make "$1" -j16 CC="${CC}"
else
	echo "Make is a required dependency, please install Make"
fi
cp ../makefiles/"$1" ../run_"$1"/executable
cd ../run_"$1" || {
	echo "cd failed, exiting now"
	exit 1
}

if [[ "$1" == "wing" ]]; then

	echo '#################################'
	echo ' WARNING : THIS CODE IS VERY SLOW'
	echo '#################################'

elif [[ "$1" == "elbow" ]]; then
	echo "test"
elif [[ "$1" == "timoshenkobeam" ]]; then

	echo '#################################################'
	echo ' WARNING : This simulates 5000 physical'
	echo ' seconds to reach equlibrium. This large'
	echo ' number is necessary because the damping rate nu'
	echo ' is very small'
	echo '#################################################'
fi

# Finally run the code
echo ""
echo '################################################'
echo "Running code now...results in run_$1 folder"
echo '################################################'
echo ""
./executable
