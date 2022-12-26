%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEX SCRIPT
% 
% This script will compike all the mex files needed for the NAOMI simulator
%
% 2017 - Adam Charles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ./MEX/
fprintf('Compiling array_SubMod.cpp and array_SubModTest.cpp for faster array substitution...\n')
mex -largeArrayDims array_SubMod.cpp '-compatibleArrayDims'
mex -largeArrayDims array_SubModTest.cpp '-compatibleArrayDims'
fprintf('Compiling array_SubSub.cpp and array_SubSubTest.cpp for faster array summing...\n')
mex -largeArrayDims array_SubSub.cpp '-compatibleArrayDims'
mex -largeArrayDims array_SubSubTest.cpp '-compatibleArrayDims'
fprintf('Compiling dendrite_dijkstra_cpp.cpp for faster dendrite growth...\n')
mex -largeArrayDims dendrite_dijkstra_cpp.cpp '-compatibleArrayDims'
fprintf('Compiling dendrite_randomwalk_cpp.cpp for faster dendrite growth...\n')
mex -largeArrayDims dendrite_randomwalk_cpp.cpp '-compatibleArrayDims'
cd ..
fprintf('done.\n')