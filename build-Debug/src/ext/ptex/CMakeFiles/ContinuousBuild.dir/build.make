# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zcx/Master/Programmings/pbrt-v3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zcx/Master/Programmings/pbrt-v3/build-Debug

# Utility rule file for ContinuousBuild.

# Include any custom commands dependencies for this target.
include src/ext/ptex/CMakeFiles/ContinuousBuild.dir/compiler_depend.make

# Include the progress variables for this target.
include src/ext/ptex/CMakeFiles/ContinuousBuild.dir/progress.make

src/ext/ptex/CMakeFiles/ContinuousBuild:
	cd /home/zcx/Master/Programmings/pbrt-v3/build-Debug/src/ext/ptex && /usr/bin/ctest -D ContinuousBuild

ContinuousBuild: src/ext/ptex/CMakeFiles/ContinuousBuild
ContinuousBuild: src/ext/ptex/CMakeFiles/ContinuousBuild.dir/build.make
.PHONY : ContinuousBuild

# Rule to build all files generated by this target.
src/ext/ptex/CMakeFiles/ContinuousBuild.dir/build: ContinuousBuild
.PHONY : src/ext/ptex/CMakeFiles/ContinuousBuild.dir/build

src/ext/ptex/CMakeFiles/ContinuousBuild.dir/clean:
	cd /home/zcx/Master/Programmings/pbrt-v3/build-Debug/src/ext/ptex && $(CMAKE_COMMAND) -P CMakeFiles/ContinuousBuild.dir/cmake_clean.cmake
.PHONY : src/ext/ptex/CMakeFiles/ContinuousBuild.dir/clean

src/ext/ptex/CMakeFiles/ContinuousBuild.dir/depend:
	cd /home/zcx/Master/Programmings/pbrt-v3/build-Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zcx/Master/Programmings/pbrt-v3 /home/zcx/Master/Programmings/pbrt-v3/src/ext/ptex /home/zcx/Master/Programmings/pbrt-v3/build-Debug /home/zcx/Master/Programmings/pbrt-v3/build-Debug/src/ext/ptex /home/zcx/Master/Programmings/pbrt-v3/build-Debug/src/ext/ptex/CMakeFiles/ContinuousBuild.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/ext/ptex/CMakeFiles/ContinuousBuild.dir/depend

