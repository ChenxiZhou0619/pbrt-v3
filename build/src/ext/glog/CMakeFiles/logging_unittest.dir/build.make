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
CMAKE_BINARY_DIR = /home/zcx/Master/Programmings/pbrt-v3/build

# Include any dependencies generated for this target.
include src/ext/glog/CMakeFiles/logging_unittest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/ext/glog/CMakeFiles/logging_unittest.dir/compiler_depend.make

# Include the progress variables for this target.
include src/ext/glog/CMakeFiles/logging_unittest.dir/progress.make

# Include the compile flags for this target's objects.
include src/ext/glog/CMakeFiles/logging_unittest.dir/flags.make

src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o: src/ext/glog/CMakeFiles/logging_unittest.dir/flags.make
src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o: /home/zcx/Master/Programmings/pbrt-v3/src/ext/glog/src/logging_unittest.cc
src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o: src/ext/glog/CMakeFiles/logging_unittest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zcx/Master/Programmings/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o"
	cd /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o -MF CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o.d -o CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o -c /home/zcx/Master/Programmings/pbrt-v3/src/ext/glog/src/logging_unittest.cc

src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.i"
	cd /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zcx/Master/Programmings/pbrt-v3/src/ext/glog/src/logging_unittest.cc > CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.i

src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.s"
	cd /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zcx/Master/Programmings/pbrt-v3/src/ext/glog/src/logging_unittest.cc -o CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.s

# Object files for target logging_unittest
logging_unittest_OBJECTS = \
"CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o"

# External object files for target logging_unittest
logging_unittest_EXTERNAL_OBJECTS =

src/ext/glog/logging_unittest: src/ext/glog/CMakeFiles/logging_unittest.dir/src/logging_unittest.cc.o
src/ext/glog/logging_unittest: src/ext/glog/CMakeFiles/logging_unittest.dir/build.make
src/ext/glog/logging_unittest: src/ext/glog/libglog.a
src/ext/glog/logging_unittest: /usr/lib/libunwind.so
src/ext/glog/logging_unittest: src/ext/glog/CMakeFiles/logging_unittest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zcx/Master/Programmings/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable logging_unittest"
	cd /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/logging_unittest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/ext/glog/CMakeFiles/logging_unittest.dir/build: src/ext/glog/logging_unittest
.PHONY : src/ext/glog/CMakeFiles/logging_unittest.dir/build

src/ext/glog/CMakeFiles/logging_unittest.dir/clean:
	cd /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog && $(CMAKE_COMMAND) -P CMakeFiles/logging_unittest.dir/cmake_clean.cmake
.PHONY : src/ext/glog/CMakeFiles/logging_unittest.dir/clean

src/ext/glog/CMakeFiles/logging_unittest.dir/depend:
	cd /home/zcx/Master/Programmings/pbrt-v3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zcx/Master/Programmings/pbrt-v3 /home/zcx/Master/Programmings/pbrt-v3/src/ext/glog /home/zcx/Master/Programmings/pbrt-v3/build /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog /home/zcx/Master/Programmings/pbrt-v3/build/src/ext/glog/CMakeFiles/logging_unittest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/ext/glog/CMakeFiles/logging_unittest.dir/depend

