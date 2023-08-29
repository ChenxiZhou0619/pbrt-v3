# Install script for directory: /home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/zcx/Master/Programmings/pbrt-v3/build/src/ext/openexr/IlmBase/IlmThread/libIlmThread.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThreadPool.h"
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThread.h"
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThreadSemaphore.h"
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThreadMutex.h"
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThreadNamespace.h"
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThreadExport.h"
    "/home/zcx/Master/Programmings/pbrt-v3/src/ext/openexr/IlmBase/IlmThread/IlmThreadForward.h"
    )
endif()

