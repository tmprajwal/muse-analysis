# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)
# Needed for cmake 2.6
IF(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# Set project name
PROJECT(tracker)




# Set the library version
SET(tracker_MAJOR_VERSION 1)
SET(tracker_MINOR_VERSION 1)
SET(tracker_PATCH_VERSION 0)
SET(tracker_VERSION "${tracker_MAJOR_VERSION}.${tracker_MINOR_VERSION}.${tracker_PATCH_VERSION}")
LIST(APPEND tracker_LIBRARY_PROPERTIES VERSION "${tracker_VERSION}" SOVERSION "${tracker_MAJOR_VERSION}")


# find source and header files, set paths, includes...
SET(LIBNAME "tracker")


FILE(GLOB SRCS src/*.cpp)
FILE(GLOB HEADERS src/include/*.h ${CMAKE_SOURCE_DIR}/include/muserawtree.h)

SET(INCLUDES ${ROOT_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/src/include" "${CMAKE_SOURCE_DIR}/include" ${GENFIT_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${INCLUDES})

# create root dictionary
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/include/LinkDef.hh)
SET(DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS})
TARGET_LINK_LIBRARIES(${LIBNAME}  musetree Plugin "-L${ROOT_LIBRARY_DIR}" ${ROOT_LIBRARIES} Eve ${GENFIT_LIBRARIES})

#TARGET_COMPILE_FEATURES(${LIBNAME} PRIVATE cxx_auto_type cxx_range_for cxx_decltype_auto)
TARGET_COMPILE_FEATURES(${LIBNAME} PRIVATE cxx_auto_type cxx_range_for)
SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES ${tracker_LIBRARY_PROPERTIES})
#TARGET_COMPILE_FEATURES(${LIBNAME} PRIVATE cxx_trailing_return_types)

# install the library
INSTALL(TARGETS ${LIBNAME} RUNTIME DESTINATION ${BINDIR} LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})

INSTALL (FILES init/tracker.xml DESTINATION "../shared/init")
