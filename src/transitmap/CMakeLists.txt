file(GLOB_RECURSE transitmap_SRC *.cpp)

set(transitmap_main TransitMapMain.cpp)

list(REMOVE_ITEM transitmap_SRC ${transitmap_main})
list(REMOVE_ITEM transitmap_SRC TestMain.cpp)

include_directories(
        ${LOOM_INCLUDE_DIR}
)

find_package(Freetype REQUIRED)
include_directories(${FREETYPE_INCLUDE_DIRS})

add_subdirectory(output/protobuf)
add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(transitmap ${transitmap_main})
add_library(transitmap_dep ${transitmap_SRC})
target_link_libraries(transitmap_dep PUBLIC ${FREETYPE_LIBRARIES})

if (Protobuf_FOUND)
	add_dependencies(transitmap_dep proto)
	target_link_libraries(transitmap transitmap_dep shared_dep dot_dep util proto ${Protobuf_LIBRARIES})
else()
	target_link_libraries(transitmap transitmap_dep shared_dep dot_dep util)
endif()
