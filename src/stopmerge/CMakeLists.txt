file(GLOB_RECURSE stopmerge_SRC *.cpp)

set(stopmerge_main StopMergeMain.cpp)

list(REMOVE_ITEM stopmerge_SRC ${stopmerge_main})

configure_file(
  "_config.h.in"
  "_config.h"
)

add_executable(stopmerge ${stopmerge_main})
add_library(stopmerge_dep ${stopmerge_SRC})

target_link_libraries(stopmerge stopmerge_dep shared_dep dot_dep util)
