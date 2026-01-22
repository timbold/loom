file(GLOB_RECURSE transitmap_SRC *.cpp)

set(transitmap_main TransitMapMain.cpp)

list(REMOVE_ITEM transitmap_SRC ${transitmap_main})
list(REMOVE_ITEM transitmap_SRC TestMain.cpp)

include_directories(
	${LOOM_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(transitmap ${transitmap_main})
add_library(transitmap_dep ${transitmap_SRC})

find_path(SQLITE3_INCLUDE_DIR sqlite3.h)
find_library(SQLITE3_LIBRARY sqlite3)
if (SQLITE3_INCLUDE_DIR AND SQLITE3_LIBRARY)
	target_include_directories(transitmap_dep SYSTEM PRIVATE ${SQLITE3_INCLUDE_DIR})
	target_link_libraries(transitmap_dep ${SQLITE3_LIBRARY})
	target_link_libraries(transitmap ${SQLITE3_LIBRARY})
else()
	message(FATAL_ERROR "sqlite3 library not found")
endif()

find_package(PNG)
if (PNG_FOUND)
	target_include_directories(transitmap_dep SYSTEM PRIVATE ${PNG_INCLUDE_DIRS})
	target_link_libraries(transitmap_dep ${PNG_LIBRARIES})
	target_link_libraries(transitmap ${PNG_LIBRARIES})
else()
	message(FATAL_ERROR "libpng not found")
endif()

find_package(JPEG)
if (JPEG_FOUND)
	target_include_directories(transitmap_dep SYSTEM PRIVATE ${JPEG_INCLUDE_DIR})
	target_link_libraries(transitmap_dep ${JPEG_LIBRARIES})
	target_link_libraries(transitmap ${JPEG_LIBRARIES})
else()
	message(FATAL_ERROR "libjpeg not found")
endif()

target_link_libraries(transitmap transitmap_dep shared_dep dot_dep util)
