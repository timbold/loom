file(GLOB_RECURSE transitmap_SRC *.cpp)

set(transitmap_main TransitMapMain.cpp)

list(REMOVE_ITEM transitmap_SRC ${transitmap_main})
list(REMOVE_ITEM transitmap_SRC TestMain.cpp)

include_directories(
        ${LOOM_INCLUDE_DIR}
)

find_package(Freetype)
find_package(ICU REQUIRED COMPONENTS uc i18n)
if (FREETYPE_FOUND)
  include_directories(${FREETYPE_INCLUDE_DIRS})
endif()

add_subdirectory(label)

add_subdirectory(output/protobuf)
add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(transitmap ${transitmap_main})
add_library(transitmap_dep ${transitmap_SRC})
add_dependencies(transitmap_dep label_fonts)
target_link_libraries(transitmap_dep PUBLIC shared_dep util dot_dep)
target_link_libraries(transitmap_dep PUBLIC ICU::uc ICU::i18n)

if (FREETYPE_FOUND)
  target_link_libraries(transitmap_dep PUBLIC ${FREETYPE_LIBRARIES})
  target_compile_definitions(transitmap_dep PUBLIC LOOM_HAVE_FREETYPE)
endif()

if (Protobuf_FOUND)
  add_dependencies(transitmap_dep proto)
  target_link_libraries(transitmap_dep PUBLIC proto)
endif()

target_link_libraries(transitmap transitmap_dep)
