find_package(ICU REQUIRED COMPONENTS uc i18n)

file(GLOB_RECURSE shared_SRC *.cpp)

include_directories(
	${LOOM_INCLUDE_DIR}
	SYSTEM ${GUROBI_INCLUDE_DIR}
	SYSTEM ${GLPK_INCLUDE_DIR}
	SYSTEM ${COIN_INCLUDE_DIR}
)

add_subdirectory(tests)

add_library(shared_dep ${shared_SRC})
target_link_libraries(shared_dep PUBLIC ICU::uc ICU::i18n)
