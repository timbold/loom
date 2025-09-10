include_directories(
	${LOOM_INCLUDE_DIR}
)

add_executable(transitmapTest TestMain.cpp SanitizeSvgTest.cpp DirMarkerTest.cpp DropOverlappingStationsTest.cpp ArrowHeadDirectionTest.cpp BgMapTest.cpp LandmarkProjectionTest.cpp LandmarkSizeTest.cpp)
target_link_libraries(transitmapTest transitmap_dep)
