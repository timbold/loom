set(FONT_FILE "${CMAKE_SOURCE_DIR}/data/fonts/TT Norms Pro Regular.otf")
set(GENERATED_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/tt_norms_pro_regular.h")

add_custom_command(
  OUTPUT ${GENERATED_HEADER}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND xxd -i "${FONT_FILE}" > ${GENERATED_HEADER}
  COMMAND sed -i -e 's/^unsigned char/static const unsigned char/' -e 's/^unsigned int/static const unsigned int/' ${GENERATED_HEADER}
  DEPENDS ${FONT_FILE}
  VERBATIM
)

add_custom_target(label_fonts ALL DEPENDS ${GENERATED_HEADER})
