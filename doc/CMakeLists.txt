cmake_minimum_required(VERSION 2.7)

include(UseLATEX.cmake)

SET(MyDocImages
  hhblits-hhsearch.png
)

ADD_LATEX_DOCUMENT(hhsuite-userguide.tex
	DEFAULT_PDF
	BIBFILES hhsuite-userguide.bib
	IMAGES ${MyDocImages}
)

#add_custom_command(TARGET pdf POST_BUILD
#	COMMAND ${CMAKE_COMMAND} -E copy
#	${CMAKE_CURRENT_BINARY_DIR}/hhsuite-userguide.pdf
#	${CMAKE_BINARY_DIR}/hhsuite-userguide.pdf
#)
