set(SOURCES
    ndct.cpp
    fdct.cpp
    sdct.cpp
)

set(library dct)
add_library(${library} STATIC)
target_include_directories(${library} PRIVATE ${dct_SOURCE_DIR}/include)
target_sources(${library} PRIVATE ${SOURCES})
