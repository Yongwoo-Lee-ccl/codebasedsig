set(Header api.h    rng.h   parm.h  matrix.h    rm.h    common.h    decoding.h)
set(C_FILE          rng.c           matrix.c    rm.c    common.c    decoding.c)

set(CMAKE_BUILD_TYPE Debug)

add_library(FILE ${Header} ${C_FILE})