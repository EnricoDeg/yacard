# Find MPI library
find_package(MPI REQUIRED COMPONENTS C)
find_library(CONFIG_LIB config)
find_path(CONFIG_INCLUDE "libconfig.h")

include_directories(${CONFIG_INCLUDE})

set(SOURCE_EXE  check.c
                config_parser.c
                utils.c
                domdcomp.c
                numerics.c
                main.c)

add_executable(canard ${SOURCE_EXE})
target_link_libraries(canard m)
target_link_libraries(canard "${CONFIG_LIB}")

# Install executables.
install (TARGETS
  canard # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

