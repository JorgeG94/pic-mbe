add_subdirectory(subdir)
target_sources(
  ${main_lib} PRIVATE sample_module.f90 pic_fragment.f90 pic_mbe.f90
                      pic_mpi.f90 pic_mpi_algorithms.f90)
