cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(tree_code)
enable_language(Fortran)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(trees.exe trees.F90 num_pars.F90 defined_types.F90 kind_numbers.F90 parameter_modules.F90   memory_modules.F90 tree_routines.F90 modified_merger_tree.F90 deltcrit.F90 memory.F90 sigmacdm_spline.F90 interp.F90 locate.F90 hyperbolic.F90 split_PCH.F90 ran3.F90 spline.F90 make_tree.F90 indexxx.F90 transfer_function.F90 unresolved_mass.F90 parameters.F90 growth_rate.F90)
