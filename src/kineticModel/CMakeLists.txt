file(GLOB lib_source_files
  kineticModel.C
  relaxationTime/*.C
)

add_library(kineticModel SHARED ${lib_source_files})

#components of twoPhaseSystem
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")

#openfoam includes
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/lnInclude")


set_target_properties (kineticModel PROPERTIES LIBRARY_OUTPUT_DIRECTORY "$ENV{FOAM_USER_LIBBIN}")
