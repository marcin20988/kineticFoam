file(GLOB lib_source_files
  kineticFluidModel.C
  relaxationTime/fullModel.C
  relaxationTime/relaxationTime.C
)

add_library(kineticFluidModel SHARED ${lib_source_files})

#component of kineticModel
include_directories ("./relaxationTime/")

#components of twoPhaseSystem
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/phaseIncompressibleTurbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/phaseIncompressibleTurbulenceModels")

#openfoam includes
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/PhaseIncompressibleTurbulenceModel")


set_target_properties (kineticFluidModel PROPERTIES LIBRARY_OUTPUT_DIRECTORY "$ENV{FOAM_USER_LIBBIN}")
