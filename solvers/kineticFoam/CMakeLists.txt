add_executable(kineticFoam kineticFoam.C)

# twoPhaseEuler libraries
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/interfacialModels/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/phaseIncompressibleTurbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/PhaseIncompressibleTurbulenceModel")
# OF libraries
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/lnInclude")

include_directories ("${PROJECT_SOURCE_DIR}/src/kineticFluidModel/")
include_directories ("${PROJECT_SOURCE_DIR}/src/kineticFluidModel/relaxationTime/")




target_link_openfoam_libraries(kineticFoam
    fluidThermophysicalModels
    specie
    turbulenceModels
    incompressibleTurbulenceModels
    incompressibleTransportModels
    compressibleEulerianInterfacialModels
    meshTools
    finiteVolume
    compressibleTwoPhaseSystem
    phaseIncompressibleTurbulenceModels
    OpenFOAM)

  target_link_libraries(kineticFoam kineticFluidModel)

set_target_properties(kineticFoam PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$ENV{FOAM_USER_APPBIN}")
