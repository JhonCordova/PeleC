set(PELE_PHYSICS_EOS_MODEL GammaLaw)
set(PELE_PHYSICS_CHEMISTRY_MODEL Null)
set(PELE_PHYSICS_TRANSPORT_MODEL Constant)
set(PELE_PHYSICS_ENABLE_SOOT OFF)
set(PELE_PHYSICS_ENABLE_SPRAY OFF)
set(PELE_PHYSICS_SPRAY_FUEL_NUM 0)
include(BuildExeAndLib)

target_sources(${pele_exe_name}
  PRIVATE
  unit-tests-main.cpp
  test-config.cpp
  test-odtles-local-engine.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTDiffusion.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTEventSampler.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTKernel.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTLineGeometry.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTLineState.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTManager.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTMomentExtractor.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTParams.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTReconcile.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTStepper.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTThermoBridge.cpp
  ${CMAKE_SOURCE_DIR}/Source/ODTLES/ODTTripletMap.cpp
  )

if(PELE_ENABLE_CUDA)
  set_source_files_properties(
    unit-tests-main.cpp
    test-config.cpp
    test-odtles-local-engine.cpp
    PROPERTIES LANGUAGE CUDA)
endif()

find_package(GTest REQUIRED)
target_link_libraries(${pele_exe_name} PRIVATE GTest::gtest)
