
cmake_minimum_required(VERSION 3.10)
project(OpenMMDihedralPlugin)

set(OPENMM_DIR "$ENV{CONDA_PREFIX}" CACHE PATH "Path to OpenMM installation")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

link_directories("${OPENMM_DIR}/lib" )

include_directories("${OPENMM_DIR}/include")
include_directories("${OPENMM_DIR}/include/openmm/reference")  
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/kernels/Reference)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_FILES
    src/DihedralForce.cpp
    src/DihedralForceImpl.cpp
    src/kernels/Reference/ReferenceDihedralKernelFactory.cpp
    src/kernels/Reference/ReferenceDihedralKernels.cpp
)

# Build the shared plugin
add_library(DihedralPlugin SHARED ${SOURCE_FILES})

# Link against the OpenMM base library
target_link_libraries(DihedralPlugin PRIVATE ${OPENMM_DIR}/lib/libOpenMM.so)

# Install plugin library to custom folder
install(TARGETS DihedralPlugin DESTINATION lib/plugins)


# ---- Build Python wheel and install with pip ----
set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")

add_custom_target(BuildPythonWheel ALL
    COMMAND "${Python_EXECUTABLE}" setup.py bdist_wheel
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python"
    DEPENDS DihedralPlugin
    COMMENT "Building Python wheel for dihedralplugin"
)


# --- Check for pip availability and version ---
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pip --version
    OUTPUT_VARIABLE PIP_VERSION_OUTPUT
    ERROR_VARIABLE PIP_VERSION_ERROR
    RESULT_VARIABLE PIP_CHECK_RESULT
)

set(PIP_OK TRUE)
if (NOT PIP_CHECK_RESULT EQUAL 0)
    message(WARNING "pip not found for ${Python_EXECUTABLE}. Skipping wheel installation.")
    set(PIP_OK FALSE)
else()
    string(REGEX MATCH "[0-9]+\\.[0-9]+" PIP_VERSION "${PIP_VERSION_OUTPUT}")
    message(STATUS "Found pip version ${PIP_VERSION}")

    if (PIP_VERSION VERSION_LESS "20.0")
        message(WARNING  "pip version ${PIP_VERSION} is too old (<20.0). Skipping installation.")
        set(PIP_OK FALSE)
    endif()
endif()

# --- Conditionally add pip install target ---
if (PIP_OK)
   add_custom_target(InstallPythonWheel ALL
    COMMAND "${Python_EXECUTABLE}" -m pip install --force-reinstall dist/*.whl
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python"
    DEPENDS BuildPythonWheel
    COMMENT "Installing Python wheel"
    )
else()
    message(STATUS " Skipping pip install step due to missing or outdated pip.")
endif()

