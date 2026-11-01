cmake_minimum_required(VERSION 3.20)

project(SCATMECH_GUI_Native VERSION 7.23 LANGUAGES CXX)

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
foreach(config Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER "${config}" config_upper)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin")
endforeach()

add_library(scatmech STATIC
    scatmech/code/allrough.cpp
    scatmech/code/askuser.cpp
    scatmech/code/axifree.cpp
    scatmech/code/axipart1.cpp
    scatmech/code/axipart2.cpp
    scatmech/code/axisym.cpp
    scatmech/code/bobvlieg1.cpp
    scatmech/code/bobvlieg2.cpp
    scatmech/code/bobvlieg3.cpp
    scatmech/code/brdf.cpp
    scatmech/code/coatedmie.cpp
    scatmech/code/crossgrating.cpp
    scatmech/code/crossgrating2.cpp
    scatmech/code/crossrcw.cpp
    scatmech/code/crough.cpp
    scatmech/code/dielfunc.cpp
    scatmech/code/diffuse.cpp
    scatmech/code/facet.cpp
    scatmech/code/fft.cpp
    scatmech/code/filmtran.cpp
    scatmech/code/finiteaperture.cpp
    scatmech/code/firstdiffuse.cpp
    scatmech/code/flake.cpp
    scatmech/code/focussedbeam.cpp
    scatmech/code/fresnel.cpp
    scatmech/code/gcross.cpp
    scatmech/code/grating.cpp
    scatmech/code/inherit.cpp
    scatmech/code/instrument.cpp
    scatmech/code/jmatrix.cpp
    scatmech/code/jvector.cpp
    scatmech/code/lambert.cpp
    scatmech/code/local.cpp
    scatmech/code/matrixmath.cpp
    scatmech/code/matrixmath2.cpp
    scatmech/code/miescat.cpp
    scatmech/code/models.cpp
    scatmech/code/mueller.cpp
    scatmech/code/nsphere.cpp
    scatmech/code/oasphere.cpp
    scatmech/code/onelayer.cpp
    scatmech/code/phasefunction.cpp
    scatmech/code/polydisperse.cpp
    scatmech/code/random.cpp
    scatmech/code/raygscat.cpp
    scatmech/code/rayinst.cpp
    scatmech/code/rayscat.cpp
    scatmech/code/raystack.cpp
    scatmech/code/rcw.cpp
    scatmech/code/reflectance.cpp
    scatmech/code/reg_brdf.cpp
    scatmech/code/reg_facet.cpp
    scatmech/code/reg_instrument.cpp
    scatmech/code/reg_lambert.cpp
    scatmech/code/reg_local.cpp
    scatmech/code/reg_rough.cpp
    scatmech/code/reg_sphrscat.cpp
    scatmech/code/rough.cpp
    scatmech/code/roughnes.cpp
    scatmech/code/scateval.cpp
    scatmech/code/scatmatrix.cpp
    scatmech/code/scattabl.cpp
    scatmech/code/sizedistribution.cpp
    scatmech/code/sphdfct.cpp
    scatmech/code/sphprt.cpp
    scatmech/code/sphrscat.cpp
    scatmech/code/stokes.cpp
    scatmech/code/subbobvlieg.cpp
    scatmech/code/subsphere.cpp
    scatmech/code/tmatrix.cpp
    scatmech/code/torrspar.cpp
    scatmech/code/transmit.cpp
    scatmech/code/two_source.cpp
    scatmech/code/twoface.cpp
    scatmech/code/urough.cpp
    scatmech/code/vector3d.cpp
    scatmech/code/zernike.cpp
    scatmech/code/zernikeexpansion.cpp
)
target_include_directories(scatmech PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/scatmech/code")
if(MSVC)
    target_compile_definitions(scatmech PUBLIC _CRT_SECURE_NO_WARNINGS NOMINMAX)
    target_compile_options(scatmech PUBLIC /utf-8 /bigobj)
endif()

add_executable(brdfprog scatmech/code/BRDFProg/BRDFProg.cpp)
add_executable(mieprog scatmech/code/MieProg/MieProg.cpp)
add_executable(reflectprog scatmech/code/ReflectProg/ReflectProg.cpp)
add_executable(rcwprog scatmech/code/RCWProg/RCWProg.cpp)
add_executable(geometry_export geometry_export.cpp)
foreach(program brdfprog mieprog reflectprog rcwprog geometry_export)
    target_link_libraries(${program} PRIVATE scatmech)
    if(MSVC)
        target_link_options(${program} PRIVATE
            /MANIFEST:EMBED
            "/MANIFESTINPUT:${CMAKE_CURRENT_SOURCE_DIR}/utf8.manifest")
        set_property(TARGET ${program} APPEND PROPERTY LINK_DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/utf8.manifest")
    endif()
endforeach()

install(TARGETS brdfprog mieprog reflectprog rcwprog geometry_export
        RUNTIME DESTINATION bin)
install(FILES SCATMECH_LICENSE.txt THIRD_PARTY_NOTICES.txt SOURCE_NOTICES.txt
        DESTINATION licenses)
