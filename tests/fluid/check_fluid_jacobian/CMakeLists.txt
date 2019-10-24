add_executable(test_fluid_jacobian
        test_fluid_jacobian.cpp
        ${MAST_ROOT_DIR}/examples/fluid/meshing/naca0012.cpp
        ${MAST_ROOT_DIR}/examples/fluid/meshing/naca0012_wing.cpp
        ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_2D.cpp
        ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_3D.cpp
        ${MAST_ROOT_DIR}/examples/fluid/meshing/cylinder.cpp
        ${MAST_ROOT_DIR}/examples/fluid/meshing/mesh_initializer.cpp
        )

target_include_directories(test_fluid_jacobian PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/tests/fluid/base
        ${PROJECT_SOURCE_DIR}/tests
        )
target_link_libraries(test_fluid_jacobian mast)

install(TARGETS test_fluid_jacobian
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/fluid/check_fluid_jacobian)