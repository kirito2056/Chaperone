qt_add_executable(chaperone-gui main.cpp)

target_compile_definitions(chaperone-gui PRIVATE
    CHAPERONE_QML_DIR="${CMAKE_CURRENT_SOURCE_DIR}/qml"
    CHAPERONE_PDB_PATH="${CMAKE_SOURCE_DIR}/data/pdb/1UBQ.pdb"
)

target_link_libraries(chaperone-gui PRIVATE
    Qt6::Quick
    Qt6::Quick3D
    chaperone_qml_module
)

qt_import_qml_plugins(chaperone-gui)
