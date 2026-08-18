cmake_minimum_required(VERSION 3.24)
project(chaperone LANGUAGES CXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Qml Quick Quick3D QmlImportScanner)

find_package(CxxQt QUIET)
if(NOT CxxQt_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        CxxQt
        GIT_REPOSITORY https://github.com/kdab/cxx-qt-cmake.git
        GIT_TAG main
    )
    FetchContent_MakeAvailable(CxxQt)
endif()

cxx_qt_import_crate(
    MANIFEST_PATH Cargo.toml
    CRATES chaperone_bridge
    LOCKED
    QT_MODULES Qt::Core Qt::Gui Qt::Qml Qt::Quick
)

cxx_qt_import_qml_module(chaperone_qml_module
    URI "Chaperone"
    SOURCE_CRATE chaperone_bridge
)

add_subdirectory(gui)
