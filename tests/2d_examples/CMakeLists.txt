# Enumerate subdirectories with CMakeLists.txt and add them.
file(GLOB CHILDREN RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*")
foreach(child ${CHILDREN})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child})
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${child}/CMakeLists.txt)
            add_subdirectory(${child})
        endif()
    endif()
endforeach()
