cmake_minimum_required(VERSION 3.12)
project(AnimalProcessingToolbox)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenCV REQUIRED)
find_package(MPI REQUIRED)

add_executable(apt main.cpp include/TransformationMatrix.h TransformationMatrix.cpp include/ImageTransformation.h ImageTransformation.cpp include/ImageDataClass.h ImageDataClass.cpp include/BorderType.h)
target_link_libraries(apt ${OpenCV_LIBS} ${MPI_LIBRARIES})

