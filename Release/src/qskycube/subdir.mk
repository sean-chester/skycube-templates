################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/qskycube/BSkyTreeEM.cc \
../src/qskycube/BSkyTreeP.cc \
../src/qskycube/Common.cc \
../src/qskycube/PivotSelection.cc \
../src/qskycube/SkyTreeUtil.cc 

CC_DEPS += \
./src/qskycube/BSkyTreeEM.d \
./src/qskycube/BSkyTreeP.d \
./src/qskycube/Common.d \
./src/qskycube/PivotSelection.d \
./src/qskycube/SkyTreeUtil.d 

OBJS += \
./src/qskycube/BSkyTreeEM.o \
./src/qskycube/BSkyTreeP.o \
./src/qskycube/Common.o \
./src/qskycube/PivotSelection.o \
./src/qskycube/SkyTreeUtil.o 


# Each subdirectory must supply rules for building sources it contributes
src/qskycube/%.o: ../src/qskycube/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -D__cplusplus=201103L -O3 -ccbin /usr/bin/gcc-5 -Xcompiler -m64 -Xcompiler -nostdlib -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mavx2 -std=c++11 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52  -odir "src/qskycube" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -D__cplusplus=201103L -O3 -ccbin /usr/bin/gcc-5 -Xcompiler -m64 -Xcompiler -nostdlib -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mavx2 -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


