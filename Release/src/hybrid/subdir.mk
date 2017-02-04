################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hybrid/mommyq.cpp 

OBJS += \
./src/hybrid/mommyq.o 

CPP_DEPS += \
./src/hybrid/mommyq.d 


# Each subdirectory must supply rules for building sources it contributes
src/hybrid/%.o: ../src/hybrid/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -D__cplusplus=201103L -O3 -ccbin /usr/bin/gcc-5 -Xcompiler -m64 -Xcompiler -nostdlib -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mavx2 -std=c++11 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52  -odir "src/hybrid" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -D__cplusplus=201103L -O3 -ccbin /usr/bin/gcc-5 -Xcompiler -m64 -Xcompiler -nostdlib -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mavx2 -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


