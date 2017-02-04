################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/main.cu 

OBJS += \
./src/main.o 

CU_DEPS += \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -D__cplusplus=201103L -O3 -ccbin /usr/bin/gcc-5 -Xcompiler -m64 -Xcompiler -nostdlib -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mavx2 -std=c++11 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -D__cplusplus=201103L -O3 -ccbin /usr/bin/gcc-5 -Xcompiler -m64 -Xcompiler -nostdlib -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mavx2 -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


