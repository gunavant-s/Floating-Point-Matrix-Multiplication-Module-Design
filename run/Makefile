COMMON_SRC		=$(wildcard ../rtl/common)
VERILOG_SRC		=$(wildcard ../rtl/*.*v)  
TESTBENCH_SRC	=$(wildcard ../testbench/*.*v)
DPI_SRC=$(wildcard ../dpi/*.c) 
DW_SRC_PATH=/mnt/apps/public/COE/synopsys_apps/syn/T-2022.03-SP4/dw/sim_ver
DW_SRC=$(shell find ${DW_SRC_PATH} -name "DW_fp*.v" -printf "%p ")
EPSILON_MULT_TEST=1.0
NUMBER_OF_TEST=4
INFO_LEVEL?=0
TIMEOUT=10000000
TEST_NUMBER=1
MYSTERY_TEST=0

syntax:
	/bin/bash -c "xmverilog $(VERILOG_SRC) -c"

load:
	/bin/bash -c ml load xcelium/23.03
	/bin/bash -c ml load synopsys/2022

all:
	/bin/bash -c "xmverilog ${TESTBENCH_SRC} +TIMEOUT=${TIMEOUT} +num_of_testcases=${NUMBER_OF_TEST} +test_mode=1 +epsilon_mult=${EPSILON_MULT_TEST}  +input_dir=../inputs +info_level=${INFO_LEVEL} +mystery_test=$(MYSTERY_TEST) $(VERILOG_SRC) +access+r"

TB:
	/bin/bash -c "xmverilog ${TESTBENCH_SRC} +TIMEOUT=${TIMEOUT} +num_of_testcases=2 +epsilon_mult=${EPSILON_MULT_TEST} +test_mode=0 +test_number=${TEST_NUMBER} +input_dir=../inputs +info_level=${INFO_LEVEL} +mystery_test=$(MYSTERY_TEST) $(VERILOG_SRC) +access+r"

build-dw:
	/bin/bash -c "vlog ${DW_SRC} "	

build:
	/bin/bash -c "vlog -sv +incdir+$(COMMON_SRC) $(TESTBENCH_SRC) $(VERILOG_SRC)"	

eval: build
	mkdir -p logs logs/test 
	/bin/bash -c "vsim -c -voptargs=+acc tb_top +TIMEOUT=${TIMEOUT} +num_of_testcases=${NUMBER_OF_TEST} +test_mode=1 +epsilon_mult=${EPSILON_MULT_TEST}  +input_dir=../inputs +info_level=${INFO_LEVEL} +mystery_test=$(MYSTERY_TEST) -do \"run -all\" | tee logs/output.log "
	grep -r "INFO" logs/output.log > logs/INFO.log
	grep -r "INFO" logs/output.log | grep  Result > logs/RESULTS.log

eval_individual: build
	mkdir -p logs logs/test 
	/bin/bash -c "vsim -c -voptargs=+acc tb_top +TIMEOUT=${TIMEOUT} +num_of_testcases=1 +epsilon_mult=${EPSILON_MULT_TEST} +test_mode=0 +test_number=${TEST_NUMBER} +input_dir=../inputs +info_level=${INFO_LEVEL} +mystery_test=$(MYSTERY_TEST) -do \"run -all\" | tee logs/output.log "
	grep -r "INFO" logs/output.log > logs/INFO.log
	grep -r "INFO" logs/output.log | grep  Result > logs/RESULTS.log

debug: build
	/bin/bash -c "vsim -voptargs=+acc tb_top +TIMEOUT=${TIMEOUT} +num_of_testcases=${NUMBER_OF_TEST} +epsilon_mult=${EPSILON_MULT_TEST}  +input_dir=../inputs +mystery_test=$(MYSTERY_TEST) -do \"log -r *; do wave.do\"" 


clean:
	rm -rf work logs results transcript
