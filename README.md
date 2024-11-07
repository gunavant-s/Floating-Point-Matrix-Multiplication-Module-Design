
# ECE464/564 Mini Project 
This document contains the instructions and commands to setup mini_project directory. In the folder tree of this mini_project, several ```Makefile```s are used to 

## Overview
- [Unzip](#unzip)
- [Start Designing](#start-designing)
- [Synthesis](#synthesis)
- [Submission](#submission)
- [Appendix](#appendix)

## Unzip
Once you have placed ```mini_project.zip``` at desired directory. Launch a terminal at that directory and use the following command to unzip.
```bash
unzip mini_project.zip
```
You should find the unzipped mini_project folder ```mini_project/```

## Start Designing
### Setup script

```mini_project/setup.sh``` is provided to load Modelsim and Synopsys

To source the script:
```bash
source setup.sh
```
This script also enables you to <kbd>Tab</kbd> complete ```make``` commands

### Mini project description

The project description document is located in ```mini_project/docs/```

### Where to put your design

A Verilog file ```mini_project/rtl/dut.sv``` is provided with all the ports already connected to the test fixture

### How to compile your design

To compile your design

Change directory to ```mini_project/run/``` 

```bash
make build-dw
make build
```

All the .sv files in ```mini_project/rtl/``` will be compiled with this command.

### How to run your design

Run with Modelsim UI 564:
```bash
make debug
```

### Evaluation Testing
To evaluate you design headless/no-gui, change directory to ```mini_project/run/```
```
make eval
```
This will produce a set of log files that will highlight the results of your design. This should only be ran as a final step before Synthesis

All log files is in the following directory ```mini_project/run/logs```

All test resutls is in the results log file ```mini_project/run/logs/RESULTS.log```

All simulation resutls is in the following log file ```mini_project/run/logs/output.log```

All simulation info is in the following log file ```mini_project/run/logs/INFO.log```

## Synthesis

Once you have a functional design, you can synthesize it in ```mini_project/synthesis/```

### Synthesis Command
The following command will synthesize your design with a default clock period of 10 ns
```bash
make all
```
### Clock Period

To run synthesis with a different clock period
```bash
make all CLOCK_PER=<YOUR_CLOCK_PERIOD>
```
For example, the following command will set the target clock period to 4 ns.

```bash
make all CLOCK_PER=10
```

## Appendix

### Directory Rundown

You will find the following directories in ```mini_project/```

* ```inputs/``` 
  * Contains the .dat files for the input SRAMs used in HW 
* ```HW_specification/```
  * Contains the HW specification document
* ```rtl/```
  * All .v files will be compiled when executing ```make vlog-v``` in ```mini_project/run/```
  * A template ```dut.v``` that interfaces with the test fixture is provided
* ```run/```
  * Contains the ```Makefile``` to compile and simulate the design
* ```scripts/```
  * Contains the python script that generates a random input/output
* ```synthesis/```
  * The directory you will use to synthesize your design
  * Synthesis reports will be exported to ```synthesis/reports/```
  * Synthesized netlist will be generated to ```synthesis/gl/```
* ```testbench/```
  * Contains the test fixture of the HW


