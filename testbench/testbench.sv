/***********************************************
*
* NAME        : testbench.sv
* AUTHOR      : Prasanth Prabu Ravichandiran
* EMAIL       : pravich2@ncsu.edu
* AFFILIATION : North Carolina State University, Raleigh, NC.
* 
* DESCRIPTION : This module is used to test the mini project for ECE564  
*               for a single matrix multiplication.  
*             
*
* NOTES:
*
* REVISION HISTORY
*    Date         Programmer        Description
*    01/22/2024   Prasanth Prabu    TB for matrix multiplication
*
************************************************/
`include "/mnt/ncsudrive/g/gsetty/fall_2024/564/hw6/mini_project/rtl/common/common.vh"
`include "/mnt/ncsudrive/g/gsetty/fall_2024/564/hw6/mini_project/testbench/sram.sv"

module tb_top();


  parameter CLK_PHASE=5;
  parameter ADDR_464=12'h000;
  parameter MAX_ROUNDS=200;
  
  // Evaluation variables
  time computeCycle;
  event computeStart;
  event computeEnd;
  event checkFinish;
  time startTime;
  time endTime;

  // Testbench control variables
  event simulationStart;
  event testStart;
  integer totalNumOfCases=0;
  integer totalNumOfPasses=0;
  real epsilon_mult=1.0;          // Overridden by Makefile
  shortreal result;
  integer info_level=0;

  
  // Testbench configuration variables 
  string input_dir;               // Overridden by Makefile
  string output_dir;              // Overridden by Makefile
  integer rounds=1;
  integer timeout=100000000;      // Overridden by Makefile 
  integer num_of_testcases = 1;   // Overridden by Makefile
  integer test_mode = 1;
  integer test_number = 2;
  integer mystery_test = 0;

  bit  [31:0 ]     mem     [int] ;

  integer num_of_result_elements = 0;
  integer num_of_matching_result_elements = 0;
  //---------------------------------------------------------------------------
  // General
  //
  reg                                   clk            ;
  reg                                   reset_n        ;
  reg                                   dut_valid        ;
  wire                                  dut_ready       ;
  
  //--------------------------------------------------------------------------
  //---------------------- sram_input ---------------------------------------------
  wire                                  dut__tb__sram_input_write_enable  ;
  wire [`SRAM_ADDR_RANGE    ]           dut__tb__sram_input_write_address ;
  wire [`SRAM_DATA_RANGE    ]           dut__tb__sram_input_write_data    ;
  wire [`SRAM_ADDR_RANGE    ]           dut__tb__sram_input_read_address  ; 
  wire [`SRAM_DATA_RANGE    ]           tb__dut__sram_input_read_data     ;
  
  //---------------------------------------------------------------------------

  //---------------------- sram_weight ---------------------------------------------
  wire                                  dut__tb__sram_weight_write_enable  ;
  wire [`SRAM_ADDR_RANGE    ]           dut__tb__sram_weight_write_address ;
  wire [`SRAM_DATA_RANGE    ]           dut__tb__sram_weight_write_data    ;
  wire [`SRAM_ADDR_RANGE    ]           dut__tb__sram_weight_read_address  ; 
  wire [`SRAM_DATA_RANGE    ]           tb__dut__sram_weight_read_data     ;
  
  //---------------------------------------------------------------------------

  //---------------------- sram_result ---------------------------------------------
  wire                                  dut__tb__sram_result_write_enable  ;
  wire [`SRAM_ADDR_RANGE    ]           dut__tb__sram_result_write_address ;
  wire [`SRAM_DATA_RANGE    ]           dut__tb__sram_result_write_data    ;
  wire [`SRAM_ADDR_RANGE    ]           dut__tb__sram_result_read_address  ; 
  wire [`SRAM_DATA_RANGE    ]           tb__dut__sram_result_read_data     ;
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //SRAM
  //sram for inputs
  sram  #(.ADDR_WIDTH   (`SRAM_ADDR_WIDTH    ),
          .DATA_WIDTH   (`SRAM_DATA_WIDTH    ))
          sram_input_mem  (
          .write_enable ( dut__tb__sram_input_write_enable         ),
          .write_address( dut__tb__sram_input_write_address        ),
          .write_data   ( dut__tb__sram_input_write_data           ), 
          .read_address ( dut__tb__sram_input_read_address         ),
          .read_data    ( tb__dut__sram_input_read_data            ),
          .clk          ( clk                                     )
         );

  //sram for weights
  sram  #(.ADDR_WIDTH   (`SRAM_ADDR_WIDTH    ),
          .DATA_WIDTH   (`SRAM_DATA_WIDTH    ))
          sram_weight_mem  (
          .write_enable ( dut__tb__sram_weight_write_enable         ),
          .write_address( dut__tb__sram_weight_write_address        ),
          .write_data   ( dut__tb__sram_weight_write_data           ), 
          .read_address ( dut__tb__sram_weight_read_address         ),
          .read_data    ( tb__dut__sram_weight_read_data            ),
          .clk          ( clk                                     )
         );

  //sram for inputs
  sram  #(.ADDR_WIDTH   (`SRAM_ADDR_WIDTH    ),
          .DATA_WIDTH   (`SRAM_DATA_WIDTH    ))
          sram_result_mem  (
          .write_enable ( dut__tb__sram_result_write_enable         ),
          .write_address( dut__tb__sram_result_write_address        ),
          .write_data   ( dut__tb__sram_result_write_data           ), 
          .read_address ( dut__tb__sram_result_read_address         ),
          .read_data    ( tb__dut__sram_result_read_data            ),
          .clk          ( clk                                     )
         );
		 
//---------------------------------------------------------------------------
// DUT 
//---------------------------------------------------------------------------
  MyDesign dut(
//---------------------------------------------------------------------------
//System signals
          .reset_n                    (reset_n                      ),  
          .clk                        (clk                          ),

//---------------------------------------------------------------------------
//Control signals
          .dut_valid                  (dut_valid                    ), 
          .dut_ready                  (dut_ready                    ),

//---------------------------------------------------------------------------
// SRAM input interface
          .dut__tb__sram_input_write_enable       (dut__tb__sram_input_write_enable     ),
          .dut__tb__sram_input_write_address      (dut__tb__sram_input_write_address    ),
          .dut__tb__sram_input_write_data         (dut__tb__sram_input_write_data       ),
          .dut__tb__sram_input_read_address       (dut__tb__sram_input_read_address     ),
          .tb__dut__sram_input_read_data          (tb__dut__sram_input_read_data        ),

//---------------------------------------------------------------------------
// SRAM weight interface
          .dut__tb__sram_weight_write_enable       (dut__tb__sram_weight_write_enable     ),
          .dut__tb__sram_weight_write_address      (dut__tb__sram_weight_write_address    ),
          .dut__tb__sram_weight_write_data         (dut__tb__sram_weight_write_data       ),
          .dut__tb__sram_weight_read_address       (dut__tb__sram_weight_read_address     ),
          .tb__dut__sram_weight_read_data          (tb__dut__sram_weight_read_data        ),

//---------------------------------------------------------------------------
// SRAM result interface
          .dut__tb__sram_result_write_enable       (dut__tb__sram_result_write_enable     ),
          .dut__tb__sram_result_write_address      (dut__tb__sram_result_write_address    ),
          .dut__tb__sram_result_write_data         (dut__tb__sram_result_write_data       ),
          .dut__tb__sram_result_read_address       (dut__tb__sram_result_read_address     ),
          .tb__dut__sram_result_read_data          (tb__dut__sram_result_read_data        )
         );

       
  //---------------------------------------------------------------------------
  //  clk
  initial 
    begin
        clk                     = 1'b0;
        forever # CLK_PHASE clk = ~clk;
    end

	initial begin
		$dumpfile("wave.fsdb");
		$dumpvars();
	end 

  //---------------------------------------------------------------------------
  // get runtime args 
  initial
  begin
    #1;
    if($value$plusargs("TIMEOUT=%d",timeout));
    if($value$plusargs("input_dir=%s",input_dir));
    if($value$plusargs("num_of_testcases=%d",num_of_testcases));
    if($value$plusargs("info_level=%d",info_level));
    if($value$plusargs("test_mode=%d",test_mode));
    if($value$plusargs("test_number=%d",test_number));
    $display("INFO: number of testcases: %d",num_of_testcases);
    if($value$plusargs("epsilon_mult=%f",epsilon_mult));
    if($value$plusargs("mystery_test=%d",mystery_test))

    repeat (5) @(posedge clk);
    ->simulationStart;
    @testStart
    wait_n_clks(timeout);
    $display("###################################");
    $display("             TIMEOUT               ");
    $display("###################################");
    $finish();
  end
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Stimulus

  task wait_n_clks;
    input integer i;
  begin
    repeat(i)
    begin
      wait(clk);
      wait(!clk);
    end
  end
  endtask

  task handshake;
  begin
    wait(!clk);
    dut_valid = 1;
    wait(clk);
    wait(!dut_ready);
    wait(!clk);
    dut_valid = 0;
    wait(clk);
    wait(dut_ready);
    wait(!clk);
    wait(clk);
  end
  endtask

  function void check_result(integer testNum, string dirname, string filename);
    shortreal epsilon;
    shortreal expected_result;
    shortreal dut_result;
    shortreal difference;
    
    mem.delete();
    $display("Received filename in check_result = %s", filename);

    epsilon = $bitstoshortreal(64'h3CB0_0000);
    $readmemh($sformatf("%s/%s_%0d/%s_%0d_result.dat",input_dir, dirname, testNum, filename, testNum), mem);
    $display($sformatf("INFO: %s/%s_%0d/%s_%0d_result.dat",input_dir, dirname, testNum, filename, testNum));

    foreach (mem[key]) begin
      num_of_result_elements ++;
      if (sram_result_mem.mem.exists(key)) begin
        expected_result = $bitstoshortreal(mem[key]&32'h7fff_ffff); //Masking sign bit
        dut_result      = $bitstoshortreal(sram_result_mem.mem[key]&32'h7fff_ffff); //Masking sign bit
        difference      = expected_result - dut_result;
        if (((-epsilon * epsilon_mult) <= difference) && (difference <= (epsilon * epsilon_mult)))  begin          
          $display("[%d] Result MATCH: expected_result = %7.20f, dut_result = %7.20f", key,expected_result, dut_result);
          num_of_matching_result_elements++;
        end
        else begin
          $display("[%d]:Difference = %7.20f, expected_result = %7.20f, dut_result = %7.20f", key, difference, expected_result, dut_result);

        end
      end
      else begin
        $display("[%d]: ERROR: SRAM result entry is missing",key);  
      end
    end

  endfunction : check_result

  task test;
    input integer testNum;
    input string dirname;
    input string filename; 
  begin
    
    $display("INFO:LVL0: ######## Running Test: %0d ########",testNum);
    wait_n_clks(10);
    sram_input_mem.loadMem($sformatf("%s/%s_%0d/%s_%0d_input.dat",input_dir, dirname,testNum, filename, testNum));
    sram_weight_mem.loadMem($sformatf("%s/%s_%0d/%s_%0d_weight.dat",input_dir, dirname, testNum, filename, testNum));
    sram_result_mem.mem.delete();
    wait_n_clks(10);
    handshake();
    wait_n_clks(10);
    check_result(.testNum(testNum), .dirname(dirname), .filename(filename));
    wait_n_clks(10);
  end
  endtask


  initial
  begin
    @(simulationStart);
    reset_n = 1;
    wait_n_clks(10);
    reset_n = 0;
    wait_n_clks(20);
    dut_valid = 0;
    wait_n_clks(20);
    reset_n = 1;
    wait_n_clks(20);
    $display("INFO: DONE WITH RESETING DUT");
    ->testStart;
    startTime=$time();

    if (!test_mode) begin
      $display("Running individual test: %d", test_number);
      test(.testNum(test_number), .dirname("test"), .filename("test")) ;
    end
    else begin
      // Running all the test cases
      for(int i=1;i<num_of_testcases+1;i++)
      begin
        test(.testNum(i), .dirname("test"), .filename("test"));
      end
      // Running mystery test cases
      if (mystery_test) begin
        for(int i=1;i<num_of_testcases+1;i++)
          begin
            test(.testNum(i), .dirname("mystery_test"), .filename("mystery_test"));
          end
      end
    end

    endTime=$time();
    if(num_of_testcases != 0)
    begin
      $display("INFO: Total number of cases       : %0d",num_of_testcases);
      $display("INFO: Total number of result pass : %0d / %0d",num_of_matching_result_elements, num_of_result_elements);
      $display("INFO: Final pass percentage       : %6.2f",(num_of_matching_result_elements * 100)/num_of_result_elements);
      $display("INFO: Final Time Result           : %0t ns",endTime-startTime);
      $display("INFO: Final Cycle Result          : %0d cycles\n",((endTime-startTime)/CLK_PHASE));
    end
    $finish();
  end
endmodule
