//---------------------------------------------------------------------------
// Filename    : dut.sv
// Author      : Gunavant Setty
// Affiliation : North Carolina State University, Raleigh, NC
// Date        : Nov 2024
// Email       : gsetty@ncsu.edu
    
// Description : This file contains the matrix multiplication of FP.
//               It writes the accumulation result back to the result SRAM (Matrix C). 
// DUT - Mini project
//---------------------------------------------------------------------------
`include "/mnt/ncsudrive/g/gsetty/fall_2024/564/hw6/mini_project/rtl/common/common.vh"
`include "/mnt/apps/public/COE/synopsys_apps/syn/T-2022.03-SP4/dw/sim_ver/DW_fp_mac.v"
`include "/mnt/apps/public/COE/synopsys_apps/syn/T-2022.03-SP4/dw/sim_ver/DW_fp_dp2.v"

module MyDesign(
//---------------------------------------------------------------------------
//System signals
  input wire reset_n                      ,
  input wire clk                          ,

//---------------------------------------------------------------------------
//Control signals
  input wire dut_valid                    ,
  output wire dut_ready                   ,

//---------------------------------------------------------------------------
//input SRAM  i.e A
  output wire                           dut__tb__sram_input_write_enable  ,
  output wire [`SRAM_ADDR_RANGE     ]   dut__tb__sram_input_write_address ,
  output wire [`SRAM_DATA_RANGE     ]   dut__tb__sram_input_write_data    ,
  output wire [`SRAM_ADDR_RANGE     ]   dut__tb__sram_input_read_address  ,
  input  wire [`SRAM_DATA_RANGE     ]   tb__dut__sram_input_read_data     ,

//weight SRAM  i.e B
  output wire                           dut__tb__sram_weight_write_enable  ,
  output wire [`SRAM_ADDR_RANGE     ]   dut__tb__sram_weight_write_address ,
  output wire [`SRAM_DATA_RANGE     ]   dut__tb__sram_weight_write_data    ,
  output wire [`SRAM_ADDR_RANGE     ]   dut__tb__sram_weight_read_address  ,
  input  wire [`SRAM_DATA_RANGE     ]   tb__dut__sram_weight_read_data     ,

//result SRAM interface i.e C
  output wire                           dut__tb__sram_result_write_enable  ,
  output wire [`SRAM_ADDR_RANGE     ]   dut__tb__sram_result_write_address ,
  output wire [`SRAM_DATA_RANGE     ]   dut__tb__sram_result_write_data    ,
  output wire [`SRAM_ADDR_RANGE     ]   dut__tb__sram_result_read_address  ,
  input  wire [`SRAM_DATA_RANGE     ]   tb__dut__sram_result_read_data

);

//q_state_output SRAM interface
	reg [15:0] sram_write_address_r ;
	reg [15:0] sram_read_address_A  ;
	reg [15:0] sram_read_address_B  ;	
		
	reg [`SRAM_DATA_RANGE 			] input_read_data; // 32-bits
	reg [`SRAM_DATA_RANGE 			] weight_read_data; // 32-bits
	reg [`SRAM_DATA_RANGE 			] accum_result; // C Result from DW_fp_mac
	reg [`SRAM_DATA_RANGE 			] mac_result_z;
	reg [`SRAM_DATA_RANGE			] sram_result_write_data_reg;
  	reg compute_complete;

// Use define to parameterize the variable sizes
  `ifndef SRAM_ADDR_WIDTH
    `define SRAM_ADDR_WIDTH 16
  `endif

  `ifndef SRAM_DATA_WIDTH
    `define SRAM_DATA_WIDTH 32
  `endif
	
	parameter [2:0]         //synopsys enum states
		S0	   = 3'b000,
		S1	   = 3'b001,
		S2     = 3'b010,
		S3 	   = 3'b011,
		S4     = 3'b100,
		S5     = 3'b101,
		S5_5   = 3'b110,
		S6     = 3'b111;

	reg [3:0] current_state, next_state;

	// Local control path variables
	reg                           set_dut_ready             ;
	reg [1:0]                     input_read_addr_sel       ;
	reg [1:0]                     weight_read_addr_sel      ;
	reg [1:0]                     temp_count_row_A_select   ; // for A count
	reg [1:0] 					  temp_count_var_A_count_select; // for increasing A until Arow
	reg [1:0] 					  temp_count_var_C_count_select; // moving to next Crow
	reg [1:0] 					  temp_count_z_select;
	reg                           write_enable_sel          ;

// Local data path variables 
	reg [`SRAM_DATA_WIDTH-1:0]  array_size; //32 bits
	
	// reg [31:0] temp_count_row_A;
	// reg [31:0] coloums_count; // for decrementing the count of number of column so that when it (ends finish) decreases then A start = 1 + row ? i.e next col of c matrix
	reg [31:0] temp_count_z; // counts while c = a*b + c, then when equal to row_B then write to Z and make it 0 again
	reg [31:0] temp_count_var_A_count; // if equals to brows then start from A1
	reg [31:0] temp_count_var_C_count; // If equals to Brows then move to next row of C, then 
	// reg [31:0] temp_accumulate;

	reg [31:0] rows_A, column_A;
	reg [31:0] rows_B, column_B;
	reg [31:0] rows_C, column_C;
	
	reg [15:0] input_read_address, weight_read_address;
	reg        write_enable_sel_r;
	reg [15:0] start_addr;
	reg [31:0] accum_result_r;
	reg [31:0] accum_result_r2;
	reg [31:0] temp;

//------------------------------------------Address Counters-----------------------------------------------------
	// This temp variable is to track A, i.e if calculation of a C is done then start from 1,
	//  and temp_count_var_C_count = Brows then A start = 1 + ROW
	always @(posedge clk) begin 
		if(!reset_n)
			temp_count_var_A_count <= 0;
		else begin
			if(temp_count_var_A_count_select == 2'b00)
				temp_count_var_A_count <= `SRAM_ADDR_WIDTH'b0;
			else if (temp_count_var_A_count_select == 2'b01)
				temp_count_var_A_count <= temp_count_var_A_count + `SRAM_ADDR_WIDTH'b1;
			else if (temp_count_var_A_count_select == 2'b10)
				temp_count_var_A_count <= temp_count_var_A_count;
			else if (temp_count_var_A_count_select == 2'b11) // when state S5
				temp_count_var_A_count <= `SRAM_ADDR_WIDTH'b01;
		end
	end

	// This temp variable is to track C rows, i.e if calculation of c row is done then move to next row
	//  Now A start = row + 1
	always @(posedge clk) begin 
		if(!reset_n)
			temp_count_var_C_count <= 0;
		else begin
			if(temp_count_var_C_count_select == 2'b00)
				temp_count_var_C_count <= `SRAM_ADDR_WIDTH'b0;
			else if (temp_count_var_C_count_select == 2'b01) begin
				temp_count_var_C_count <= temp_count_var_C_count + `SRAM_ADDR_WIDTH'b1;
			end
			else if (temp_count_var_C_count_select == 2'b10)
				temp_count_var_C_count <= temp_count_var_C_count;
			else if (temp_count_var_C_count_select == 2'b11) begin
				temp_count_var_C_count <= `SRAM_ADDR_WIDTH'b01;
			end
		end
	end
	// another counter for setting where z = a*b +c. If it reaches final value i.e Brows then Write_C then next C calculation -> S2
	// when temp_a == rowB then write to C

	always @(posedge clk) begin
		if(!reset_n) begin
			temp_count_z <= 0;
			accum_result <= 0;
		end
		else begin
			if(temp_count_z_select == 2'b00) begin
				if(temp_count_var_A_count==1 && write_enable_sel == 1) begin
					accum_result <= 0;
					accum_result_r <=0 ;
					temp_count_z <= 0;
				end
				else begin
					temp_count_z <= 0;
					// mac_result_z <= 0;
					accum_result_r <=0;
					accum_result <= 0;
				end
			end
			else if(temp_count_z_select == 2'b01) begin
				if(temp_count_var_A_count==1 && write_enable_sel == 1) begin
					accum_result <= 0;
					accum_result_r <=0 ;
					temp_count_z <= 0;
				end
				else begin
				temp_count_z <= temp_count_z + 1;
				accum_result <= mac_result_z;
				end
			end
			else if(temp_count_z_select == 2'b10) begin
				if(temp_count_var_A_count==1 && write_enable_sel == 1) begin
					accum_result <= 0;
					accum_result_r <=0 ;
					temp_count_z <= 0;
				end
				else begin
					temp_count_z <= temp_count_z;
					if(write_enable_sel_r == 1 && temp_count_var_C_count != column_B - 1 && temp_count_var_A_count==1)
						accum_result <= accum_result;//accum_result_r
					else if(write_enable_sel_r == 1 && temp_count_var_C_count == column_B -1)
						accum_result <= 0;
					else 
						accum_result <= accum_result;
					end
				end
			else
				temp_count_z <= 0;
		end
	end
	// store the accumulate in temp reg C and keep on adding 

	// SRAM A read address calculation
	always @(posedge clk) begin
		if (!reset_n) begin
			sram_read_address_A <= `SRAM_ADDR_WIDTH'b0;
		end
		else begin
			if (input_read_addr_sel == 2'b00)
				sram_read_address_A <= `SRAM_ADDR_WIDTH'b0;
			else if (input_read_addr_sel == 2'b01)
				sram_read_address_A <= sram_read_address_A + `SRAM_ADDR_WIDTH'b1;
			else if (input_read_addr_sel == 2'b10)
				sram_read_address_A <= sram_read_address_A;
			else if (input_read_addr_sel == 2'b11)
				sram_read_address_A <= sram_read_address_A-(column_B-`SRAM_ADDR_WIDTH'b01); //SRAM_ADDR_WIDTH'b10
		end
	end


	// SRAM B read address calculation
	always @(posedge clk) begin
		if (!reset_n) begin
			sram_read_address_B <= 0;
		end
		else begin
			if (weight_read_addr_sel == 2'b00)
				sram_read_address_B <= `SRAM_ADDR_WIDTH'b0;
			else if (weight_read_addr_sel == 2'b01) begin
				if(sram_read_address_B == column_B * rows_B)
					sram_read_address_B <= `SRAM_ADDR_WIDTH'b01;
				else
					sram_read_address_B <= sram_read_address_B + `SRAM_ADDR_WIDTH'b1;
			end
			else if (weight_read_addr_sel == 2'b10) 
				sram_read_address_B <= sram_read_address_B;
			else if (weight_read_addr_sel == 2'b11) //if counter_C == row_B
				sram_read_address_B <= `SRAM_ADDR_WIDTH'b01;
				// sram_read_address_B <= `SRAM_ADDR_WIDTH'b00;
		end
	end
	// for C write addr calculation
	always @(posedge clk) begin
		if (!reset_n) begin
			sram_write_address_r <= 16'hFFFF;  
		end
		
		else begin
			if(current_state == S0) begin
				sram_write_address_r <= 16'hFFFF;
			end
			if (write_enable_sel) begin
				//sram_result_write_data_reg <= mac_result_z; //accum_result
				sram_write_address_r <= sram_write_address_r + 16'b1;
			end
		end
	end


//------------------------------------------Counters ^-----------------------------------------------------

	always @(posedge clk or negedge reset_n) begin
		if(!reset_n)
			current_state <= S0;
		else
			current_state <= next_state;
	end	
//------------------------------------------------Changes state ^--------------------------------------------

	always @(*) begin
		case(current_state) //prevent unwanted latches
			S0: begin
				if(dut_valid) begin
					rows_A							= 0;
					rows_B							= 0;
					column_A				    	= 0;
					column_B 			        	= 0;
					column_C 			        	= 0;
					start_addr							= 1;
					sram_result_write_data_reg 		= 0;
         			// sram_write_address_r            = -1;
					input_read_addr_sel  			= 2'b0;
					set_dut_ready    			    = 1'b0;
					weight_read_addr_sel 			= 2'b0;
					write_enable_sel				= 0;
					temp_count_var_A_count_select 	= 2'b0;
					temp_count_var_C_count_select 	= 2'b0;
					temp_count_z_select				= 2'b0;
					next_state 						= S1; 
				end
				else begin 
					next_state = S0;
					set_dut_ready    			    = 1'b1;
                    // sram_write_address_r            = -1;
                end
			end

			S1: begin 				//get rows and col of A,B,C
				set_dut_ready       = 1'b0;
				input_read_data  = tb__dut__sram_input_read_data;
				weight_read_data = tb__dut__sram_weight_read_data;
				rows_A 	         = input_read_data[31:16];
				column_A 		 = input_read_data[15:0];
				rows_B  		 = weight_read_data[31:16];
				column_B		 = weight_read_data[15:0];
				column_C 		 = column_B;
				rows_C   		 = rows_A;

				// set_dut_ready    			  = 1'b0;
				temp_count_var_A_count_select = 2'b00;
				temp_count_var_C_count_select = 2'b10;
				temp_count_z_select 		  = 2'b0;
				input_read_addr_sel			  = 0;
				weight_read_addr_sel		  = 0;
				write_enable_sel 			  = 0;
				next_state 					  =	S2;
			end

			S2: begin  // update counters only
				set_dut_ready    			  = 1'b0;

				temp_count_var_A_count_select = 2'b10;
				temp_count_var_C_count_select = 2'b10; //00
				temp_count_z_select 		  = 0;
				input_read_addr_sel			  = 2'b01;
				weight_read_addr_sel		  = 2'b01;
				write_enable_sel 			  = 0;
		//		input_read_data  = tb__dut__sram_input_read_data;
		//		weight_read_data = tb__dut__sram_weight_read_data;
				next_state 					  =	S3;
			end

			S3: begin  // read from SRAM
				set_dut_ready    			  = 1'b0;

				temp_count_var_A_count_select = 2'b01;
				temp_count_var_C_count_select = 2'b10;
				input_read_addr_sel			  = 2'b01;
				weight_read_addr_sel		  = 2'b01;
				temp_count_z_select 		  = 2'b10;
				write_enable_sel 			  = 0;
				next_state 					  =	S4;
			end

			S4: begin
				set_dut_ready    			  = 1'b0;
				temp_count_var_A_count_select = 2'b10;
				temp_count_var_C_count_select = 2'b10;
				input_read_addr_sel			  = 2'b10;
				weight_read_addr_sel		  = 2'b10;//10
				temp_count_z_select 		  = 2'b01;
				//write_enable_sel 			  = 0;

				if(temp_count_var_A_count + 1 == column_B) begin
					write_enable_sel 			  = 1;
					if(temp_count_var_C_count == column_B -1) begin
						input_read_addr_sel			  = 2'b10;//01
						weight_read_addr_sel		  = 2'b10;
					end
					next_state = S5;
				end
				else begin
					write_enable_sel 			  = 0;
					next_state = S3;
				end
			end

			S5: begin // write into sram c
				set_dut_ready    			              = 1'b0;
				write_enable_sel 			              = 0;
				if((temp_count_var_C_count == column_B -1) && ~(sram_write_address_r == (column_C*rows_C - 1))) begin // next row
					temp_count_var_A_count_select     = 2'b11;
					temp_count_var_C_count_select     = 2'b0;

					input_read_addr_sel               = 2'b01;
					weight_read_addr_sel              = 2'b01;
					
					input_read_addr_sel	        		  = 2'b01;//01
					weight_read_addr_sel		          = 2'b01;//10
					temp_count_z_select 		          = 2'b00;//01
					next_state 					 			        = S3;
				end
				else if(sram_write_address_r == (column_C*rows_C - 1)) begin
					next_state = S6;
				end
				else if (sram_read_address_B == rows_B*column_B) begin //when 1 C calculation is done, start A and B from start
					input_read_addr_sel           = 2'b01;
					weight_read_addr_sel          = 2'b01;
				  temp_count_var_A_count_select = 2'b0;
					temp_count_var_C_count_select = 2'b01;
					temp_count_z_select 		      = 2'b0;
					next_state                    = S5_5;
				end  //only a was changing, when added, 397 to 403, new state S5_5
				else begin //if(temp_count_var_C_count < column_B) begin // in same row
					temp_count_var_A_count_select = 2'b0;
					temp_count_var_C_count_select = 2'b01;
					input_read_addr_sel			  = 2'b11;
					weight_read_addr_sel		  = 2'b01;
					
					if(temp_count_var_C_count == column_B - 1)
						temp_count_z_select 		  = 2'b10;
					else 
					   temp_count_z_select 		  = 2'b0;
						 next_state 					    = S5_5;// hold state
				end
			end
			S5_5: begin
				input_read_addr_sel			  = 2'b10;
				weight_read_addr_sel		  = 2'b10;
				next_state = S3;
			end
			S6: begin
				set_dut_ready = 1;
				next_state = S0;
			end
			default: next_state = S0;
		endcase
	end

	always@(posedge clk) begin
		write_enable_sel_r <= write_enable_sel;
	end

	assign dut__tb__sram_result_write_data =  mac_result_z;//sram_result_write_data_reg ;
	assign dut__tb__sram_input_read_address = sram_read_address_A;
	assign dut__tb__sram_weight_read_address = sram_read_address_B;

	assign dut__tb__sram_result_write_address = sram_write_address_r;
	assign dut__tb__sram_result_write_enable = write_enable_sel_r;
	
	assign dut__tb__sram_input_write_enable = 0;
	assign dut__tb__sram_input_write_address = 0;
	assign dut__tb__sram_input_write_data = 0;
	assign dut__tb__sram_weight_write_enable = 0;
	assign dut__tb__sram_weight_write_address = 0;
	assign dut__tb__sram_weight_write_data = 0;
	assign dut__tb__sram_result_read_address = 0;

// DUT ready handshake logic
	always @(posedge clk) begin : proc_compute_complete
		if(!reset_n) begin
			compute_complete <= 0;
		end else begin
			compute_complete <= (set_dut_ready == 1'b1) ? 1'b1 : 1'b0;
		end
	end
	assign dut_ready = compute_complete;


DW_fp_mac_inst
  FP_MAC (
  .inst_a(tb__dut__sram_input_read_data),
  .inst_b(tb__dut__sram_weight_read_data),
  .inst_c(accum_result),
  .inst_rnd(3'b0),
  .z_inst(mac_result_z),
  .status_inst()
);

endmodule

module DW_fp_mac_inst #(
  parameter inst_sig_width = 23,
  parameter inst_exp_width = 8,
  parameter inst_ieee_compliance = 0 // These need to be fixed to decrease error
) (
  input wire [inst_sig_width+inst_exp_width : 0] inst_a,
  input wire [inst_sig_width+inst_exp_width : 0] inst_b,
  input wire [inst_sig_width+inst_exp_width : 0] inst_c,
  input wire [2 : 0] inst_rnd,
  output wire [inst_sig_width+inst_exp_width : 0] z_inst,
  output wire [7 : 0] status_inst
);
	
  // Instance of DW_fp_mac
  DW_fp_mac #(inst_sig_width, inst_exp_width, inst_ieee_compliance) U1 (
    .a(inst_a),
    .b(inst_b),
    .c(inst_c),
    .rnd(inst_rnd),
    .z(z_inst),
    .status(status_inst)
		// status_inst doesn't need to be used at all. inst_rnd is used to set a value for rounding, but we don't need to use it for this, so you can set this to 0.
  );

endmodule: DW_fp_mac_inst