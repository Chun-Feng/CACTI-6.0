/*------------------------------------------------------------
 *                              CACTI 6.0
 *         Copyright 2007 Hewlett-Packard Development Corporation
 *                         All Rights Reserved
 *
 * Permission to use, copy, and modify this software and its documentation is
 * hereby granted only under the following terms and conditions.  Both the
 * above copyright notice and this permission notice must appear in all copies
 * of the software, derivative works or modified versions, and any portions
 * thereof, and both notices must appear in supporting documentation.
 *
 * Users of this software agree to the terms and conditions set forth herein, and
 * hereby grant back to Hewlett-Packard Company and its affiliated companies ("HP")
 * a non-exclusive, unrestricted, royalty-free right and license under any changes, 
 * enhancements or extensions  made to the core functions of the software, including 
 * but not limited to those affording compatibility with other hardware or software
 * environments, but excluding applications which incorporate this software.
 * Users further agree to use their best efforts to return to HP any such changes,
 * enhancements or extensions that they make and inform HP of noteworthy uses of
 * this software.  Correspondence should be provided to HP at:
 *
 *                       Director of Intellectual Property Licensing
 *                       Office of Strategy and Technology
 *                       Hewlett-Packard Company
 *                       1501 Page Mill Road
 *                       Palo Alto, California  94304
 *
 * This software may be distributed (but not offered for sale or transferred
 * for compensation) to third parties, provided such third parties agree to
 * abide by the terms and conditions of this notice.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND HP DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS.   IN NO EVENT SHALL HP 
 * CORPORATION BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 * SOFTWARE.
 *------------------------------------------------------------*/

#ifndef _areadef
#define _areadef
#include "basic_circuit.h"
#include "leakage.h"

double area_all_datasubarrays;
double area_all_tagsubarrays;
double area_all_tagramcells;
double faarea_all_subarrays ;

double aspect_ratio_data;
double aspect_ratio_tag;
double aspect_ratio_subbank;
double aspect_ratio_total;

/*
was: #define Widthptondiff 4.0
now: is 4* feature size, taken from the Intel 65nm process paper
*/

//v4.1: Making all constants static variables. Initially these variables are based
//off 0.8 micron process values; later on in init_tech_params function of leakage.c 
//they are scaled to input tech node parameters

//#define Widthptondiff 3.2
double Widthptondiff;
/* 
was: #define Widthtrack    3.2
now: 3.2*FEATURESIZE, i.e. 3.2*0.8
*/
//#define Widthtrack    (3.2*0.8)
double Widthtrack;
//#define Widthcontact  1.6
double Widthcontact;
//#define Wpoly         0.8
double Wpoly;
//#define ptocontact    0.4
double ptocontact;
//#define stitch_ramv   6.0 
double stitch_ramv;
//#define BitHeight16x2 33.6
//#define BitHeight1x1 (2*7.746*0.8) /* see below */
double BitHeight1x1;
//#define stitch_ramh   12.0
double stitch_ramh;
//#define BitWidth16x2  192.8
//#define BitWidth1x1	  (7.746*0.8) 
double BitWidth1x1;
/* dt: Assume that each 6-T SRAM cell is 120F^2 and has an aspect ratio of 1(width) to 2(height), than the width is 2*sqrt(60)*F */
//#define WidthNOR1     11.6
double WidthNOR1;
//#define WidthNOR2     13.6
double WidthNOR2;
//#define WidthNOR3     20.8
double WidthNOR3;
//#define WidthNOR4     28.8
double WidthNOR4;
//#define WidthNOR5     34.4
double WidthNOR5;
//#define WidthNOR6     41.6
double WidthNOR6;
//#define Predec_height1    140.8
double Predec_height1;
//#define Predec_width1     270.4
double Predec_width1;
//#define Predec_height2    140.8
double Predec_height2;
//#define Predec_width2     539.2
double Predec_width2;
//#define Predec_height3    281.6
double Predec_height3;    
//#define Predec_width3     584.0
double Predec_width3;
//#define Predec_height4    281.6
double Predec_height4;  
//#define Predec_width4     628.8
double Predec_width4; 
//#define Predec_height5    422.4
double Predec_height5; 
//#define Predec_width5     673.6
double Predec_width5;
//#define Predec_height6    422.4
double Predec_height6;
//#define Predec_width6     718.4
double Predec_width6;
//#define Wwrite		  1.2
double Wwrite;
//#define SenseampHeight    152.0
double SenseampHeight;
//#define OutdriveHeight	  200.0
double OutdriveHeight;
//#define FAOutdriveHeight  229.2
double FAOutdriveHeight;
//#define FArowWidth	  382.8
double FArowWidth;
//#define CAM2x2Height_1p	  48.8
double CAM2x2Height_1p;
//#define CAM2x2Width_1p	  44.8
double CAM2x2Width_1p;
//#define CAM2x2Height_2p   80.8 
double CAM2x2Height_2p;   
//#define CAM2x2Width_2p    76.8
double CAM2x2Width_2p;
//#define DatainvHeight     25.6
double DatainvHeight;
//#define Wbitdropv 	  30.0
double Wbitdropv;
//#define decNandWidth      34.4
double decNandWidth;
//#define FArowNANDWidth    71.2
double FArowNANDWidth;
//#define FArowNOR_INVWidth 28.0  
double FArowNOR_INVWidth;

//#define FAHeightIncrPer_first_rw_or_w_port 16.0
double FAHeightIncrPer_first_rw_or_w_port;
//#define FAHeightIncrPer_later_rw_or_w_port 16.0
double FAHeightIncrPer_later_rw_or_w_port;
//#define FAHeightIncrPer_first_r_port       12.0
double FAHeightIncrPer_first_r_port;
//#define FAHeightIncrPer_later_r_port       12.0
double FAHeightIncrPer_later_r_port;
//#define FAWidthIncrPer_first_rw_or_w_port  16.0 
double FAWidthIncrPer_first_rw_or_w_port;
//#define FAWidthIncrPer_later_rw_or_w_port  9.6
double FAWidthIncrPer_later_rw_or_w_port;
//#define FAWidthIncrPer_first_r_port        12.0
double FAWidthIncrPer_first_r_port;
//#define FAWidthIncrPer_later_r_port        9.6
double FAWidthIncrPer_later_r_port;

#define tracks_precharge_p    12
#define tracks_precharge_nx2   5 
#define tracks_outdrvselinv_p  3
#define tracks_outdrvfanand_p  6  

#define CONVERT_TO_MMSQUARE 1.0/1000000.0

void
area_subbanked (int baddr,int b0,int RWP,int ERP,int EWP,int Ndbl,int Ndwl,double Nspd,int Ntbl,int Ntwl,int Ntspd,
		double NSubbanks,input_params_t *parameters,area_type *result_subbanked,arearesult_type *result);




#define MAX_NUMBER_GATES_STAGE 20
#define MAX_NUMBER_HTREE_NODES 20

typedef struct{
	int flag_driver_exists;
	int flag_driving_decoder_output;
	int number_input_addr_bits;
	int number_gates_nand2_path;
	int number_gates_nand3_path;
	int min_number_gates;
	int number_parallel_instances_driving_1_nand2_load;
	int number_parallel_instances_driving_2_nand2_load;
	int number_parallel_instances_driving_4_nand2_load;
	int number_parallel_instances_driving_2_nand3_load;
	int number_parallel_instances_driving_8_nand3_load;
	int number_parallel_instances_nand3_path;
	double c_load_nand2_path_predecode_block_driver_output;
	double c_load_nand3_path_predecode_block_driver_output;
	double r_load_nand2_path_predecode_block_driver_output;
	double r_load_nand3_path_predecode_block_driver_output;
	double width_nand2_path_n[MAX_NUMBER_GATES_STAGE];
	double width_nand2_path_p[MAX_NUMBER_GATES_STAGE];
	double width_nand3_path_n[MAX_NUMBER_GATES_STAGE];
	double width_nand3_path_p[MAX_NUMBER_GATES_STAGE];
	double delay_nand2_path;
	double delay_nand3_path;
	powerDef power_nand2_path;
	powerDef power_nand3_path;
} predecoder_block_driver;


typedef struct{
	int flag_block_exists;
	int number_input_addr_bits;
	double c_load_predecoder_block_output;
	double r_wire_predecoder_block_output;
	int branch_effort_nand2_gate_output;
	int branch_effort_nand3_gate_output;
	int flag_two_unique_paths;
	int flag_second_level_gate;
	int number_inputs_first_level_gate;
	int number_gates_first_level_nand2_path;
	int number_gates_first_level_nand3_path;
	int number_gates_second_level;
	int min_number_gates_first_level;
	int min_number_gates_second_level;
	int number_first_level_parallel_instances_nand2;
	int number_first_level_parallel_instances_nand3;
	int number_second_level_parallel_instances;
	double width_first_level_nand2_path_n[MAX_NUMBER_GATES_STAGE];
	double width_first_level_nand2_path_p[MAX_NUMBER_GATES_STAGE];
	double width_first_level_nand3_path_n[MAX_NUMBER_GATES_STAGE];
	double width_first_level_nand3_path_p[MAX_NUMBER_GATES_STAGE];
	double width_second_level_n[MAX_NUMBER_GATES_STAGE];
	double width_second_level_p[MAX_NUMBER_GATES_STAGE];
	double delay_nand2_path;
	double delay_nand3_path;
	powerDef power_nand2_path;
	powerDef power_nand3_path;
	powerDef power_second_level;
} predecoder_block;


typedef struct{
	int flag_decoder_exists;
	int number_input_signals;
	double c_load_decoder_output;
	double r_wire_decoder_output;
	int number_gates;
	int min_number_gates;
	double width_decoder_n[MAX_NUMBER_GATES_STAGE];
	double width_decoder_p[MAX_NUMBER_GATES_STAGE];
	double delay;
	powerDef power;
} decoder;


typedef struct{
	int number_gates;
	int min_number_gates;
	double width_n[MAX_NUMBER_GATES_STAGE];
	double width_p[MAX_NUMBER_GATES_STAGE];
	double c_gate_load;
	double c_wire_load;
	double r_wire_load;
	double f_wire;
	double delay;
	powerDef power;
} addr_datain_htree_node;

typedef struct{
	int number_gates_nand_driving_inv;
	int number_gates_inv_driving_inv;
	int number_gates_inv_driving_nand;
	int number_gates_nand_driving_inv_final_seg;
	int number_gates_tristate_driver_driving_inv;
	int min_number_gates;
	int min_number_gates_tristate_driver;
	double nand_driving_inv_width_n[MAX_NUMBER_GATES_STAGE];
	double nand_driving_inv_width_p[MAX_NUMBER_GATES_STAGE];
	double inv_driving_inv_width_n[MAX_NUMBER_GATES_STAGE];
	double inv_driving_inv_width_p[MAX_NUMBER_GATES_STAGE];
	double inv_driving_nand_width_n[MAX_NUMBER_GATES_STAGE];
	double inv_driving_nand_width_p[MAX_NUMBER_GATES_STAGE];
	double nand_driving_inv_final_seg_width_n[MAX_NUMBER_GATES_STAGE];
	double nand_driving_inv_final_seg_width_p[MAX_NUMBER_GATES_STAGE];
	double tristate_driver_driving_inv_width_n[MAX_NUMBER_GATES_STAGE];
	double tristate_driver_driving_inv_width_p[MAX_NUMBER_GATES_STAGE];
	double tristate_driver_driving_inv_width_nor2_n;
	double tristate_driver_driving_inv_width_nor2_p;
	double c_gate_load_nand_driving_inv;
	double c_gate_load_inv_driving_inv;
	double c_gate_load_inv_driving_nand;
	double c_gate_load_tristate_driver_driving_inv;
	double c_gate_load_nand_driving_inv_final_seg;
	double c_wire_load;
	double r_wire_load;
	double c_wire_load_final_seg;
	double r_wire_load_final_seg;
	double delay_nand_driving_inv;
	double delay_inv_driving_inv;
	double delay_inv_driving_nand;
	double delay_nand_driving_inv_final_seg;
	double delay_tristate_driver_driving_inv;
	powerDef power_nand_driving_inv;
	powerDef power_inv_driving_inv;
	powerDef power_inv_driving_nand;
	powerDef power_nand_driving_inv_final_seg;
	powerDef power_tristate_driver_driving_inv;
	double area_nand2_driving_inv;
	double area_inv_driving_inv;
	double area_inv_driving_nand2;
	double area_nand2_driving_inv_final_seg;
	double area_tristate_driver_driving_inv;
} addr_datain_htree_at_mat_interval;



typedef struct{
	int number_gates;
	int min_number_gates;
	double width_n[MAX_NUMBER_GATES_STAGE];
	double width_p[MAX_NUMBER_GATES_STAGE];
	double c_gate_load;
	double c_wire_load;
	double r_wire_load;
	double delay;
	powerDef power;
} driver;

typedef struct{
	int number_gates;
	int min_number_gates;
	double width_n[MAX_NUMBER_GATES_STAGE];
	double width_p[MAX_NUMBER_GATES_STAGE];
	double width_nor2_n;
	double width_nor2_p;
	double c_gate_load;
	double c_wire_load;
	double r_wire_load;
	double f_wire;
	double delay;
	powerDef power;
} dataout_htree_node;

typedef struct {
	double length;
	double width_n[MAX_NUMBER_GATES_STAGE];
	double width_p[MAX_NUMBER_GATES_STAGE];
	double delay;
	powerDef power;
	int number_gates;
	int minimum_number_gates;
	double optimal_repeater_size;
	int number_repeater_gates;
	double c_gate_load;
}point_to_point_interconnect_segment;


//double length_wire_htree_node[MAX_NUMBER_GATES_STAGE];

//predecoder_block row_predec_blk_1, row_predec_blk_2, bit_mux_predec_blk_1, bit_mux_predec_blk_2,
//	senseamp_mux_predec_blk_1, senseamp_mux_predec_blk_2;

//decoder row_dec, bit_mux_dec, senseamp_mux_dec;
//row_predec_blk_driver_1
//predecoder_block_driver  row_predec_blk_driver_2, bit_mux_predec_blk_driver_1,
//	bit_mux_predec_blk_driver_2, senseamp_mux_predec_blk_driver_1, senseamp_mux_predec_blk_driver_2;

//addr_datain_htree_node addr_din_htree_node[MAX_NUMBER_HTREE_NODES], 
//horizontal_addr_din_htree_node[MAX_NUMBER_HTREE_NODES];
//dataout_htree_node dout_htree_node[MAX_NUMBER_HTREE_NODES], subarray_output_htree_node;
//dataout_htree_node  subarray_output_htree_node;

//powerDef tot_power,  tot_power_row_predecode_block_drivers,
//tot_power_bit_mux_predecode_block_drivers, tot_power_senseamp_mux_predecode_block_drivers,
//tot_power_row_predecode_blocks, tot_power_bit_mux_predecode_blocks, 
//tot_power_senseamp_mux_predecode_blocks, tot_power_row_decoders, 
//tot_power_bit_mux_decoders, tot_power_senseamp_mux_decoders;


#define NAND2_LEAK_STACK_FACTOR 0.2
#define NOR2_LEAK_STACK_FACTOR 0.2
#define NAND3_LEAK_STACK_FACTOR 0.2
#define MAX_NUMBER_ARRAY_PARTITIONS 100000


area_type subarraymem_area(int number_rows_subarray, int number_cols_subarray, int number_subarrays, int RWP,int ERP,int EWP,int 
				  NSER);

area_type area_mat(int is_tag, int number_rows_subarray, int number_cols_subarray, 
        int number_subarrays, int deg_bitline_muxing, 
        int deg_senseamp_muxing_non_associativity, 
		int Ndsam,
		int number_addr_bits_mat,
		int number_datain_bits_mat,
		int number_dataout_bits_mat,
        int number_way_select_signals_mat,
		input_params_t *parameters,
        predecoder_block* row_predec_blk_1, predecoder_block* row_predec_blk_2, 
        predecoder_block* bit_mux_predec_blk_1, predecoder_block* bit_mux_predec_blk_2, 
        predecoder_block* senseamp_mux_predec_blk_1, predecoder_block* senseamp_mux_predec_blk_2,
		predecoder_block* dummy_way_select_predec_blk_1,
        decoder *row_dec, decoder *bit_mux_dec, decoder *senseamp_mux_dec,
        predecoder_block_driver *row_predec_blk_driver_1, 
        predecoder_block_driver *row_predec_blk_driver_2,
        predecoder_block_driver *bit_mux_predec_blk_driver_1,
        predecoder_block_driver *bit_mux_predec_blk_driver_2,
        predecoder_block_driver *senseamp_mux_predec_blk_driver_1,
        predecoder_block_driver *senseamp_mux_predec_blk_driver_2, 
		predecoder_block_driver *way_select_driver_1, 
        dataout_htree_node *subarray_output_htree_node);

area_type area_single_bank(int number_rows_subarray, int is_tag, input_params_t *parameters, arearesult_type *result, 
		    int number_horizontal_htree_nodes, int number_vertical_htree_nodes, 
			int number_mats_horizontal_direction, int number_mats_vertical_direction,
			int number_addr_bits_mat, int number_way_select_signals_mat, int tagbits, 
			int number_datain_bits_mat, int number_dataout_bits_mat,
			addr_datain_htree_node *ptr_horizontal_addr_din_htree_node,
			addr_datain_htree_node *ptr_addr_din_htree_node,
			dataout_htree_node *ptr_dout_htree_node, point_to_point_interconnect_segment 
			*ptr_horizontal_addr_intcnt_segment_within_bank, point_to_point_interconnect_segment
			*ptr_horizontal_datain_intcnt_segment_within_bank,
            predecoder_block* row_predec_blk_1, predecoder_block* row_predec_blk_2, 
            predecoder_block* bit_mux_predec_blk_1, predecoder_block* bit_mux_predec_blk_2, 
            predecoder_block* senseamp_mux_predec_blk_1, predecoder_block* senseamp_mux_predec_blk_2,
            decoder *row_dec, decoder *bit_mux_dec, decoder *senseamp_mux_dec,
            predecoder_block_driver *row_predec_blk_driver_1, 
            predecoder_block_driver *row_predec_blk_driver_2,
            predecoder_block_driver *bit_mux_predec_blk_driver_1,
            predecoder_block_driver *bit_mux_predec_blk_driver_2,
            predecoder_block_driver *senseamp_mux_predec_blk_driver_1,
            predecoder_block_driver *senseamp_mux_predec_blk_driver_2,
            powerDef *tot_power, powerDef * tot_power_row_predecode_block_drivers,
            powerDef *tot_power_bit_mux_predecode_block_drivers,
            powerDef *tot_power_senseamp_mux_predecode_block_drivers,
            powerDef *tot_power_row_predecode_blocks, powerDef *tot_power_bit_mux_predecode_blocks,
            powerDef *tot_power_senseamp_mux_predecode_blocks,
            powerDef *tot_power_row_decoders, powerDef *tot_power_bit_mux_decoders,
            powerDef *tot_power_senseamp_mux_decoders);

area_type area_all_banks(int uca_banks, double bank_height, double bank_width, 
						 int number_bits_routed_to_bank, double *length_htree_route_to_bank,
						 int number_mats_vertical_direction);


#define MAX_STORED_SOLUTIONS 500

typedef struct {
	int tag_array_index;
	int data_array_index;
	double access_time;
	double cycle_time;
	double area;
	double efficiency;
	powerDef total_power;
} solution;


typedef struct {
	int Ndwl;
	int Ndbl;
	double Nspd;
	int deg_bitline_muxing; /* Ndcm */
	int Ndsam;
	double access_time;
	double cycle_time;
	double multisubbank_interleave_cycle_time;
	double area_ram_cells;
	double area;
	powerDef power;
	double delay_senseamp_mux_decoder;
	double delay_before_subarray_output_driver;
	double delay_from_subarray_output_driver_to_output;
}mem_array;


#define MAX_NUMBER_GATES_HORIZONTAL_ADDR_INTERCONNECT_SEGMENT 1000
//mem_array tag_arr[MAX_NUMBER_ARRAY_PARTITIONS];
//mem_array data_arr[MAX_NUMBER_ARRAY_PARTITIONS];
//solution best_solution2[1000000];
#endif

double HTREE_NODES_AT_MAT_INTERVALS;
double VERTICAL_HTREE_WIRES_OVER_THE_ARRAY;
double BROADCAST_ADDR_DATAIN_OVER_VERTICAL_HTREES;

addr_datain_htree_at_mat_interval horizontal_addr_datain_htree_at_mat_interval, 
vertical_addr_datain_htree_at_mat_interval;

