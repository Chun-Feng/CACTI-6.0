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


#ifndef _cacti_interface
#define _cacti_interface
/*dt: currently we assume that everything is counted in bytes*/
#define CHUNKSIZE 8

typedef struct {
    double dynamic;
    double leakage;
} powerComponents;


typedef struct {
    powerComponents readOp;
    powerComponents writeOp;
} powerDef;

enum wire_type {
    Global /* gloabl wires with repeaters */, 
    Global_10 /* 10% delay penalty */,
    Global_20 /* 20% delay penalty */,
    Global_30 /* 30% delay penalty */,
    Low_swing /* differential low power wires with high area overhead */,
    Semi_global /* midlevel wires with repeaters*/, 
    Transmission /* tranmission lines with high area overhead */, 
    Optical /* optical wires */, 
    Invalid_wtype
};


typedef struct {
    double height;
    double width;
    double area;
    double scaled_area;
}area_type;

typedef struct {
        area_type dataarray_area,datapredecode_area;
        area_type datacolmuxpredecode_area,datacolmuxpostdecode_area;
        area_type datawritesig_area;
        area_type tagarray_area,tagpredecode_area;
        area_type tagcolmuxpredecode_area,tagcolmuxpostdecode_area;
        area_type tagoutdrvdecode_area;
        area_type tagoutdrvsig_area;
        double totalarea, subbankarea;
        double total_dataarea;
        double total_tagarea;
        double max_efficiency, efficiency;
        double max_aspect_ratio_total, aspect_ratio_total;
        double perc_data, perc_tag, perc_cont, sub_eff, total_eff;
}arearesult_type;

typedef struct  input_parameter{
    long int cache_size;
    int number_of_sets;
    int block_size;
    int associativity;
    int tag_associativity, data_associativity;
    int rw_ports;
    int excl_read_ports;
    int excl_write_ports;
    int single_ended_read_ports;
    int uca_banks; /* UCA banks */
    int fully_assoc;
    double fudgefactor;
    double tech_size;
    double vdd_periph_global;
    int sequential_access;
    int access_mode;
    int fast_access;
    int force_tag; /* to model user defined tag size (set this to 1) */
    int tag_size;
    int output_width;
    int pure_sram;
    int dram;
    int nuca;
    double max_freq;
    int delay_wt;
    int dynamic_power_wt;
    int leakage_power_wt;
    int cycle_time_wt;
    int area_wt;
    int delay_dev;
    int dynamic_power_dev;
    int leakage_power_dev;
    int cycle_time_dev;
    int area_dev;
  
    int delay_wt_nuca;
    int dynamic_power_wt_nuca;
    int leakage_power_wt_nuca;
    int cycle_time_wt_nuca;
    int area_wt_nuca;
    int delay_dev_nuca;
    int dynamic_power_dev_nuca;
    int leakage_power_dev_nuca;
    int cycle_time_dev_nuca;
    int area_dev_nuca;
    int temp;
    int nuca_bank_count;
    int force_nuca_bank;
    int force_wiretype;
    enum wire_type wire_inter_mats;
    int cores;
    int cache_level; // 0 - L2, 1 - L3
    int print_detail; // 1 - detailed output
    int ed; // optimize for ed (1) or ed2 (2) or none (0)
} input_params_t;

typedef struct {
   int subbanks;
   double access_time,cycle_time;
   double senseext_scale;
   powerDef total_power;
   int best_Ndwl,best_Ndbl, best_data_deg_bitline_muxing, best_Ndsam;
   double max_leakage_power, max_access_time, max_cycle_time, max_dynamic_power;
   double min_leakage_power, min_access_time, min_cycle_time, min_dynamic_power;
   double best_Nspd;
   int best_Ntwl,best_Ntbl, best_tag_deg_bitline_muxing, best_Ntsam;
   double best_Ntspd;
   int best_muxover;
   powerDef total_routing_power;
   powerDef total_power_without_routing, total_power_allbanks;
   double subbank_address_routing_delay;
   powerDef subbank_address_routing_power;
   double decoder_delay_data,decoder_delay_tag;
   powerDef decoder_power_data,decoder_power_tag;
   double dec_data_driver,dec_data_3to8,dec_data_inv;
   double dec_tag_driver,dec_tag_3to8,dec_tag_inv;
   double wordline_delay_data,wordline_delay_tag;
   powerDef wordline_power_data,wordline_power_tag;
   double bitline_delay_data,bitline_delay_tag;
   powerDef bitline_power_data,bitline_power_tag;
   double sense_amp_delay_data,sense_amp_delay_tag;
   powerDef sense_amp_power_data,sense_amp_power_tag;
   double total_out_driver_delay_data;
   powerDef total_out_driver_power_data;
   double compare_part_delay;
   double drive_mux_delay;
   double selb_delay;
   powerDef compare_part_power, drive_mux_power, selb_power;
   double data_output_delay;
   powerDef data_output_power;
   double drive_valid_delay;
   powerDef drive_valid_power;
   double precharge_delay;
   int data_nor_inputs;
   int tag_nor_inputs;
} result_type;


typedef struct{
    int Ndwl;
    int Ndbl;
    double Nspd;
    int deg_bitline_muxing;
    int Ndsam;
    enum wire_type wt;
    int number_activated_mats_horizontal_direction;
    double delay_route_to_bank;
    double delay_addr_din_horizontal_htree;
    double delay_addr_din_vertical_htree;
    double delay_row_predecode_driver_and_block;
    double delay_row_decoder;
    double delay_bitlines;
    double delay_sense_amp;
    double delay_subarray_output_driver;
    double delay_bit_mux_predecode_driver_and_block;
    double delay_bit_mux_decoder;
    double delay_senseamp_mux_predecode_driver_and_block;
    double delay_senseamp_mux_decoder;
    double delay_dout_vertical_htree;
    double delay_dout_horizontal_htree;
    double delay_comparator;
    double access_time;
    double cycle_time;
    double multisubbank_interleave_cycle_time;
    powerDef power_routing_to_bank;
    powerDef power_addr_horizontal_htree;
    powerDef power_datain_horizontal_htree;
    powerDef power_dataout_horizontal_htree;
    powerDef power_addr_vertical_htree;
    powerDef power_datain_vertical_htree;
    powerDef power_row_predecoder_drivers;
    powerDef power_row_predecoder_blocks;
    powerDef power_row_decoders;
    powerDef power_bit_mux_predecoder_drivers;
    powerDef power_bit_mux_predecoder_blocks;
    powerDef power_bit_mux_decoders;
    powerDef power_senseamp_mux_predecoder_drivers;
    powerDef power_senseamp_mux_predecoder_blocks;
    powerDef power_senseamp_mux_decoders;
    powerDef power_bitlines;
    powerDef power_sense_amps;
    powerDef power_output_drivers_at_subarray;
    powerDef power_dataout_vertical_htree;
    powerDef power_comparators;
    powerDef total_power;
    double area;
    double all_banks_height;
    double all_banks_width;
    double bank_height;
    double bank_width;
    double subarray_memory_cell_area_height;
    double subarray_memory_cell_area_width;
    double mat_height;
    double mat_width;
    double routing_area_height_within_bank;
    double routing_area_width_within_bank;
    double area_efficiency;
    double refresh_power;
    double dram_refresh_period;
    double dram_array_availability;

    double delay_before_subarray_output_driver;
    double delay_from_subarray_output_driver_to_output;
}results_mem_array;


typedef struct{
    results_mem_array tag_array;
    results_mem_array data_array;
    double access_time;
    double cycle_time;
    double area;
    double area_efficiency;
    powerDef power;
    input_params_t *params;
    double cache_ht;
    double cache_len;
    char file_n[100]; /* .fig/.png file name to output cache 
                         organization */
}final_results;

typedef struct uca_org{
    results_mem_array *tag_array;
    results_mem_array *data_array;
    double access_time;
    double cycle_time;
    double area;
    double area_efficiency;
    powerDef power;
    input_params_t *params;
    double cache_ht;
    double cache_len;
    char file_n[100]; /* .fig/.png file name to output cache 
                         organization */
} uca_org_t;

typedef struct pda_res{
    powerComponents power;
    area_type area_stats;
    double delay;
    double cycle_time;
}pda_res_t;

// The following global variables are temporarily used 
// to get case study results

/* linked list to store results */
typedef struct uca_res_lentry{
    uca_org_t fin_res;
    struct uca_res_lentry *next_lentry;
}uca_res_lentry_t;

enum router_type {
    Low_power /* less number of crossbar ports */,
    Normal,
    High_freq
};

typedef struct wire_stats {
    enum wire_type wt;
    pda_res_t wire_pda;
    /* full swing */
    double repeater_size;
    double repeater_spacing;
    double wire_spacing;
    double wire_width;
    double wire_length;
    /* low swing */
    pda_res_t transmitter;
    pda_res_t l_wire;
    pda_res_t sense_amp;
    int nsense; /* no. of sense amp connected to the wire */
}wire_stats_t;

typedef struct router_stats {
    pda_res_t router_pda;
    int flit_size;
    int buffer_ent;
    int vc_count;
    double max_cyc;
    double cycle_time;
    double simple_buffer_read;
    double simple_buffer_write;
    uca_org_t ubuffer;
//    results_mem_array buffer;
    pda_res_t crossbar;
    pda_res_t arbiter;
    enum router_type rt;
} router_stats_t;

typedef struct nuca_org {
    int size;
    /* area, power, access time, and cycle time stats */
    pda_res_t nuca_pda;
    pda_res_t bank_pda;
//    pda_res_t router_pda;
    pda_res_t wire_pda;
    wire_stats_t h_wire;
    wire_stats_t v_wire;
    router_stats_t router;
    /* for particular network configuration
     * calculated based on a cycle accurate
     * simulation Ref: CACTI 6 - Tech report
     */
    double contention;

    /* grid network stats */
    double avg_hops;
    int rows;
    int columns;
    input_params_t *params;
    int bank_count;
} nuca_org_t;

typedef struct nuca_res_lentry{
    nuca_org_t nres;
    struct nuca_res_lentry *next_lentry;
}nuca_res_lentry_t;

final_results cacti_interface(
        int cache_size,
        int line_size,
        int associativity,
        int rw_ports,
        int excl_read_ports,
        int excl_write_ports,
        int single_ended_read_ports,
        int banks,
        double tech_node,
        int output_width,
        int specific_tag,
        int tag_width,
        int access_mode,
        int pure_sram,
        int dram,
        int delay_wt,
        int dynamic_power_wt,
        int leakage_power_wt,
        int cycle_time_wt,
        int temp,
        int sram_cell_and_wordline_tech_flavor_in, 
        int periph_global_tech_flavor_in, 
        int interconnect_projection_type_in,
        int wire_inside_mat_type_in, 
        int wire_outside_mat_type_in, 
        int HTREE_NODES_AT_MAT_INTERVALS_in, 
        int VERTICAL_HTREE_WIRES_OVER_THE_ARRAY_in,
        int BROADCAST_ADDR_DATAIN_OVER_VERTICAL_HTREES_in, 
        double MIN_PERCENT_WITHIN_BEST_AREA_in, 
        double MIN_PERCENT_WITHIN_BEST_DELAY_in,
        int gl,
        int start,
        int contr_lat);
#endif
