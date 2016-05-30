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

#include "areadef.h"

area_type
subarraymem_area(int number_rows_subarray, int number_cols_subarray, int number_subarrays, int RWP,int ERP,int EWP,int 
				  NSER)	/* returns area of subarray */
{
  area_type memarea;
 
  memarea.height = height_cell * number_rows_subarray;
  memarea.width = width_cell * number_cols_subarray;
  memarea.area = memarea.height * memarea.width;
  return (memarea);

}

//Function width_transistor_after_folding returns the width of a transistor after folding.
//Inputs: Input width of the transistor, maximum width that a transistor can have before
//it gets folded 
//Output: Width of the folded transistor
double width_transistor_after_folding(double input_width, double threshold_folding_width)
{
	double total_diff_width, width_poly, poly_to_poly_spacing_with_in_between_contact;
	int number_folded_transistors;
	if(input_width <= 0){
		return(0);
	}
	number_folded_transistors = (int) (ceil(input_width / threshold_folding_width));
	poly_to_poly_spacing_with_in_between_contact = WIDTHPOLYCONTACT + 2 * SPACINGPOLYTOCONTACT;
	width_poly = FEATURESIZE;
    total_diff_width = number_folded_transistors * width_poly + (number_folded_transistors + 1) * 
		poly_to_poly_spacing_with_in_between_contact;
	return(total_diff_width);
}


//Function area_sense_amplifier computes the area of a sense amplifier.
//Inputs: pointer to the sense amp data structure, pitch at which sense amp has to be laid out
//Output: Area of the sense amplifier
area_type area_sense_amplifier(double width_latch_pmos, double width_latch_nmos,
							   double width_enable, double width_iso, double pitch_sense_amp)
{
	area_type sense_amp;
	double height_pmos_transistors, height_nmos_transistors;

	//First compute the height occupied by all PMOS transistors. Calculate the height of each PMOS transistor
	//after folding.  

	height_pmos_transistors = 
		2 * width_transistor_after_folding(width_latch_pmos, pitch_sense_amp) + 
		width_transistor_after_folding(width_iso, pitch_sense_amp) +
		2 * MIN_GAP_BET_SAME_TYPE_DIFFS;

	//Now compute the height occupied by all NMOS transistors
	height_nmos_transistors = 
		2 * width_transistor_after_folding(width_latch_nmos, pitch_sense_amp) +
		width_transistor_after_folding(width_enable, pitch_sense_amp) +
		2 * MIN_GAP_BET_SAME_TYPE_DIFFS;

	//Calculate total height by considering gap between the p and n diffusion areas. 
	sense_amp.height = height_pmos_transistors + height_nmos_transistors + MIN_GAP_BET_P_AND_N_DIFFS;
	sense_amp.width = pitch_sense_amp;
	sense_amp.area = sense_amp.height * sense_amp.width;
	return(sense_amp);
}




//v5.0: This function calculates the area occupied by certain simple CMOS gates.
//Assumptions: The gate area is calculated using simple equations. The desired height of each gate
//is taken as input. The max width of PMOS and NMOS transistors that can fit within this height is 
//calculated. If the input width of transistors exceeds these max width, the transistors are
//folded and the area is then calculated. Calculation of the area makes use of simple design rules 
//like width of poly contact, minimum poly-to-poly spacing, max height of n/p diffusion etc. 
//BOUNDARY CASES: CHECK WHAT HAPPENS FOR 1-INPUT NOR/NAND GATE

/*area_type gatearea(int gatetype, int numberofinputs, double widthpmos, double widthnmos, 
				   double height_gate)
{
	double height_transistor_region, ratio_p_to_n;
	double width_folded_pmos, width_folded_nmos;
	double poly_to_poly_spacing_with_in_between_contact;
	int number_folded_pmos, number_folded_nmos;
	double total_ndiff_width, total_pdiff_width;
	area_type gate;
	
	height_transistor_region = height_gate - 2 * HPOWERRAIL;
	ratio_p_to_n = widthpmos / (widthpmos + widthnmos);
	width_folded_pmos = ratio_p_to_n * (height_transistor_region - MIN_GAP_BET_P_AND_N_DIFFS);
	width_folded_nmos = (1 - ratio_p_to_n) * (height_transistor_region - MIN_GAP_BET_P_AND_N_DIFFS);
	number_folded_pmos = (int) (ceil(widthpmos / width_folded_pmos));
    number_folded_nmos = (int) (ceil(widthnmos / width_folded_nmos));
	poly_to_poly_spacing_with_in_between_contact = WIDTHPOLYCONTACT + 2 * SPACINGPOLYTOCONTACT;
	if((number_folded_pmos == 0) && (number_folded_nmos == 0)){
		gate.width = 0;
		gate.height = 0;
		return(gate);
	}

	switch (gatetype)
    {
		case INV:
				//diff width = 
				total_ndiff_width = number_folded_nmos * FEATURESIZE +
									(number_folded_nmos + 1) * poly_to_poly_spacing_with_in_between_contact;
				total_pdiff_width = number_folded_pmos * FEATURESIZE +
									(number_folded_pmos + 1) * poly_to_poly_spacing_with_in_between_contact;
				//The source diffusions connect to VDD/GND rails using multiple contacts typically.
				//Presently assuming that only one contact contributes to the total p/n diffusion width
				break;

		case NOR:
				if(number_folded_nmos == 1){
					total_ndiff_width = numberofinputs * FEATURESIZE + 
						(numberofinputs + 1) * poly_to_poly_spacing_with_in_between_contact;
				}
				else{
						total_ndiff_width = numberofinputs * number_folded_nmos * FEATURESIZE +
							(numberofinputs * number_folded_nmos + 1) * poly_to_poly_spacing_with_in_between_contact;
				}

				if(number_folded_pmos == 1){
					total_pdiff_width = numberofinputs * FEATURESIZE + 
						(numberofinputs - 2) * SPACINGPOLYTOPOLY +
						2 * poly_to_poly_spacing_with_in_between_contact; 
				}
				else{
					total_pdiff_width = numberofinputs * number_folded_pmos * FEATURESIZE +
							(numberofinputs * number_folded_pmos + 1) * poly_to_poly_spacing_with_in_between_contact;
				}		
				break;

		case NAND:
				if(number_folded_pmos == 1){
					total_pdiff_width = numberofinputs * FEATURESIZE + 
						(numberofinputs + 1) * poly_to_poly_spacing_with_in_between_contact;
				}
				else{
						total_pdiff_width = numberofinputs * number_folded_pmos * FEATURESIZE +
							(numberofinputs * number_folded_pmos + 1) * poly_to_poly_spacing_with_in_between_contact;
				}

				if(number_folded_nmos == 1){
					total_ndiff_width = numberofinputs * FEATURESIZE + 
						(numberofinputs - 2) * SPACINGPOLYTOPOLY +
						2 * poly_to_poly_spacing_with_in_between_contact; 
				}
				else{
					total_ndiff_width = numberofinputs * number_folded_nmos * FEATURESIZE +
							(numberofinputs * number_folded_nmos + 1) * poly_to_poly_spacing_with_in_between_contact;
				}
				break;
		default:
				 printf ("Unknown gate type: %d", gatetype);
				 exit(1);
     }
	gate.width = MAX(total_ndiff_width, total_pdiff_width);
	if(width_folded_nmos > widthnmos){//means that the height of the gate can 
		//be made smaller than the input height specified, so calculate the height of the gate.
		gate.height = widthnmos + widthpmos + MIN_GAP_BET_P_AND_N_DIFFS + 2 * HPOWERRAIL;
	}
	else{
		gate.height = height_gate;
	}
	gate.area = gate.width * gate.height;
 return(gate);
}*/

double width_diffusion(int numberofinputs, int number_folded_transistors)
{
	double total_diff_width, width_poly;

	width_poly = FEATURESIZE;
	total_diff_width = 2 * (WIDTHPOLYCONTACT + 2 * SPACINGPOLYTOCONTACT) + 
				numberofinputs * width_poly + (numberofinputs - 1) * SPACINGPOLYTOPOLY;
	if(number_folded_transistors > 1){
		total_diff_width += (number_folded_transistors - 1) * (WIDTHPOLYCONTACT + 2 * SPACINGPOLYTOCONTACT +
			numberofinputs * width_poly + (numberofinputs - 1) * SPACINGPOLYTOPOLY);
	}
	return(total_diff_width);
}


area_type gatearea(int gatetype, int numberofinputs, double widthpmos, double widthnmos, 
				   double height_gate)
{
	double height_transistor_region, ratio_p_to_n;
	double width_folded_pmos, width_folded_nmos;
	double poly_to_poly_spacing_with_in_between_contact;
	int number_folded_pmos, number_folded_nmos;
	double total_ndiff_width, total_pdiff_width;
	area_type gate;
	
	height_transistor_region = height_gate - 2 * HPOWERRAIL;
	ratio_p_to_n = widthpmos / (widthpmos + widthnmos);
	width_folded_pmos = ratio_p_to_n * (height_transistor_region - MIN_GAP_BET_P_AND_N_DIFFS);
	width_folded_nmos = (1 - ratio_p_to_n) * (height_transistor_region - MIN_GAP_BET_P_AND_N_DIFFS);
	number_folded_pmos = (int) (ceil(widthpmos / width_folded_pmos));
    number_folded_nmos = (int) (ceil(widthnmos / width_folded_nmos));
	poly_to_poly_spacing_with_in_between_contact = WIDTHPOLYCONTACT + 2 * SPACINGPOLYTOCONTACT;
	if((number_folded_pmos == 0) && (number_folded_nmos == 0)){
		gate.width = 0;
		gate.height = 0;
		return(gate);
	}

	switch (gatetype)
    {
		case INV:
				total_ndiff_width = width_diffusion(1, number_folded_nmos);
				total_pdiff_width = width_diffusion(1, number_folded_pmos);
				break;

		case NOR:
				total_ndiff_width = width_diffusion(1, number_folded_nmos);
				total_pdiff_width = width_diffusion(numberofinputs, number_folded_pmos);
				break;

		case NAND:
				total_ndiff_width = width_diffusion(numberofinputs, number_folded_nmos);
				total_pdiff_width = width_diffusion(1, number_folded_pmos);
				break;
		default:
				 printf ("Unknown gate type: %d", gatetype);
				 exit(1);
     }
	gate.width = MAX(total_ndiff_width, total_pdiff_width);
	if(width_folded_nmos > widthnmos){//means that the height of the gate can 
		//be made smaller than the input height specified, so calculate the height of the gate.
		gate.height = widthnmos + widthpmos + MIN_GAP_BET_P_AND_N_DIFFS + 2 * HPOWERRAIL;
	}
	else{
		gate.height = height_gate;
	}
	gate.area = gate.width * gate.height;
 return(gate);
}



//v5.0: Making the area model sensitive to the the widths that have been computed by the 
//delay model. 

area_type area_predecoder_block_driver(predecoder_block_driver *ptr_predec_blk_driver,
									   predecoder_block *ptr_predec_blk)
{
	int i;
	area_type predecoder_block_driver, inv;
	double cumulative_area_nand2_path, cumulative_area_nand3_path;
	double leakage_nand2_path, leakage_nand3_path;

	predecoder_block_driver.width = 0;
	predecoder_block_driver.height = 0;
	predecoder_block_driver.area = 0;

	cumulative_area_nand2_path = 0;
	cumulative_area_nand3_path = 0;
	leakage_nand2_path = 0;
	leakage_nand3_path = 0;
			
	if(ptr_predec_blk_driver->flag_driver_exists){//First check whether a predecoder block driver is needed
		for(i = 0; i < ptr_predec_blk_driver->number_gates_nand2_path; ++i){
			inv = gatearea(INV, 1, ptr_predec_blk_driver->width_nand2_path_p[i],
				ptr_predec_blk_driver->width_nand2_path_n[i], DEFAULTHEIGHTCELL); 
			cumulative_area_nand2_path += inv.area;	
			leakage_nand2_path += cmos_ileakage(ptr_predec_blk_driver->width_nand2_path_n[i],
				ptr_predec_blk_driver->width_nand2_path_p[i], temper);
		}
		cumulative_area_nand2_path *= (ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load +
			ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load +
			ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load);
		leakage_nand2_path *= (ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load +
			ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load +
			ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load);
		for(i = 0; i < ptr_predec_blk_driver->number_gates_nand3_path; ++i){
			inv = gatearea(INV, 1, ptr_predec_blk_driver->width_nand3_path_p[i], 
				ptr_predec_blk_driver->width_nand3_path_n[i], DEFAULTHEIGHTCELL); 
			cumulative_area_nand3_path += inv.area;	
			leakage_nand3_path += cmos_ileakage(ptr_predec_blk_driver->width_nand3_path_n[i],
				ptr_predec_blk_driver->width_nand3_path_p[i], temper);
		}
		cumulative_area_nand3_path *= (ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load +
			ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load);
		leakage_nand3_path *= (ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load +
			ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load);
		ptr_predec_blk_driver->power_nand2_path.readOp.leakage = leakage_nand2_path * vdd_periph_global;
		ptr_predec_blk_driver->power_nand3_path.readOp.leakage = leakage_nand3_path * vdd_periph_global;
		predecoder_block_driver.area = cumulative_area_nand2_path + cumulative_area_nand3_path;
	}		
	return(predecoder_block_driver);
}


area_type area_predecoder(predecoder_block *ptr_predec_blk)
{
	int i;
	double cumulative_area_first_level_nand2_path, cumulative_area_first_level_nand3_path, 
		cumulative_area_first_level, cumulative_area_second_level;
	double leakage_first_level_nand2_path, leakage_first_level_nand3_path, leakage_second_level;
	area_type predec_blk, nand2, nand3, inv;
	ptr_predec_blk->number_first_level_parallel_instances_nand2 = 0;
	ptr_predec_blk->number_first_level_parallel_instances_nand3 = 0;
	ptr_predec_blk->number_second_level_parallel_instances = 0;
	cumulative_area_first_level_nand2_path = 0;
	cumulative_area_first_level_nand3_path = 0;
	cumulative_area_second_level = 0;
	leakage_first_level_nand2_path = 0;
	leakage_first_level_nand3_path = 0;
	leakage_second_level = 0;
	predec_blk.width = 0;
	predec_blk.height = 0;
	predec_blk.area = 0;

	if(ptr_predec_blk->flag_block_exists){//First check whether a predecoder block is needed
		nand2 = gatearea(NAND, 2, ptr_predec_blk->width_first_level_nand2_path_p[0], 
				ptr_predec_blk->width_first_level_nand2_path_n[0], DEFAULTHEIGHTCELL); 
		nand3 = gatearea(NAND, 3, ptr_predec_blk->width_first_level_nand3_path_p[0], 
				ptr_predec_blk->width_first_level_nand3_path_n[0], DEFAULTHEIGHTCELL); 
		cumulative_area_first_level_nand2_path += nand2.area;
		cumulative_area_first_level_nand3_path += nand3.area;
		leakage_first_level_nand2_path += cmos_ileakage(ptr_predec_blk->width_first_level_nand2_path_n[0],
				ptr_predec_blk->width_first_level_nand2_path_p[0], temper) * NAND2_LEAK_STACK_FACTOR;
		leakage_first_level_nand3_path += cmos_ileakage(ptr_predec_blk->width_first_level_nand3_path_n[0],
				ptr_predec_blk->width_first_level_nand3_path_p[0], temper) * NAND3_LEAK_STACK_FACTOR;


		if(ptr_predec_blk->number_input_addr_bits == 1){//2 NAND2 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand2 = 2;
			ptr_predec_blk->number_second_level_parallel_instances = 0;
		}

		if(ptr_predec_blk->number_input_addr_bits == 2){//4 NAND2 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand2 = 4;
			ptr_predec_blk->number_second_level_parallel_instances = 0;
		}

		if(ptr_predec_blk->number_input_addr_bits == 3){//8 NAND3 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand3 = 8;
			ptr_predec_blk->number_second_level_parallel_instances = 0;
		}

		if(ptr_predec_blk->number_input_addr_bits == 4){//4 + 4 NAND2 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand2 = 8;
			ptr_predec_blk->number_second_level_parallel_instances = 16;
		}

		if(ptr_predec_blk->number_input_addr_bits == 5){//4 NAND2 gates, 8 NAND3 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand2 = 4;
			ptr_predec_blk->number_first_level_parallel_instances_nand3 = 8;
		    ptr_predec_blk->number_second_level_parallel_instances = 32;
		}

		if(ptr_predec_blk->number_input_addr_bits == 6){//8 + 8 NAND3 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand3 = 16;
			ptr_predec_blk->number_second_level_parallel_instances = 64;
		}

		if(ptr_predec_blk->number_input_addr_bits == 7){//4 + 4 NAND2 gates, 8 NAND3 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand2 = 8;
			ptr_predec_blk->number_first_level_parallel_instances_nand3 = 8;
			ptr_predec_blk->number_second_level_parallel_instances = 128;
		}

		if(ptr_predec_blk->number_input_addr_bits == 8){//4 NAND2 gates, 8 + 8 NAND3 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand2 = 4;
			ptr_predec_blk->number_first_level_parallel_instances_nand3 = 16;
			ptr_predec_blk->number_second_level_parallel_instances = 256;
		}

		if(ptr_predec_blk->number_input_addr_bits == 9){//8 + 8 + 8 NAND3 gates
			ptr_predec_blk->number_first_level_parallel_instances_nand3 = 24;
			ptr_predec_blk->number_second_level_parallel_instances = 512;
		}

		for(i = 1; i < ptr_predec_blk->number_gates_first_level_nand2_path; ++i){
			inv = gatearea(INV, 1, ptr_predec_blk->width_first_level_nand2_path_p[i], 
				ptr_predec_blk->width_first_level_nand2_path_n[i], DEFAULTHEIGHTCELL); 
			cumulative_area_first_level_nand2_path += inv.area;	
			leakage_first_level_nand2_path += cmos_ileakage(ptr_predec_blk->width_first_level_nand2_path_n[i],
				ptr_predec_blk->width_first_level_nand2_path_p[i], temper);
		}
		cumulative_area_first_level_nand2_path *= ptr_predec_blk->number_first_level_parallel_instances_nand2;
		leakage_first_level_nand2_path *= ptr_predec_blk->number_first_level_parallel_instances_nand2;

		for(i = 1; i < ptr_predec_blk->number_gates_first_level_nand3_path; ++i){
			inv = gatearea(INV, 1, ptr_predec_blk->width_first_level_nand3_path_p[i], 
				ptr_predec_blk->width_first_level_nand3_path_n[i], DEFAULTHEIGHTCELL); 
			cumulative_area_first_level_nand3_path += inv.area;	
			leakage_first_level_nand3_path += cmos_ileakage(ptr_predec_blk->width_first_level_nand3_path_n[i],
				ptr_predec_blk->width_first_level_nand3_path_p[i], temper);
		}
		cumulative_area_first_level_nand3_path *= ptr_predec_blk->number_first_level_parallel_instances_nand3;
		leakage_first_level_nand3_path *= ptr_predec_blk->number_first_level_parallel_instances_nand3;
		cumulative_area_first_level = cumulative_area_first_level_nand2_path +
			cumulative_area_first_level_nand3_path;
		
		cumulative_area_second_level = 0;
		if(ptr_predec_blk->flag_second_level_gate){
			if(ptr_predec_blk->flag_second_level_gate == 2){
				nand2 = gatearea(NAND, 2, ptr_predec_blk->width_second_level_p[0], 
					ptr_predec_blk->width_second_level_n[0], DEFAULTHEIGHTCELL); 
				cumulative_area_second_level += nand2.area;
				leakage_second_level += cmos_ileakage(ptr_predec_blk->width_second_level_n[0],
					ptr_predec_blk->width_second_level_p[0], temper) *
					NAND2_LEAK_STACK_FACTOR;
			}
			else{
				if(ptr_predec_blk->flag_second_level_gate == 3){
					nand3 = gatearea(NAND, 3, ptr_predec_blk->width_second_level_p[0],
						ptr_predec_blk->width_second_level_n[0], DEFAULTHEIGHTCELL); 
					cumulative_area_second_level += nand3.area;
					leakage_second_level += cmos_ileakage(ptr_predec_blk->width_second_level_n[0],
						ptr_predec_blk->width_second_level_p[0], temper) *
						NAND3_LEAK_STACK_FACTOR;
				}
			}
		}
		for(i = 1; i < ptr_predec_blk->number_gates_second_level; ++i){
			inv = gatearea(INV, 1, ptr_predec_blk->width_second_level_p[i], 
				ptr_predec_blk->width_second_level_n[i], DEFAULTHEIGHTCELL); 
			cumulative_area_second_level += inv.area;	
			leakage_second_level += cmos_ileakage(ptr_predec_blk->width_second_level_n[i],
					ptr_predec_blk->width_second_level_p[i], temper);
		}
		cumulative_area_second_level *= ptr_predec_blk->number_second_level_parallel_instances;
		leakage_second_level *= ptr_predec_blk->number_second_level_parallel_instances;

		ptr_predec_blk->power_nand2_path.readOp.leakage = leakage_first_level_nand2_path * vdd_periph_global;
		ptr_predec_blk->power_nand3_path.readOp.leakage = leakage_first_level_nand3_path * vdd_periph_global;
		ptr_predec_blk->power_second_level.readOp.leakage = leakage_second_level * vdd_periph_global; 
		predec_blk.area = cumulative_area_first_level + cumulative_area_second_level;
	}
	return(predec_blk);
}

area_type area_decoder(decoder *ptr_dec)
{
	int i;
	area_type dec, nand2, nand3, inv;
	double cumulative_area, cumulative_power;
	cumulative_area = 0;
	dec.width = 0;
	dec.height = 0;
	dec.area = 0;
	cumulative_power = 0;
	if(ptr_dec->flag_decoder_exists){//First check if this decoder exists
		nand2 = gatearea(NAND, 2, ptr_dec->width_decoder_p[0],ptr_dec->width_decoder_n[0], 
			4 * height_cell); //The height of a row-decoder-driver cell is fixed to be 4 * height_cell;
		nand3 = gatearea(NAND, 3, ptr_dec->width_decoder_p[0], ptr_dec->width_decoder_n[0],
			4 * height_cell); 
		if(ptr_dec->number_input_signals == 2){
			cumulative_area += nand2.area;	
			cumulative_power += cmos_ileakage(ptr_dec->width_decoder_n[0],ptr_dec->width_decoder_p[0], temper) * 
				NAND2_LEAK_STACK_FACTOR;
		}
		else{
			if(ptr_dec->number_input_signals == 3){
				cumulative_area += nand3.area;		
				cumulative_power += cmos_ileakage(ptr_dec->width_decoder_n[0], 
					ptr_dec->width_decoder_p[0], temper) * NAND3_LEAK_STACK_FACTOR;
			}
		}

		for(i = 1; i < ptr_dec->number_gates; ++i){
			inv = gatearea(INV, 1, ptr_dec->width_decoder_p[i], ptr_dec->width_decoder_n[i], 
				4 * height_cell); 
			cumulative_area += inv.area;	
			cumulative_power += cmos_ileakage(ptr_dec->width_decoder_n[i], 
				ptr_dec->width_decoder_p[i], temper);
		}
		ptr_dec->power.readOp.leakage = cumulative_power * vdd_periph_global;
		dec.area = cumulative_area;
	}
	return(dec);
}

area_type area_comparators(int tagbits, int A)
{
	area_type comparator, inv, nand2;
	double cumulative_area;

	comparator.width = 0;
	comparator.height = 0;
	comparator.area = 0;
	tagbits = tagbits / 4;
	cumulative_area = 0;
	//Calculate area of the inverters in the timing chain
	inv = gatearea(INV, 1, Wcompinvp1, Wcompinvn1, DEFAULTHEIGHTCELL); 
	cumulative_area += inv.area;
	inv = gatearea(INV, 1, Wcompinvp2, Wcompinvn2, DEFAULTHEIGHTCELL); 
	cumulative_area += inv.area;
	inv = gatearea(INV, 1, Wcompinvp3, Wcompinvn3, DEFAULTHEIGHTCELL); 
	cumulative_area += inv.area;
	inv = gatearea(INV, 1, Wevalinvp, Wevalinvn, DEFAULTHEIGHTCELL); 
	cumulative_area += inv.area;
	nand2 = gatearea(NAND, 2, 0, Wcompn, DEFAULTHEIGHTCELL); 
	cumulative_area += nand2.area * 2 * tagbits;
	inv = gatearea(INV, 1, Wcompp, 0, DEFAULTHEIGHTCELL); 
	cumulative_area += inv.area;
	cumulative_area *= 4 * A;
	comparator.area = cumulative_area;
	return(comparator);
}


area_type bit_mux_sense_amp_precharge_sa_mux_area(int number_cols_subarray, int deg_bitline_muxing,
												  int Ndsam, double subarray_mem_cell_area_width,
												  int RWP, int ERP, int EWP)
{
	int number_sense_amps_subarray;
	double column_mux_transistor_height,  sense_amp_mux_height, precharge_circuitry_height,
		height_bit_mux_decode_output_wires, height_senseamp_mux_decode_output_wires;
	double cumulative_height;
	area_type sense_amp_area, bit_mux_sense_amp_precharge_sa_mux;
	
	cumulative_height = 0;
	precharge_circuitry_height = 0;
	column_mux_transistor_height = 0;
	height_bit_mux_decode_output_wires = 0;
	sense_amp_mux_height = 0;
	height_senseamp_mux_decode_output_wires = 0;
	bit_mux_sense_amp_precharge_sa_mux.height = 0;
	bit_mux_sense_amp_precharge_sa_mux.width = 0;
	bit_mux_sense_amp_precharge_sa_mux.area = 0;

	precharge_circuitry_height = width_transistor_after_folding(width_pmos_bitline_precharge, width_cell / (2 *(RWP + ERP))) + 
		width_transistor_after_folding(width_pmos_bitline_equalization, width_cell / (RWP + ERP));
	cumulative_height += precharge_circuitry_height;

	if(deg_bitline_muxing > 1){
		column_mux_transistor_height = width_transistor_after_folding(width_nmos_bit_mux, width_cell / (2 *(RWP + ERP + RWP)));
		height_bit_mux_decode_output_wires = deg_bitline_muxing * wire_inside_mat_pitch * (RWP + ERP + RWP);
		cumulative_height += column_mux_transistor_height + height_bit_mux_decode_output_wires;
	}

	number_sense_amps_subarray = number_cols_subarray / deg_bitline_muxing;
	sense_amp_area = area_sense_amplifier(WsenseP, WsenseN, WsenseEn, Wiso, width_cell * deg_bitline_muxing / (RWP + ERP));
	cumulative_height += sense_amp_area.height;
	if(Ndsam > 1){
		sense_amp_mux_height = width_transistor_after_folding(width_nmos_sense_amp_mux, 
			width_cell * Ndsam / (RWP + ERP));
		height_senseamp_mux_decode_output_wires = Ndsam * wire_inside_mat_pitch * (RWP + ERP);
		cumulative_height += sense_amp_mux_height + height_senseamp_mux_decode_output_wires;
	}

	bit_mux_sense_amp_precharge_sa_mux.height = cumulative_height;
	bit_mux_sense_amp_precharge_sa_mux.width = subarray_mem_cell_area_width;
	bit_mux_sense_amp_precharge_sa_mux.area = bit_mux_sense_amp_precharge_sa_mux.height * 
		bit_mux_sense_amp_precharge_sa_mux.width;
	return(bit_mux_sense_amp_precharge_sa_mux);

}

area_type subarray_output_driver_area(int number_cols_subarray, int deg_bitline_muxing,
									  int Ndsam, double subarray_mem_cell_area_width,
									  dataout_htree_node *ptr_htree_node)
{
	int number_output_drivers_subarray, number_gates_subarray_output_driver, j;
	double subarray_output_driver_area;
	area_type subarray_output_driver, subarray_output_driver_nand2, subarray_output_driver_nor2,
		subarray_output_driver_inv;

	subarray_output_driver.width = 0;
	subarray_output_driver.height = 0;
	subarray_output_driver.area = 0;

	number_output_drivers_subarray = number_cols_subarray / (deg_bitline_muxing * Ndsam);
	number_gates_subarray_output_driver = ptr_htree_node->number_gates;
	subarray_output_driver = gatearea(INV, 1, 
		ptr_htree_node->width_p[number_gates_subarray_output_driver-1],
		ptr_htree_node->width_n[number_gates_subarray_output_driver-1], DEFAULTHEIGHTCELL);
    subarray_output_driver_nand2 =  gatearea(NAND, 2,
		ptr_htree_node->width_p[number_gates_subarray_output_driver-2], 
		ptr_htree_node->width_n[number_gates_subarray_output_driver-2], DEFAULTHEIGHTCELL);
	subarray_output_driver_nor2 =  gatearea(NOR, 2, ptr_htree_node->width_nor2_p, 
		ptr_htree_node->width_nor2_n, DEFAULTHEIGHTCELL);
	subarray_output_driver_area = (subarray_output_driver_nand2.area + 
		subarray_output_driver_nor2.area + 	subarray_output_driver.area) * 
		number_output_drivers_subarray;
	for(j = number_gates_subarray_output_driver - 3; j >= 0; --j){
		subarray_output_driver_inv = gatearea(INV, 1, ptr_htree_node->width_p[j],
			ptr_htree_node->width_n[j], DEFAULTHEIGHTCELL);
		subarray_output_driver_area += subarray_output_driver_inv.area * number_output_drivers_subarray;
	}

	subarray_output_driver.height = subarray_output_driver_area / subarray_mem_cell_area_width;
	subarray_output_driver.width = subarray_mem_cell_area_width;
	subarray_output_driver.area = subarray_output_driver_area;
	return(subarray_output_driver);
}


area_type area_htree_node(addr_datain_htree_node *ptr_htree_node)
{
	int i;
	double cumulative_area;
	area_type inv, nand2, htree_node;
	cumulative_area = 0;	
	htree_node.width = 0;
	htree_node.height = 0;
	htree_node.area = 0;
	nand2 = gatearea(NAND, 2, ptr_htree_node->width_p[0], ptr_htree_node->width_n[0], 
		DEFAULTHEIGHTCELL);
	cumulative_area += nand2.area;	
	for(i = 1; i < ptr_htree_node->number_gates; ++i){
		inv = gatearea(INV, 1,  ptr_htree_node->width_p[i], ptr_htree_node->width_n[i],
			DEFAULTHEIGHTCELL);
		cumulative_area += inv.area;
	}
	htree_node.area = cumulative_area;
	return(htree_node);
}


area_type area_dataout_htree_node(dataout_htree_node *ptr_htree_node)
{
	int i;
	double cumulative_area;
	area_type inv, nand2, nor2, htree_node;
	cumulative_area = 0;	
	htree_node.width = 0;
	htree_node.height = 0;
	htree_node.area = 0;
	i = ptr_htree_node->number_gates - 1;
	inv = gatearea(INV, 1,  ptr_htree_node->width_p[i], ptr_htree_node->width_n[i], DEFAULTHEIGHTCELL);
	cumulative_area += inv.area;	
	i = ptr_htree_node->number_gates - 2;
	nand2 = gatearea(NAND, 2, ptr_htree_node->width_p[i], ptr_htree_node->width_n[i], DEFAULTHEIGHTCELL);
	cumulative_area += nand2.area;	
	nor2 = gatearea(NOR, 2, ptr_htree_node->width_nor2_p, ptr_htree_node->width_nor2_n, DEFAULTHEIGHTCELL);
	cumulative_area += nor2.area;	
	for(i = 0; i <= ptr_htree_node->number_gates - 3; ++i){
		inv = gatearea(INV, 1,  ptr_htree_node->width_p[i], ptr_htree_node->width_n[i], 
			DEFAULTHEIGHTCELL);
		cumulative_area += inv.area;
	}
	htree_node.area = cumulative_area;
	return(htree_node);
}

void
calculate_area_drivers_addr_datain_htree_at_mat_interval(addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval)
{
	area_type nand2, inv, nor2;
	int i; 

	//Calculate area of NAND2 driving inverter load driver
	nand2 = gatearea(NAND, 2, ptr_htree_at_mat_interval->nand_driving_inv_width_p[0], 
		ptr_htree_at_mat_interval->nand_driving_inv_width_n[0], DEFAULTHEIGHTCELL);
	ptr_htree_at_mat_interval->area_nand2_driving_inv += nand2.area;
	for(i = 0; i < ptr_htree_at_mat_interval->number_gates_nand_driving_inv; ++i){
		inv = gatearea(INV, 1,  ptr_htree_at_mat_interval->nand_driving_inv_width_p[i], 
			ptr_htree_at_mat_interval->nand_driving_inv_width_n[i], DEFAULTHEIGHTCELL);
		ptr_htree_at_mat_interval->area_nand2_driving_inv += inv.area;
	}

	//Calculate area of inverter driving inverter load driver
	for(i = 0; i < ptr_htree_at_mat_interval->number_gates_inv_driving_inv; ++i){
		inv = gatearea(INV, 1,  ptr_htree_at_mat_interval->inv_driving_inv_width_p[i], 
			ptr_htree_at_mat_interval->inv_driving_inv_width_n[i], DEFAULTHEIGHTCELL);
		ptr_htree_at_mat_interval->area_inv_driving_inv += inv.area;
	}

	//Calculate area of inverter driving NAND2 load driver
	for(i = 0; i < ptr_htree_at_mat_interval->number_gates_inv_driving_inv; ++i){
		inv = gatearea(INV, 1,  ptr_htree_at_mat_interval->inv_driving_inv_width_p[i], 
			ptr_htree_at_mat_interval->inv_driving_inv_width_n[i], DEFAULTHEIGHTCELL);
		ptr_htree_at_mat_interval->area_inv_driving_nand2 += inv.area;
	}

	//Calculate area of NAND2 driving inverter load driver in final segment of H-tree
	nand2 = gatearea(NAND, 2, ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[0], 
		ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0], DEFAULTHEIGHTCELL);
	ptr_htree_at_mat_interval->area_nand2_driving_inv_final_seg += nand2.area;
	for(i = 0; i < ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg; ++i){
		inv = gatearea(INV, 1,  ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[i], 
			ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[i], DEFAULTHEIGHTCELL);
		ptr_htree_at_mat_interval->area_nand2_driving_inv_final_seg += inv.area;
	}

	//Calculate area of tristate driver driving inverter load 
	i = ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 1;
	inv = gatearea(INV, 1,  ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[i], ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[i], DEFAULTHEIGHTCELL);
	ptr_htree_at_mat_interval->area_tristate_driver_driving_inv += inv.area;	
	i = ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 2;
	nand2 = gatearea(NAND, 2, ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[i], ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[i], DEFAULTHEIGHTCELL);
	ptr_htree_at_mat_interval->area_tristate_driver_driving_inv += nand2.area;	
	nor2 = gatearea(NOR, 2, ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p, ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p, DEFAULTHEIGHTCELL);
	ptr_htree_at_mat_interval->area_tristate_driver_driving_inv += nor2.area;	
	for(i = 0; i <= ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 3; ++i){
		inv = gatearea(INV, 1,  ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[i], 
			ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[i], DEFAULTHEIGHTCELL);
		ptr_htree_at_mat_interval->area_tristate_driver_driving_inv += inv.area;
	}
}


area_type area_addr_datain_htree_node_at_mat_interval(addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval, 
													  int number_mats_to_cover_in_this_htree_segment)
{
	double cumulative_area;
	int i; 
	area_type addr_datain_htree_node_at_mat_interval_node;

	cumulative_area = 0;
	if(number_mats_to_cover_in_this_htree_segment == 1){//final segment in H-tree - going over half a 
		//subarray
		cumulative_area += ptr_htree_at_mat_interval->area_nand2_driving_inv_final_seg;
    }
    else{
	    for(i = 0; i < number_mats_to_cover_in_this_htree_segment; ++i){
		    if(i==0){//Driver structure is NAND2 driving inverter load
			    cumulative_area += ptr_htree_at_mat_interval->area_nand2_driving_inv;
		    }
		    else if(i==number_mats_to_cover_in_this_htree_segment - 1){//Driver structure is inverter
			    //driving NAND2 load
			    cumulative_area += ptr_htree_at_mat_interval->area_inv_driving_nand2;
		    }
		    else{//Driver structure is inverter driving inverter load
			    cumulative_area += ptr_htree_at_mat_interval->area_inv_driving_inv;
		    }
	    }
    }
    addr_datain_htree_node_at_mat_interval_node.height = 0;
    addr_datain_htree_node_at_mat_interval_node.width = 0;
    addr_datain_htree_node_at_mat_interval_node.area = cumulative_area;
    return(addr_datain_htree_node_at_mat_interval_node);
}

area_type area_dataout_vertical_htree_at_mat_interval(addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval,
													  int number_mats_to_cover_in_this_htree_segment)
{
	double cumulative_area;
	area_type vertical_dataout_htree_node_at_mat_interval;
	int i;

	cumulative_area = 0;
	for(i = 0; i < number_mats_to_cover_in_this_htree_segment; ++i){
		if(i==0){//Driver structure is  tristate driver driving inverter load
			cumulative_area += ptr_htree_at_mat_interval->area_tristate_driver_driving_inv;
		}
		else{//Driver structure is inverter driving inverter load
			cumulative_area += ptr_htree_at_mat_interval->area_inv_driving_inv;
		}
	}
	vertical_dataout_htree_node_at_mat_interval.height = 0;
	vertical_dataout_htree_node_at_mat_interval.width = 0;
	vertical_dataout_htree_node_at_mat_interval.area = cumulative_area;
	return(vertical_dataout_htree_node_at_mat_interval);
}



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
        dataout_htree_node *subarray_output_htree_node)
{
	 area_type subarray_mem_cell_area;
	 double height_horizontal_non_cell_area_within_mat, width_vertical_non_cell_area_within_mat;
	 double area_mat_center_circuitry, area_rectangle_center_mat, area_row_decoder,
		 width_row_decoder;

	 area_type row_predecoder_block_driver_1,  row_predecoder_block_driver_2, bit_mux_predecoder_block_driver_1, 
		 bit_mux_predecoder_block_driver_2, senseamp_mux_predecoder_block_driver_1, senseamp_mux_predecoder_block_driver_2, 
		 way_sel_driver_1, row_predecoder_block_1, row_predecoder_block_2, bit_mux_predecoder_block_1, bit_mux_predecoder_block_2, 
		 senseamp_mux_predecoder_block_1, senseamp_mux_predecoder_block_2, row_decoder, bit_mux_decoder, senseamp_mux_decoder;

	 int branch_effort_predecoder_block_1_output, branch_effort_predecoder_block_2_output;

	 double   width_row_predecode_output_wires, height_addr_datain_wires_within_mat, height_bit_mux_sense_amp_precharge_sa_mux,
		 height_subarray_output_driver;
	 double area_efficiency_mat;

	area_type data_mat, bit_mux_sense_amp_precharge_sa_mux, subarray_output_driver;

	int RWP, ERP, EWP, NSER;

	 RWP = parameters->rw_ports;
	 ERP = parameters->excl_read_ports;
	 EWP = parameters->excl_write_ports;
	 NSER = parameters->single_ended_read_ports;//Need to think whether NSER has an impact on the area
	 //model anymore. Presently NSER is being used in the SRAM cell width calculation. 

	 subarray_mem_cell_area = subarraymem_area(number_rows_subarray, number_cols_subarray, number_subarrays,
		 parameters->rw_ports, parameters->excl_read_ports, parameters->excl_write_ports,
		 parameters->single_ended_read_ports);

	 row_predecoder_block_driver_1 = area_predecoder_block_driver(row_predec_blk_driver_1, row_predec_blk_1);
	 row_predecoder_block_driver_2 = area_predecoder_block_driver(row_predec_blk_driver_2, row_predec_blk_2);
	 bit_mux_predecoder_block_driver_1 = area_predecoder_block_driver(bit_mux_predec_blk_driver_1, bit_mux_predec_blk_1);
	 bit_mux_predecoder_block_driver_2 = area_predecoder_block_driver(bit_mux_predec_blk_driver_2, bit_mux_predec_blk_2);
	 senseamp_mux_predecoder_block_driver_1 = area_predecoder_block_driver(senseamp_mux_predec_blk_driver_1, senseamp_mux_predec_blk_1);
	 senseamp_mux_predecoder_block_driver_2 = area_predecoder_block_driver(senseamp_mux_predec_blk_driver_2, senseamp_mux_predec_blk_2);
     way_sel_driver_1 =  area_predecoder_block_driver(way_select_driver_1, dummy_way_select_predec_blk_1);

	 row_predecoder_block_1 = area_predecoder(row_predec_blk_1);
	 row_predecoder_block_2 = area_predecoder(row_predec_blk_2);
	 bit_mux_predecoder_block_1 = area_predecoder(bit_mux_predec_blk_1);
	 bit_mux_predecoder_block_2 = area_predecoder(bit_mux_predec_blk_2);
	 senseamp_mux_predecoder_block_1 = area_predecoder(senseamp_mux_predec_blk_1);
	 senseamp_mux_predecoder_block_2 = area_predecoder(senseamp_mux_predec_blk_2);
	
	 row_decoder = area_decoder(row_dec);
	 bit_mux_decoder = area_decoder(bit_mux_dec);
	 senseamp_mux_decoder = area_decoder(senseamp_mux_dec);
	 area_row_decoder = row_decoder.area * number_rows_subarray * (RWP + ERP + EWP);
	 width_row_decoder = area_row_decoder / subarray_mem_cell_area.height;
	 

	 bit_mux_sense_amp_precharge_sa_mux = bit_mux_sense_amp_precharge_sa_mux_area(
		 number_cols_subarray, deg_bitline_muxing, Ndsam, subarray_mem_cell_area.width,
		 RWP, ERP, EWP);
	 height_bit_mux_sense_amp_precharge_sa_mux = bit_mux_sense_amp_precharge_sa_mux.height;
	 subarray_output_driver = subarray_output_driver_area(number_cols_subarray, deg_bitline_muxing,
		 Ndsam, subarray_mem_cell_area.width, subarray_output_htree_node);
	 height_subarray_output_driver = subarray_output_driver.height * (RWP + ERP);
	
			
	 branch_effort_predecoder_block_1_output = (int) (pow(2, row_predec_blk_2->number_input_addr_bits));
	 branch_effort_predecoder_block_2_output = (int) (pow(2, row_predec_blk_1->number_input_addr_bits));
	 width_row_predecode_output_wires = (branch_effort_predecoder_block_1_output + branch_effort_predecoder_block_2_output) * 
		 wire_inside_mat_pitch * (RWP + ERP + EWP);
	 
	 height_horizontal_non_cell_area_within_mat = 2 * (height_bit_mux_sense_amp_precharge_sa_mux + 
		 height_subarray_output_driver);
	 width_vertical_non_cell_area_within_mat = MAX(width_row_predecode_output_wires, 
		 2 * width_row_decoder);
	 if(!VERTICAL_HTREE_WIRES_OVER_THE_ARRAY){
         /* low swing wires will take more area */
         if (parameters->wire_inter_mats == Low_swing) {
             height_addr_datain_wires_within_mat = (number_addr_bits_mat + number_way_select_signals_mat +
                     number_datain_bits_mat / 2 + number_dataout_bits_mat / 2) * 4*wire_inside_mat_pitch * (
                         RWP + ERP + EWP);
             height_horizontal_non_cell_area_within_mat = 2 * height_bit_mux_sense_amp_precharge_sa_mux +
             height_addr_datain_wires_within_mat;
//                 MAX(height_addr_datain_wires_within_mat, 2 *  height_subarray_output_driver);
         }
         else {
             height_addr_datain_wires_within_mat = (number_addr_bits_mat + number_way_select_signals_mat +
                     number_datain_bits_mat / 2 + number_dataout_bits_mat / 2) * wire_inside_mat_pitch * (
                         RWP + ERP + EWP);
             height_horizontal_non_cell_area_within_mat = 2 * height_bit_mux_sense_amp_precharge_sa_mux +
//                 MAX(height_addr_datain_wires_within_mat, 2 *  height_subarray_output_driver);
             height_addr_datain_wires_within_mat;
         }
	 }

	 area_rectangle_center_mat = height_horizontal_non_cell_area_within_mat *
		 width_vertical_non_cell_area_within_mat;
     area_mat_center_circuitry = (row_predecoder_block_driver_1.area  + 
			 bit_mux_predecoder_block_driver_1.area + senseamp_mux_predecoder_block_driver_1.area  +
			 way_sel_driver_1.area + row_predecoder_block_driver_2.area +
			 bit_mux_predecoder_block_driver_2.area + senseamp_mux_predecoder_block_driver_2.area + 
			 row_predecoder_block_1.area + bit_mux_predecoder_block_1.area + 
			 senseamp_mux_predecoder_block_1.area + row_predecoder_block_2.area +
			 bit_mux_predecoder_block_2.area +  senseamp_mux_predecoder_block_2.area + 
			 bit_mux_decoder.area +  senseamp_mux_decoder.area) * (RWP + ERP + EWP);

	
	data_mat.height = 2 * subarray_mem_cell_area.height + height_horizontal_non_cell_area_within_mat;
	data_mat.width = 2 * subarray_mem_cell_area.width + width_vertical_non_cell_area_within_mat;
	data_mat.area = data_mat.height * data_mat.width + area_mat_center_circuitry;
	area_efficiency_mat = subarray_mem_cell_area.area * 4 * 100 / data_mat.area;
	data_mat.width = data_mat.area / data_mat.height;
	return(data_mat);
}


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
            powerDef *tot_power_senseamp_mux_decoders)
{
	area_type horizontal_addr_din_htree_node, vertical_addr_din_htree_node, dout_htree_node,
		bank, tag_comparators;
	double area_horizontal_addr_datain_routing, width_horizontal_addr_datain_routing,
		height_horizontal_addr_datain_routing, area_vertical_addr_datain_routing, area_dataout_routing, 
		width_vertical_addr_datain_routing,
		width_vertical_dataout_routing, height_vertical_addr_datain_routing, height_vertical_dataout_routing,
		area_horizontal_interconnect_routing, width_horizontal_interconnect_routing,
		height_horizontal_interconnect_routing, area_tag_comparators,
		area_active_horizontal_addr_datain_dataout_routing, height_active_area_horizontal_addr_datain_routing,
		height_wire_area_horizontal_addr_datain_dataout_routing, width_wire_area_vertical_addr_datain_routing, 
		width_wire_area_vertical_dataout_routing;
	double asp_ratio_temp;
	int htree_seg_multiplier, i;
	int number_mats, number_mats_to_cover_in_this_htree_segment;
	int RWP, ERP, EWP, NSER;

	 RWP = parameters->rw_ports;
	 ERP = parameters->excl_read_ports;
	 EWP = parameters->excl_write_ports;
	 NSER = parameters->single_ended_read_ports;

    area_horizontal_addr_datain_routing = 0;
	area_vertical_addr_datain_routing = 0;
	area_dataout_routing = 0;
	area_active_horizontal_addr_datain_dataout_routing = 0;
	width_horizontal_addr_datain_routing = 0;
	width_vertical_addr_datain_routing = 0;
	width_vertical_dataout_routing = 0;
	height_active_area_horizontal_addr_datain_routing = 0;
	height_wire_area_horizontal_addr_datain_dataout_routing = 0;
	height_horizontal_addr_datain_routing = 0;
    height_vertical_addr_datain_routing = 0;
	height_vertical_dataout_routing = 0;
	area_horizontal_interconnect_routing = 0;
	width_horizontal_interconnect_routing = 0;
	height_horizontal_interconnect_routing = 0;
	area_tag_comparators = 0;
	width_wire_area_vertical_addr_datain_routing = 0;
	width_wire_area_vertical_dataout_routing = 0;

	number_mats = number_mats_horizontal_direction * number_mats_vertical_direction;
	
	if(HTREE_NODES_AT_MAT_INTERVALS){
		calculate_area_drivers_addr_datain_htree_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval);
		calculate_area_drivers_addr_datain_htree_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval);
	}
	htree_seg_multiplier = 1;
	for(i = 0; i < number_horizontal_htree_nodes; ++i){
		horizontal_addr_din_htree_node = area_htree_node(&ptr_horizontal_addr_din_htree_node[i]);
		if(HTREE_NODES_AT_MAT_INTERVALS){
			number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_horizontal_htree_nodes - i) + 0.5)) / 2, 1);
			horizontal_addr_din_htree_node = area_addr_datain_htree_node_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment);
		}
		area_active_horizontal_addr_datain_dataout_routing += 
			(number_addr_bits_mat + number_way_select_signals_mat) * horizontal_addr_din_htree_node.area *
			htree_seg_multiplier * (RWP + ERP + EWP) + number_datain_bits_mat * number_mats_horizontal_direction * 
			horizontal_addr_din_htree_node.area * (RWP + EWP) + number_dataout_bits_mat * 
			number_mats_horizontal_direction * horizontal_addr_din_htree_node.area * (RWP + ERP);
		height_wire_area_horizontal_addr_datain_dataout_routing += 
			(number_addr_bits_mat + number_way_select_signals_mat) * wire_outside_mat_pitch *
			(RWP + ERP + EWP) + (number_datain_bits_mat * number_mats_horizontal_direction *
			wire_outside_mat_pitch / htree_seg_multiplier) * (RWP + EWP) + (number_dataout_bits_mat *
			number_mats_horizontal_direction * wire_outside_mat_pitch / htree_seg_multiplier) * 
			(RWP + ERP);
		htree_seg_multiplier *= 2;
	}
	
	
	//If tag array, add area occupied by comparators
	if(is_tag){
		tag_comparators = area_comparators(tagbits, parameters->tag_associativity);
		area_tag_comparators = tag_comparators.area * (RWP + ERP);
	}
	area_active_horizontal_addr_datain_dataout_routing += area_tag_comparators;

	width_horizontal_addr_datain_routing = number_mats_horizontal_direction * width_mat;
	height_active_area_horizontal_addr_datain_routing = area_active_horizontal_addr_datain_dataout_routing / 
			width_horizontal_addr_datain_routing;
	height_horizontal_addr_datain_routing = MAX(height_active_area_horizontal_addr_datain_routing,
		height_wire_area_horizontal_addr_datain_dataout_routing);

	
	htree_seg_multiplier = 1; 
	for(i = 0; i < number_vertical_htree_nodes; ++i){
		vertical_addr_din_htree_node = area_htree_node(&ptr_addr_din_htree_node[i]);
		if(HTREE_NODES_AT_MAT_INTERVALS){
			number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_vertical_htree_nodes - i) + 0.5)) / 2, 1);
			vertical_addr_din_htree_node = 
				area_addr_datain_htree_node_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment);
		}
		area_vertical_addr_datain_routing += ((number_addr_bits_mat + number_datain_bits_mat + 
			number_way_select_signals_mat) * number_mats_horizontal_direction *	
			htree_seg_multiplier * vertical_addr_din_htree_node.area) * (RWP + EWP);
         if (parameters->wire_inter_mats == Low_swing) {
             width_wire_area_vertical_addr_datain_routing = ((number_addr_bits_mat + number_way_select_signals_mat + 
                         number_datain_bits_mat) * number_mats_horizontal_direction * 4*wire_outside_mat_pitch) * 
                 (RWP + EWP);
         }
         else {
             width_wire_area_vertical_addr_datain_routing = ((number_addr_bits_mat + number_way_select_signals_mat + 
                         number_datain_bits_mat) * number_mats_horizontal_direction * wire_outside_mat_pitch) * 
                 (RWP + EWP);
         }
             htree_seg_multiplier *= 2;
	}
	if(number_mats_vertical_direction > 1){
		area_vertical_addr_datain_routing *= 2; 
	}
	
	width_vertical_addr_datain_routing = number_mats_horizontal_direction * width_mat;
	height_vertical_addr_datain_routing = area_vertical_addr_datain_routing / width_vertical_addr_datain_routing;

	if(!VERTICAL_HTREE_WIRES_OVER_THE_ARRAY){//vertical H-tree wires don't go over the array
		height_vertical_addr_datain_routing = 0;
	}

	htree_seg_multiplier = MAX(((int)(pow(2.0, number_vertical_htree_nodes - 1) + 0.5)) / 2, 1);
	for(i = 0; i < number_vertical_htree_nodes - 1; ++i){
		dout_htree_node = area_dataout_htree_node(&ptr_dout_htree_node[i]);
		if(HTREE_NODES_AT_MAT_INTERVALS){
			number_mats_to_cover_in_this_htree_segment = MAX((int)(pow(2.0, i) + 0.5), 1);
			dout_htree_node = 
				area_dataout_vertical_htree_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment);
		}
		area_dataout_routing += (number_dataout_bits_mat * number_mats_horizontal_direction * 
			htree_seg_multiplier * dout_htree_node.area) * (RWP + ERP);
         if (parameters->wire_inter_mats == Low_swing) {
             width_wire_area_vertical_dataout_routing += (number_dataout_bits_mat * 
                     number_mats_horizontal_direction * 4*wire_outside_mat_pitch) * (RWP + ERP);
         }
         else {
             width_wire_area_vertical_dataout_routing += (number_dataout_bits_mat * 
                     number_mats_horizontal_direction * 4*wire_outside_mat_pitch) * (RWP + ERP);
         }
		htree_seg_multiplier /= 2;
	}

	if(number_mats_vertical_direction > 1){
		area_dataout_routing *= 2; 
	}

	width_vertical_dataout_routing = number_mats_horizontal_direction * width_mat;
	height_vertical_dataout_routing = area_dataout_routing / width_vertical_dataout_routing;

	if(!VERTICAL_HTREE_WIRES_OVER_THE_ARRAY){//vertical H-tree wires don't go over the array
		height_vertical_dataout_routing = 0;
	}

	bank.height = number_mats_vertical_direction * height_mat + 
		height_horizontal_addr_datain_routing + height_vertical_addr_datain_routing + height_vertical_dataout_routing;
	bank.width = number_mats_horizontal_direction * width_mat;

	if(!VERTICAL_HTREE_WIRES_OVER_THE_ARRAY){
		bank.width = number_mats_horizontal_direction * width_mat + 
			width_wire_area_vertical_addr_datain_routing + width_wire_area_vertical_dataout_routing;	
	}

	result->totalarea = bank.height * bank.width;
	asp_ratio_temp = bank.height / bank.width;
	result->aspect_ratio_total =(asp_ratio_temp  > 1.0) ? (asp_ratio_temp) : 1.0 / asp_ratio_temp;

	tot_power_row_predecode_block_drivers->readOp.leakage = 
		(row_predec_blk_driver_1->power_nand2_path.readOp.leakage + 
		row_predec_blk_driver_1->power_nand3_path.readOp.leakage + 
		row_predec_blk_driver_2->power_nand2_path.readOp.leakage + 
		row_predec_blk_driver_2->power_nand3_path.readOp.leakage) * number_mats;

	tot_power_bit_mux_predecode_block_drivers->readOp.leakage =
		(bit_mux_predec_blk_driver_1->power_nand2_path.readOp.leakage + 
		bit_mux_predec_blk_driver_1->power_nand3_path.readOp.leakage +
		bit_mux_predec_blk_driver_2->power_nand2_path.readOp.leakage + 
		bit_mux_predec_blk_driver_2->power_nand3_path.readOp.leakage) * number_mats;

	tot_power_senseamp_mux_predecode_block_drivers->readOp.leakage = 
		(senseamp_mux_predec_blk_driver_1->power_nand2_path.readOp.leakage + 
		senseamp_mux_predec_blk_driver_1->power_nand3_path.readOp.leakage + 
		senseamp_mux_predec_blk_driver_2->power_nand2_path.readOp.leakage + 
		senseamp_mux_predec_blk_driver_2->power_nand3_path.readOp.leakage) * number_mats;

	tot_power->readOp.leakage += tot_power_row_predecode_block_drivers->readOp.leakage +  
		tot_power_bit_mux_predecode_block_drivers->readOp.leakage + 
		tot_power_senseamp_mux_predecode_block_drivers->readOp.leakage;

	tot_power_row_predecode_blocks->readOp.leakage = (row_predec_blk_1->power_nand2_path.readOp.leakage + 
		row_predec_blk_1->power_nand3_path.readOp.leakage + row_predec_blk_1->power_second_level.readOp.leakage +
		row_predec_blk_2->power_nand2_path.readOp.leakage+ row_predec_blk_2->power_nand3_path.readOp.leakage + 
		row_predec_blk_2->power_second_level.readOp.leakage) * number_mats;

	tot_power_bit_mux_predecode_blocks->readOp.leakage = (bit_mux_predec_blk_1->power_nand2_path.readOp.leakage + 
		bit_mux_predec_blk_1->power_nand3_path.readOp.leakage + bit_mux_predec_blk_1->power_second_level.readOp.leakage +
		bit_mux_predec_blk_2->power_nand2_path.readOp.leakage + bit_mux_predec_blk_2->power_nand3_path.readOp.leakage + 
		bit_mux_predec_blk_2->power_second_level.readOp.leakage) * number_mats;

	tot_power_senseamp_mux_predecode_blocks->readOp.leakage = (senseamp_mux_predec_blk_1->power_nand2_path.readOp.leakage + 
		senseamp_mux_predec_blk_1->power_nand3_path.readOp.leakage + senseamp_mux_predec_blk_1->power_second_level.readOp.leakage +
		senseamp_mux_predec_blk_2->power_nand2_path.readOp.leakage + senseamp_mux_predec_blk_2->power_nand3_path.readOp.leakage + 
		senseamp_mux_predec_blk_2->power_second_level.readOp.leakage) * number_mats;

	tot_power->readOp.leakage += tot_power_row_predecode_blocks->readOp.leakage + 
		tot_power_bit_mux_predecode_blocks->readOp.leakage + tot_power_senseamp_mux_predecode_blocks->readOp.leakage;

	tot_power_row_decoders->readOp.leakage = row_dec->power.readOp.leakage * number_rows_subarray * 4 * 
		number_mats;
	tot_power_bit_mux_decoders->readOp.leakage = bit_mux_dec->power.readOp.leakage * number_mats;
	tot_power_senseamp_mux_decoders->readOp.leakage = senseamp_mux_dec->power.readOp.leakage * number_mats;

	tot_power->readOp.leakage += tot_power_row_decoders->readOp.leakage +
		tot_power_bit_mux_decoders->readOp.leakage + tot_power_senseamp_mux_decoders->readOp.leakage;
	return(bank);
}


area_type area_all_banks(int uca_banks, double bank_height, double bank_width, 
						 int number_bits_routed_to_bank, double *length_htree_route_to_bank,
						 int number_mats_vertical_direction)
{
	double pitch_of_routed_tracks_one_bank, all_banks_height, all_banks_width,
		curr_pitch_routed_tracks, curr_htree_vertical_seg_length, 
		curr_htree_horizontal_seg_length;
	area_type all_banks;
	int i, number_segments_htree, first_seg_vertical, uca_banks_horizontal_direction,
		uca_banks_vertical_direction, curr_horizontal_routed_tracks_instances, 
		curr_vertical_routed_tracks_instances, curr_seg_vertical,
		number_vertical_segments_htree, number_horizontal_segments_htree;
	
	*length_htree_route_to_bank = 0;
	all_banks_height = 0;
	all_banks_width = 0;
	pitch_of_routed_tracks_one_bank = (number_bits_routed_to_bank) * wire_outside_mat_pitch; 
	if(uca_banks > 1){
		number_segments_htree = (int)(logbasetwo((double)(uca_banks)));
		if(number_segments_htree%2 == 0){
			first_seg_vertical = 0;
		}
		else{
			first_seg_vertical = 1;
		}
		if(logbasefour(uca_banks) == floor(logbasefour((double) (uca_banks)))){
			uca_banks_horizontal_direction = (int)(sqrt(uca_banks));
		}
		else{
			uca_banks_horizontal_direction = (int)(sqrt(uca_banks / 2)) * 2;
		}
		uca_banks_vertical_direction = uca_banks / uca_banks_horizontal_direction;
		all_banks_height = uca_banks_vertical_direction * bank_height;
		all_banks_width = uca_banks_horizontal_direction * bank_width;
		if(first_seg_vertical){
			number_vertical_segments_htree = (int) (number_segments_htree / 2) + 1;
		}
		else{
			number_vertical_segments_htree = number_segments_htree / 2;
		}
		number_horizontal_segments_htree = number_segments_htree - number_vertical_segments_htree;

		if(number_mats_vertical_direction == 1){
			number_vertical_segments_htree = number_vertical_segments_htree - 1;
			number_segments_htree = number_segments_htree - 1;
		}
	
		curr_pitch_routed_tracks = uca_banks * pitch_of_routed_tracks_one_bank;
		curr_horizontal_routed_tracks_instances = 1;
		curr_vertical_routed_tracks_instances = 1;
		curr_seg_vertical = first_seg_vertical;
		curr_htree_vertical_seg_length = uca_banks_vertical_direction * bank_height / 2;
		curr_htree_horizontal_seg_length = uca_banks_horizontal_direction * bank_width / 2;
		for(i = 0; i < number_segments_htree; ++i){
			if(curr_seg_vertical){
				all_banks_width += curr_pitch_routed_tracks * curr_vertical_routed_tracks_instances;
				*length_htree_route_to_bank += curr_htree_vertical_seg_length;
				curr_seg_vertical = 0;
				curr_htree_vertical_seg_length = curr_htree_vertical_seg_length / 2;
				curr_vertical_routed_tracks_instances = curr_vertical_routed_tracks_instances * 2;
			}
			else{
				all_banks_height += curr_pitch_routed_tracks * curr_horizontal_routed_tracks_instances;
				curr_seg_vertical = 1;
				*length_htree_route_to_bank += curr_htree_horizontal_seg_length;
				curr_htree_horizontal_seg_length = curr_htree_horizontal_seg_length / 2;
				curr_horizontal_routed_tracks_instances = curr_horizontal_routed_tracks_instances * 2;
			}
			curr_pitch_routed_tracks = curr_pitch_routed_tracks / 2;
		}
	}
	else{//uca_banks == 1;
		all_banks_width = bank_width;
		all_banks_height = bank_height;
	}
	all_banks.height = all_banks_height;
	all_banks.width = all_banks_width;
	return(all_banks);
}
