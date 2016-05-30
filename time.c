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

#include "time.h"

#ifdef __linux__
pthread_mutex_t com_mut_var = PTHREAD_MUTEX_INITIALIZER;
#endif
#define USE_GLOBAL

//#define DNUCA

int MIN_BANKSIZE=65536; /* *2B */
int NUCA = 0;
float cumm_per[1024];
int start_index = 0;
double FREQUENCY = 5; //GHz

int CONTR_2_BANK_LAT = 0; /* in cycles */
//double contention[4][7];
int cont_stats[2 /*l2 or l3*/][5/* cores */][ROUTER_TYPES][7 /*banks*/][8 /* cycle time */];
int core_in = 0;
//router_stats_t router[ROUTER_TYPES];
router_stats_t router_s[ROUTER_TYPES];

void initialize_addr_datain_htree(
        addr_datain_htree_node *ptr_htree_node, 
        int number_htree_nodes,
        double htree_seg_length,  
        double *delay, powerDef *power, 
        double *length_wire_htree_node)
{
    int i, j, multiplier;
	double input_gate_cap, F, c_load;

    if(number_htree_nodes == 0){
        number_htree_nodes = 1;
    }
    *delay = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;

    for(i = 0; i < MAX_NUMBER_HTREE_NODES; ++i){
		ptr_htree_node[i].power.readOp.dynamic = 0;
        ptr_htree_node[i].power.readOp.leakage = 0;
        ptr_htree_node[i].power.writeOp.dynamic = 0;
        ptr_htree_node[i].power.writeOp.leakage = 0;
        for(j = 0; j < MAX_NUMBER_GATES_STAGE; ++j){
            ptr_htree_node[i].width_n[j] = 0;
            ptr_htree_node[i].width_p[j] = 0;
         }
    }

    multiplier = 1;
	for(i = number_htree_nodes - 1; i >= 0; --i){
		ptr_htree_node[i].min_number_gates = 2;
		length_wire_htree_node[i] = multiplier * htree_seg_length;
        ptr_htree_node[i].r_wire_load = wire_outside_mat_r_per_micron * length_wire_htree_node[i];
        ptr_htree_node[i].c_wire_load = wire_outside_mat_c_per_micron * length_wire_htree_node[i];
		ptr_htree_node[i].f_wire = ptr_htree_node[i].r_wire_load * ptr_htree_node[i].c_wire_load / (2 * kinv);
		multiplier *= 2;
	}

	//Calculate width of input gate of each H-tree segment
	ptr_htree_node[0].width_n[0] = 2 * minimum_width_nmos;
	ptr_htree_node[0].width_p[0] = ptr_htree_node[0].width_n[0];
	for(i = number_htree_nodes - 1; i > 0; --i){
		if(i == number_htree_nodes - 1){
			ptr_htree_node[i].c_gate_load = gatecap(minimum_width_nmos + minimum_width_pmos, 0);
		}
        else{
			ptr_htree_node[i].c_gate_load = 2 * gatecap(ptr_htree_node[i+1].width_n[0] +
				ptr_htree_node[i+1].width_p[0], 0);
		}
		c_load = ptr_htree_node[i].c_gate_load + ptr_htree_node[i].c_wire_load;
		F = (ptr_htree_node[i-1].r_wire_load * gnand2 * c_load / (2 * kinv * fopt)) * 
			(1 + sqrt(1 + 2 * fopt / ptr_htree_node[i-1].f_wire));
		if(F <= pow(fopt, 2)){
			//Include a two-stage buffer
			F = pow(fopt, 2);
		}
		input_gate_cap = gnand2 * c_load / F;
		ptr_htree_node[i].width_n[0] = (input_gate_cap / 2) / gatecap(1, 0);
		if(ptr_htree_node[i].width_n[0] < 2 * minimum_width_nmos){
			ptr_htree_node[i].width_n[0] = 2 * minimum_width_nmos;
		}
		if(ptr_htree_node[i].width_n[0] > MAX_NMOS_WIDTH){
			//Two-stage buffer with 2nd stage composed of MAX_NMOS_WIDTH. First stage of buffer
			//is MAX_NMOS_WIDTH divided by fopt. 
			c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
			input_gate_cap = gnand2 * c_load / F;
			ptr_htree_node[i].width_n[0] = (input_gate_cap / 2) / gatecap(1, 0);
		}
		ptr_htree_node[i].width_p[0] = ptr_htree_node[i].width_n[0];
	}
	if(number_htree_nodes > 1){
		ptr_htree_node[0].c_gate_load = gatecap(ptr_htree_node[1].width_n[0] + ptr_htree_node[1].width_p[0], 0);
	}
	else{
		ptr_htree_node[0].c_gate_load = gatecap(minimum_width_nmos + minimum_width_pmos, 0);
	}
}


void initialize_addr_datain_htree_with_nodes_at_mat_interval(
        addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval, 
        double *delay, powerDef *power, double mat_dimension)
{
    int j;
	//double c_input;

    *delay = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;

   for(j = 0; j < MAX_NUMBER_GATES_STAGE; ++j){
	   ptr_htree_at_mat_interval->nand_driving_inv_width_n[j] = 0;
	   ptr_htree_at_mat_interval->nand_driving_inv_width_p[j] = 0;
	   ptr_htree_at_mat_interval->inv_driving_inv_width_n[j] = 0;
	   ptr_htree_at_mat_interval->inv_driving_inv_width_p[j] = 0;
	   ptr_htree_at_mat_interval->inv_driving_nand_width_n[j] = 0;
	   ptr_htree_at_mat_interval->inv_driving_nand_width_p[j] = 0;
	   ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j] = 0;
	   ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j] = 0;
	   ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] = 0;
	   ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j] = 0;
   }

   ptr_htree_at_mat_interval->delay_nand_driving_inv = 0;
   ptr_htree_at_mat_interval->delay_inv_driving_inv = 0;
   ptr_htree_at_mat_interval->delay_inv_driving_nand = 0;
   ptr_htree_at_mat_interval->delay_tristate_driver_driving_inv = 0;
   ptr_htree_at_mat_interval->delay_nand_driving_inv_final_seg = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.leakage = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv.writeOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv.writeOp.leakage = 0;
   ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.leakage = 0;
   ptr_htree_at_mat_interval->power_inv_driving_inv.writeOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_inv_driving_inv.writeOp.leakage = 0;
   ptr_htree_at_mat_interval->power_inv_driving_nand.readOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_inv_driving_nand.readOp.leakage = 0;
   ptr_htree_at_mat_interval->power_inv_driving_nand.writeOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_inv_driving_nand.writeOp.leakage = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.writeOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.writeOp.leakage = 0;
   ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage = 0;
   ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.writeOp.dynamic = 0;
   ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.writeOp.leakage = 0;
   ptr_htree_at_mat_interval->area_nand2_driving_inv = 0;
   ptr_htree_at_mat_interval->area_inv_driving_inv = 0;
   ptr_htree_at_mat_interval->area_inv_driving_nand2 = 0;
   ptr_htree_at_mat_interval->area_nand2_driving_inv_final_seg = 0;
   ptr_htree_at_mat_interval->area_tristate_driver_driving_inv = 0;
   ptr_htree_at_mat_interval->min_number_gates = 2;
   ptr_htree_at_mat_interval->min_number_gates_tristate_driver = 4;
   ptr_htree_at_mat_interval->r_wire_load = wire_outside_mat_r_per_micron * mat_dimension;
   ptr_htree_at_mat_interval->c_wire_load = wire_outside_mat_c_per_micron * mat_dimension;
   ptr_htree_at_mat_interval->r_wire_load_final_seg = wire_outside_mat_r_per_micron * mat_dimension / 2;
   ptr_htree_at_mat_interval->c_wire_load_final_seg = wire_outside_mat_c_per_micron * mat_dimension / 2;

   //if(INPUT_WIRE_TO_INPUT_GATE_CAP_RATIO == 0){
	   ptr_htree_at_mat_interval->nand_driving_inv_width_n[0] = 2 * minimum_width_nmos;
	   ptr_htree_at_mat_interval->nand_driving_inv_width_p[0] = 2 * minimum_width_nmos;
	   ptr_htree_at_mat_interval->inv_driving_inv_width_n[0] = minimum_width_nmos;
	   ptr_htree_at_mat_interval->inv_driving_inv_width_p[0] = minimum_width_pmos;
	   ptr_htree_at_mat_interval->inv_driving_nand_width_n[0] = minimum_width_nmos;
	   ptr_htree_at_mat_interval->inv_driving_nand_width_p[0] = minimum_width_pmos;
	   ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0] = 2 * minimum_width_nmos;
	   ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[0] = 2 * minimum_width_nmos;
	   ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[0] = minimum_width_nmos;
	   ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[0] = minimum_width_pmos;
   /*}
   else{
	  ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv_final_seg =
		gatecap(minimum_width_nmos + minimum_width_pmos, 0);
	  c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv_final_seg +
		   ptr_htree_at_mat_interval->c_wire_load_final_seg;
	  F = (ptr_htree_at_mat_interval->r_wire_load * gnand2 * c_load / (2 * kinv * fopt)) * 
		  (1 + sqrt(1 + 2 * fopt / ptr_htree_at_mat_interval->f_wire));
	  input_gate_cap = gnand2 * c_load / F;
	  ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0] = 
		  (input_gate_cap / 2) / gatecap(1, 0);
	  if(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0] < 2 * minimum_width_nmos){
		  ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0] = 2 * minimum_width_nmos;
	  }
	   ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[0] = 
		   ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0];

	   c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_nand +
		   ptr_htree_at_mat_interval->c_wire_load;
	   F = (ptr_htree_at_mat_interval->c_wire_load  / ptr_htree_at_mat_interval->c_wire_load) * 
				(1 + sqrt(1 + 2 * fopt / ptr_htree_at_mat_interval->f_wire)) / 
				(sqrt(1 + 2 * fopt / ptr_htree_at_mat_interval->f_wire) - 1);
	   input_gate_cap = c_load / F;
	   ptr_htree_at_mat_interval->inv_driving_nand_width_n[0] = (1.0 / 3.0) * input_gate_cap / gatecap(1, 0);
	   if(ptr_htree_at_mat_interval->inv_driving_nand_width_n[0] < 2 * minimum_width_nmos){
		   ptr_htree_at_mat_interval->inv_driving_nand_width_n[0] = 2 * minimum_width_nmos;
	   }
	   ptr_htree_at_mat_interval->inv_driving_nand_width_p[0] =
		   2 * ptr_htree_at_mat_interval->inv_driving_nand_width_n[0];

	   ptr_htree_at_mat_interval->inv_driving_inv_width_n[0] = 
		   ptr_htree_at_mat_interval->inv_driving_nand_width_n[0];
	   ptr_htree_at_mat_interval->inv_driving_inv_width_p[0] =
		   ptr_htree_at_mat_interval->inv_driving_nand_width_p[0];
	   ptr_htree_at_mat_interval->nand_driving_inv_width_n[0] = 
		   (input_gate_cap / 2) / gatecap(1, 0);
	   if(ptr_htree_at_mat_interval->nand_driving_inv_width_n[0] < 2 * minimum_width_nmos){
		   ptr_htree_at_mat_interval->nand_driving_inv_width_n[0] = 2 * minimum_width_nmos;
	   }
	   ptr_htree_at_mat_interval->nand_driving_inv_width_p[0] = 
		   ptr_htree_at_mat_interval->nand_driving_inv_width_n[0];
	   ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[0] =
		   ptr_htree_at_mat_interval->inv_driving_inv_width_n[0];
	   ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[0] =
		   ptr_htree_at_mat_interval->inv_driving_inv_width_p[0];
   }*/
 
   ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv = 
	   gatecap(ptr_htree_at_mat_interval->inv_driving_inv_width_n[0] +
	   ptr_htree_at_mat_interval->inv_driving_inv_width_p[0], 0);
   ptr_htree_at_mat_interval->c_gate_load_inv_driving_inv = 
	    gatecap(ptr_htree_at_mat_interval->inv_driving_inv_width_n[0] +
	   ptr_htree_at_mat_interval->inv_driving_inv_width_p[0], 0);
   ptr_htree_at_mat_interval->c_gate_load_inv_driving_nand = 2 * 
	   gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[0] +
	   ptr_htree_at_mat_interval->nand_driving_inv_width_p[0], 0);
   ptr_htree_at_mat_interval->c_gate_load_tristate_driver_driving_inv = 
	    gatecap(ptr_htree_at_mat_interval->inv_driving_inv_width_n[0] +
		ptr_htree_at_mat_interval->inv_driving_inv_width_p[0], 0);
	ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv_final_seg =
		gatecap(minimum_width_nmos + minimum_width_pmos, 0);
}



void compute_widths_addr_datain_htree(addr_datain_htree_node *ptr_htree_node, int number_htree_nodes)
{
    int i, j;
    double c_load, c_input, F, f;
    if(number_htree_nodes == 0){
        number_htree_nodes = 1;
    }
    for(i = 0; i < number_htree_nodes; ++i){
        c_load = ptr_htree_node[i].c_gate_load + ptr_htree_node[i].c_wire_load;
        F =  gnand2 * c_load / gatecap(ptr_htree_node[i].width_n[0] + ptr_htree_node[i].width_p[0], 0);
        ptr_htree_node[i].number_gates = (int) (log(F) / log(fopt) + 0.5);
        f = pow(F, 1.0 / ptr_htree_node[i].number_gates);
        if(ptr_htree_node[i].number_gates%2 != 0){
            --ptr_htree_node[i].number_gates;
        }
		if(ptr_htree_node[i].number_gates < ptr_htree_node[i].min_number_gates){
            ptr_htree_node[i].number_gates = ptr_htree_node[i].min_number_gates;
        }
        f = pow(F, 1.0 / ptr_htree_node[i].number_gates);

		j = ptr_htree_node[i].number_gates - 1;
		c_input = c_load / f;
		ptr_htree_node[i].width_n[j] = (1.0 / 3.0) * c_input / gatecap(1, 0);
		if(ptr_htree_node[i].width_n[j] < minimum_width_nmos){
			ptr_htree_node[i].width_n[j] = minimum_width_nmos;
		}
		ptr_htree_node[i].width_p[j]  = 2 * ptr_htree_node[i].width_n[j];

		if(ptr_htree_node[i].width_n[j] > MAX_NMOS_WIDTH){
			c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
			F =  gnand2 * c_load / gatecap(ptr_htree_node[i].width_n[0] +
				ptr_htree_node[i].width_p[0], 0);
			ptr_htree_node[i].number_gates = (int) (log(F) / log(fopt) + 0.5) + 1;
			if(ptr_htree_node[i].number_gates%2 != 0){
				--ptr_htree_node[i].number_gates;
			}
			if(ptr_htree_node[i].number_gates < ptr_htree_node[i].min_number_gates){
				ptr_htree_node[i].number_gates = ptr_htree_node[i].min_number_gates;
			}
			f = pow(F, 1.0 / ptr_htree_node[i].number_gates);
			j = ptr_htree_node[i].number_gates - 1;
			ptr_htree_node[i].width_n[j] = MAX_NMOS_WIDTH;
			ptr_htree_node[i].width_p[j]  = 2 * ptr_htree_node[i].width_n[j];
		}
		
		for(j = ptr_htree_node[i].number_gates - 2; j > 0; --j){
            ptr_htree_node[i].width_n[j] = ptr_htree_node[i].width_n[j+1] / f;
            if(ptr_htree_node[i].width_n[j] < minimum_width_nmos){
                ptr_htree_node[i].width_n[j] = minimum_width_nmos;
            }
            ptr_htree_node[i].width_p[j]  = 2 * ptr_htree_node[i].width_n[j];
        }
    }
}


void compute_widths_addr_datain_htree_with_nodes_at_mat_interval(
	addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval)
{
    int j;
    double c_load, F, f, H, B, G, c_input;

	//Compute widths for NAND2 driving inverter load driver
	c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv  + 
		ptr_htree_at_mat_interval->c_wire_load;
	F =  gnand2 * c_load / gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[0] + 
		ptr_htree_at_mat_interval->nand_driving_inv_width_p[0], 0);
	ptr_htree_at_mat_interval->number_gates_nand_driving_inv = (int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_nand_driving_inv);
	if(ptr_htree_at_mat_interval->number_gates_nand_driving_inv%2 != 0){
		++ptr_htree_at_mat_interval->number_gates_nand_driving_inv;
	}
	if(ptr_htree_at_mat_interval->number_gates_nand_driving_inv < ptr_htree_at_mat_interval->min_number_gates){
		ptr_htree_at_mat_interval->number_gates_nand_driving_inv = ptr_htree_at_mat_interval->min_number_gates;
	}
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_nand_driving_inv);

	j = ptr_htree_at_mat_interval->number_gates_nand_driving_inv - 1;
	c_input = c_load / f;
	ptr_htree_at_mat_interval->nand_driving_inv_width_n[j] = (1.0 / 3.0) * (c_input / gatecap(1, 0));
	ptr_htree_at_mat_interval->nand_driving_inv_width_p[j] = 2 *
		ptr_htree_at_mat_interval->nand_driving_inv_width_n[j];
	for(j = ptr_htree_at_mat_interval->number_gates_nand_driving_inv - 2; j > 0; --j){
		ptr_htree_at_mat_interval->nand_driving_inv_width_n[j] =
			ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] / f;
		if(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j] < minimum_width_nmos){
			ptr_htree_at_mat_interval->nand_driving_inv_width_n[j] = minimum_width_nmos;
		}
		ptr_htree_at_mat_interval->nand_driving_inv_width_p[j]  = 2 * 
			ptr_htree_at_mat_interval->nand_driving_inv_width_n[j];
	}

	//Compute widths for inverter driving inverter load driver
	c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_inv  + 
		ptr_htree_at_mat_interval->c_wire_load;
	F =  c_load / gatecap(ptr_htree_at_mat_interval->inv_driving_inv_width_n[0] + 
		ptr_htree_at_mat_interval->inv_driving_inv_width_p[0], 0);
	ptr_htree_at_mat_interval->number_gates_inv_driving_inv = (int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_inv_driving_inv);
	if(ptr_htree_at_mat_interval->number_gates_inv_driving_inv%2 != 0){
		++ptr_htree_at_mat_interval->number_gates_inv_driving_inv;
	}
	if(ptr_htree_at_mat_interval->number_gates_inv_driving_inv < ptr_htree_at_mat_interval->min_number_gates){
		ptr_htree_at_mat_interval->number_gates_inv_driving_inv = ptr_htree_at_mat_interval->min_number_gates;
	}
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_inv_driving_inv);

	j = ptr_htree_at_mat_interval->number_gates_inv_driving_inv - 1;
	c_input = c_load / f;
	ptr_htree_at_mat_interval->inv_driving_inv_width_n[j] = (1.0 / 3.0) * (c_input / gatecap(1, 0));
	ptr_htree_at_mat_interval->inv_driving_inv_width_p[j] = 2 *
		ptr_htree_at_mat_interval->inv_driving_inv_width_n[j];
	for(j = ptr_htree_at_mat_interval->number_gates_inv_driving_inv - 2; j > 0; --j){
		c_input = gatecap(ptr_htree_at_mat_interval->inv_driving_inv_width_n[j], 0);
		ptr_htree_at_mat_interval->inv_driving_inv_width_n[j] =
			ptr_htree_at_mat_interval->inv_driving_inv_width_n[j+1] / f;
		if(ptr_htree_at_mat_interval->inv_driving_inv_width_n[j] < minimum_width_nmos){
			ptr_htree_at_mat_interval->inv_driving_inv_width_n[j] = minimum_width_nmos;
		}
		ptr_htree_at_mat_interval->inv_driving_inv_width_p[j]  = 2 * 
			ptr_htree_at_mat_interval->inv_driving_inv_width_n[j];
	}
	
	//Compute widths for inverter driving NAND2 load driver
	c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_nand  + 
		ptr_htree_at_mat_interval->c_wire_load;
	F =  c_load / gatecap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[0] + 
		ptr_htree_at_mat_interval->inv_driving_nand_width_p[0], 0);
	ptr_htree_at_mat_interval->number_gates_inv_driving_nand = (int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_inv_driving_nand);
	if(ptr_htree_at_mat_interval->number_gates_inv_driving_nand%2 != 0){
		++ptr_htree_at_mat_interval->number_gates_inv_driving_nand;
	}
	if(ptr_htree_at_mat_interval->number_gates_inv_driving_nand < ptr_htree_at_mat_interval->min_number_gates){
		ptr_htree_at_mat_interval->number_gates_inv_driving_nand = ptr_htree_at_mat_interval->min_number_gates;
	}
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_inv_driving_nand);
	j = ptr_htree_at_mat_interval->number_gates_inv_driving_nand - 1;
	c_input = c_load / f;
	ptr_htree_at_mat_interval->inv_driving_nand_width_n[j] = (1.0 / 3.0) * (c_input / gatecap(1, 0));
	ptr_htree_at_mat_interval->inv_driving_nand_width_p[j] = 2 *
		ptr_htree_at_mat_interval->inv_driving_nand_width_n[j];
	for(j = ptr_htree_at_mat_interval->number_gates_inv_driving_nand - 2; j > 0; --j){
		c_input = gatecap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], 0);
		ptr_htree_at_mat_interval->inv_driving_nand_width_n[j] =
			ptr_htree_at_mat_interval->inv_driving_nand_width_n[j+1] / f;
		if(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j] < minimum_width_nmos){
			ptr_htree_at_mat_interval->inv_driving_nand_width_n[j] = minimum_width_nmos;
		}
		ptr_htree_at_mat_interval->inv_driving_nand_width_p[j]  = 2 * 
			ptr_htree_at_mat_interval->inv_driving_nand_width_n[j];
	}

	//Compute widths for NAND2 driving inverter load driver of final segment in H-tree
	c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv_final_seg  + 
		ptr_htree_at_mat_interval->c_wire_load_final_seg;
	F =  gnand2 * c_load / gatecap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[0] + 
		ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[0] , 0);
	ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg = (int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg);
	if(ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg%2 != 0){
		++ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg;
	}
	if(ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg < ptr_htree_at_mat_interval->min_number_gates){
		ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg = ptr_htree_at_mat_interval->min_number_gates;
	}
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg);
	j = ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg - 1;
	c_input = c_load / f;
	ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j] = (1.0 / 3.0) * (c_input / gatecap(1, 0));
	ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j] = 2 *
		ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j];
	for(j = ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg - 2; j > 0; --j){
		c_input = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], 0);
		ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j] =
			ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j+1] / f;
		if(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j] < minimum_width_nmos){
			ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j] = minimum_width_nmos;
		}
		ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j]  = 2 * 
			ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j];
	}

	//Compute widths for tristate driver driving inverter load (in dataout H-tree)
	B = (gnand2 * gpmos + gnor2 * gnmos) / (gnand2 * gpmos);
    G = gnand2 * gpmos;
	c_load = ptr_htree_at_mat_interval->c_gate_load_tristate_driver_driving_inv + 
		ptr_htree_at_mat_interval->c_wire_load;
	H = c_load / gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[0] + 
		ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[0] , 0);
	F = G * B * H;
	ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv = 
		(int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv);
	if(ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv%2 != 0){
		++ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv;
	}
	if(ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv < ptr_htree_at_mat_interval->min_number_gates_tristate_driver){
            ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv = 
				ptr_htree_at_mat_interval->min_number_gates_tristate_driver;
	}
	f = pow(F, 1.0 / ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv);
	j = ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 1;
	ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j] = gpmos * c_load /(f * gatecap(1, 0));
	if(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j] < minimum_width_pmos){
		ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j] = minimum_width_pmos;
	}
	ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] = ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j] / 2;

	 j = ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 2;
	 c_input = gnand2 * 
		 gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0) / f;
	 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] = 
		 (c_input / gatecap(1, 0)) / 2;
	 if(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] < 2 * minimum_width_nmos){
		 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] = 
			 2 * minimum_width_nmos;
	 }
	 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j] =
		 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j];
	 c_input = gnor2 * 
		 gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1], 0) / f;
	 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n  = 
		 (1.0 / 5.0) * c_input / gatecap(1, 0);
	 if(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n  < minimum_width_nmos){
		 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n  = minimum_width_nmos;
	 }
	 ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p =
		 4 * ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n;


	 for(j = ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 3; j >= 1; --j){
			c_input = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1] +
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n + 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p, 0) / f;
            ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] = 
				(1.0 / 3.0) * (c_input / gatecap(1, 0));
            if(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] < minimum_width_nmos){
                ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j] = minimum_width_nmos;
            }
            ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j]  = 
				2 * ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j];
        }
}


void initialize_dataout_htree(
        dataout_htree_node *ptr_htree_node, 
        int number_htree_nodes,
        double *delay, powerDef *power,
        double *length_wire_htree_node)
{
    int i, j;
	double input_gate_cap, F, c_load, nand2_gate_cap;

    *delay = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;
    for(i = 0; i < MAX_NUMBER_HTREE_NODES; ++i){
		 ptr_htree_node[i].power.readOp.dynamic = 0;
		 ptr_htree_node[i].power.readOp.leakage = 0;
		 ptr_htree_node[i].power.writeOp.dynamic = 0;
		 ptr_htree_node[i].power.writeOp.leakage = 0;
        for(j = 0; j < MAX_NUMBER_GATES_STAGE; ++j){
            ptr_htree_node[i].width_n[j] = 0;
            ptr_htree_node[i].width_p[j] = 0;
        }
    }

	for(i = number_htree_nodes - 2; i >= 0; --i){
        ptr_htree_node[i].min_number_gates = 2;
        ptr_htree_node[i].r_wire_load = wire_outside_mat_r_per_micron * length_wire_htree_node[i];
        ptr_htree_node[i].c_wire_load = wire_outside_mat_c_per_micron * length_wire_htree_node[i];
		ptr_htree_node[i].f_wire = ptr_htree_node[i].r_wire_load * ptr_htree_node[i].c_wire_load / (2 * kinv);
    }

	for(i = 0; i < number_htree_nodes - 1; ++i){
		if(i == number_htree_nodes - 2){
			ptr_htree_node[i].width_n[0] = 2 * minimum_width_nmos;
			ptr_htree_node[i].c_gate_load =  gatecap(ptr_htree_node[i-1].width_n[0] +
					ptr_htree_node[i-1].width_p[0] + ptr_htree_node[i-1].width_nor2_n + ptr_htree_node[i-1].width_nor2_p , 0);
		}
		else{
			if(i == 0){
				ptr_htree_node[i].c_gate_load = gatecap(minimum_width_nmos + minimum_width_pmos, 0);
			}
			else{
				ptr_htree_node[i].c_gate_load =  gatecap(ptr_htree_node[i-1].width_n[0] +
					ptr_htree_node[i-1].width_p[0] + ptr_htree_node[i-1].width_nor2_n + ptr_htree_node[i-1].width_nor2_p , 0);
			}
			c_load = ptr_htree_node[i].c_gate_load + ptr_htree_node[i].c_wire_load;
			F = (ptr_htree_node[i+1].r_wire_load * gnand2 * gpmos * c_load / (2 * kinv * fopt)) *
				(1 + sqrt( 1 + 2 * fopt / ptr_htree_node[i+1].f_wire));
			if(F <= pow(fopt, 2)){
				//Assume a two-stage tristate driver
				F = pow(fopt, 2);
			}
			input_gate_cap = gnand2 * gpmos * c_load / F;
            ptr_htree_node[i].width_n[0] = (input_gate_cap / 2) / gatecap(1, 0);
            if(ptr_htree_node[i].width_n[0] < 2 * minimum_width_nmos){
				ptr_htree_node[i].width_n[0] = 2 * minimum_width_nmos;
			}
			if(ptr_htree_node[i].width_n[0] > MAX_NMOS_WIDTH){
				//Two-stage buffer with 2nd stage composed of MAX_PMOS_WIDTH. First stage NAND2 width
				//calculated as follows
				c_load = gatecap(2 * MAX_NMOS_WIDTH, 0);
				input_gate_cap = gnand2 * c_load / fopt;
				ptr_htree_node[i].width_n[0] = (input_gate_cap / 2) / gatecap(1, 0);
			}
		}
		ptr_htree_node[i].width_p[0] = ptr_htree_node[i].width_n[0];
		nand2_gate_cap = gatecap(ptr_htree_node[i].width_n[0] + ptr_htree_node[i].width_p[0], 0);
		input_gate_cap = gnor2 * gnmos * nand2_gate_cap / (gnand2 * gpmos);
		ptr_htree_node[i].width_nor2_n = (1.0 / 5.0) * input_gate_cap / gatecap(1, 0);
        if(ptr_htree_node[i].width_nor2_n < minimum_width_nmos){
            ptr_htree_node[i].width_nor2_n = minimum_width_nmos;
        }
        ptr_htree_node[i].width_nor2_p = 4 * ptr_htree_node[i].width_nor2_n;
	}
}


void compute_widths_dataout_htree(dataout_htree_node *ptr_htree_node, int number_htree_nodes)
{

    int i, j;
    double G, c_load, H, f, B, F, c_input;
    B = 1;
    G = gnand2 * gpmos;
    for(i = number_htree_nodes - 2; i >= 0; --i){
        c_load = ptr_htree_node[i].c_gate_load + ptr_htree_node[i].c_wire_load;
        H = c_load / gatecap(ptr_htree_node[i].width_n[0] + ptr_htree_node[i].width_p[0] , 0);
        F = G * B * H;
        ptr_htree_node[i].number_gates = (int) (log(F) / log(fopt) + 0.5);
        f = pow(F, 1.0 / ptr_htree_node[i].number_gates);
        if(ptr_htree_node[i].number_gates%2 != 0){
            --ptr_htree_node[i].number_gates;
        }
        if(ptr_htree_node[i].number_gates < ptr_htree_node[i].min_number_gates){
            ptr_htree_node[i].number_gates = ptr_htree_node[i].min_number_gates;
        }
        f = pow(F, 1.0 / ptr_htree_node[i].number_gates);
        j = ptr_htree_node[i].number_gates - 1;
        ptr_htree_node[i].width_p[j] = gpmos * c_load /(f * gatecap(1, 0));
        if(ptr_htree_node[i].width_p[j] < minimum_width_pmos){
            ptr_htree_node[i].width_p[j] = minimum_width_pmos;
        }
		ptr_htree_node[i].width_n[j] = ptr_htree_node[i].width_p[j] / 2;

		if(ptr_htree_node[i].width_p[j] > MAX_PMOS_WIDTH){
			c_load = gatecap(MAX_PMOS_WIDTH + MAX_PMOS_WIDTH / 2, 0);
			F =  gnand2 * c_load / gatecap(ptr_htree_node[i].width_n[0] +
				ptr_htree_node[i].width_p[0] , 0);
			ptr_htree_node[i].number_gates = (int) (log(F) / log(fopt) + 0.5) + 1;
			if(ptr_htree_node[i].number_gates%2 != 0){
				--ptr_htree_node[i].number_gates;
			}
			if(ptr_htree_node[i].number_gates < ptr_htree_node[i].min_number_gates){
				ptr_htree_node[i].number_gates = ptr_htree_node[i].min_number_gates;
			}
			f = pow(F, 1.0 / ptr_htree_node[i].number_gates);
			j = ptr_htree_node[i].number_gates - 1;
			ptr_htree_node[i].width_p[j] = MAX_PMOS_WIDTH;
			ptr_htree_node[i].width_n[j] = ptr_htree_node[i].width_p[j] / 2;
		}
 
        for(j = ptr_htree_node[i].number_gates - 2; j >= 1; --j){
			if(j == ptr_htree_node[i].number_gates - 2){
				c_input = gatecap(ptr_htree_node[i].width_p[j], 0);
			}
			else{
				c_input = gatecap(ptr_htree_node[i].width_n[j+1] + 
					ptr_htree_node[i].width_p[j+1], 0) / f;
			}
            ptr_htree_node[i].width_n[j] = (1.0 / 3.0) * (c_input / gatecap(1, 0));
            if(ptr_htree_node[i].width_n[j] < minimum_width_nmos){
                ptr_htree_node[i].width_n[j] = minimum_width_nmos;
            }
            ptr_htree_node[i].width_p[j]  = 2 * ptr_htree_node[i].width_n[j];
        }
    }
}

void
compute_widths_subarray_output_driver(dataout_htree_node *ptr_htree_node, int number_htree_nodes)
{
    int j;
    double G, H, F, f, c_load, c_input, B, nand2_gate_cap;
   
    ptr_htree_node->min_number_gates = 2;
    ptr_htree_node->c_gate_load = gatecap(2 * minimum_width_nmos + 2 * minimum_width_nmos +
		minimum_width_nmos + 4 * minimum_width_nmos, 0);//Minimum-size NAND2 + Minimum-size NOR2
    ptr_htree_node->r_wire_load = wire_outside_mat_r_per_micron * height_subarray;
    ptr_htree_node->c_wire_load = wire_outside_mat_c_per_micron * height_subarray;
    ptr_htree_node->power.readOp.dynamic = 0;
    ptr_htree_node->power.readOp.leakage = 0;
    ptr_htree_node->power.writeOp.dynamic = 0;
    ptr_htree_node->power.writeOp.leakage = 0;
    c_load = ptr_htree_node->c_gate_load + ptr_htree_node->c_wire_load;

	ptr_htree_node->width_n[0] = 2 * minimum_width_nmos;
    ptr_htree_node->width_p[0] = ptr_htree_node->width_n[0];
	nand2_gate_cap = gatecap(ptr_htree_node->width_n[0] + ptr_htree_node->width_p[0], 0);
	c_input = gnor2 * gnmos * nand2_gate_cap / (gnand2 * gpmos);
	ptr_htree_node->width_nor2_n = (1.0 / 5.0) * c_input / gatecap(1, 0);
	if(ptr_htree_node->width_nor2_n < minimum_width_nmos){
		ptr_htree_node->width_nor2_n = minimum_width_nmos;
	}
	ptr_htree_node->width_nor2_p = 4 * ptr_htree_node->width_nor2_n;


	B = 1;
    G = gnand2 * gpmos;
	c_load = ptr_htree_node->c_gate_load + ptr_htree_node->c_wire_load;
	H = c_load / gatecap(ptr_htree_node->width_n[0] + ptr_htree_node->width_p[0] , 0);
	F = G * B * H;
	ptr_htree_node->number_gates = (int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_htree_node->number_gates);
	if(ptr_htree_node->number_gates%2 != 0){
		--ptr_htree_node->number_gates;
	}
	if(ptr_htree_node->number_gates < ptr_htree_node->min_number_gates){
		ptr_htree_node->number_gates = ptr_htree_node->min_number_gates;
	}
	f = pow(F, 1.0 / ptr_htree_node->number_gates);
	j = ptr_htree_node->number_gates - 1;
	ptr_htree_node->width_p[j] = gpmos * c_load /(f * gatecap(1, 0));
	if(ptr_htree_node->width_p[j] < minimum_width_pmos){
		ptr_htree_node->width_p[j] = minimum_width_pmos;
	}
	ptr_htree_node->width_n[j] = ptr_htree_node->width_p[j] / 2;

	if(ptr_htree_node->width_p[j] > MAX_PMOS_WIDTH){
		c_load = gatecap(MAX_PMOS_WIDTH + MAX_PMOS_WIDTH / 2, 0);
		F =  gnand2 * c_load / gatecap(ptr_htree_node->width_n[0] +	ptr_htree_node->width_p[0], 0);
		ptr_htree_node->number_gates = (int) (log(F) / log(fopt) + 0.5) + 1;
		if(ptr_htree_node->number_gates%2 != 0){
			--ptr_htree_node->number_gates;
		}
		if(ptr_htree_node->number_gates < ptr_htree_node->min_number_gates){
			ptr_htree_node->number_gates = ptr_htree_node->min_number_gates;
		}
		f = pow(F, 1.0 / ptr_htree_node->number_gates);
		j = ptr_htree_node->number_gates - 1;
		ptr_htree_node->width_p[j] = MAX_PMOS_WIDTH;
		ptr_htree_node->width_n[j] = ptr_htree_node->width_p[j] / 2;
	}

	for(j = ptr_htree_node->number_gates - 2; j >= 1; --j){
		if(j == ptr_htree_node->number_gates - 2){
			c_input = gatecap(ptr_htree_node->width_p[j], 0);
		}
		else{
			c_input = gatecap(ptr_htree_node->width_n[j+1] + ptr_htree_node->width_p[j+1], 0) / f;
		}
		ptr_htree_node->width_n[j] = (1.0 / 3.0) * (c_input / gatecap(1, 0));
		if(ptr_htree_node->width_n[j] < minimum_width_nmos){
			ptr_htree_node->width_n[j] = minimum_width_nmos;
		}
		ptr_htree_node->width_p[j]  = 2 * ptr_htree_node->width_n[j];
	}
}


void delay_addr_datain_htree(addr_datain_htree_node *ptr_htree_node, int number_htree_nodes, 
        int hhtree_nodes, double inrisetime, double *outrisetime, double *delay, 
        input_params_t *parameters, area_type *data_array)
{
    int i, j;
    double rd, c_load, c_intrinsic, tf, this_delay;
//    pda_res_t ls_pda;
    wire_stats_t wst;
    this_delay = 0;
    *delay = 0;
    double len1, size1, tc, nand_d = 0;
    
 
    if (parameters->wire_inter_mats == Low_swing) {
        rd = (data_array->height/2 + data_array->width);
        wst.wt = Low_swing;
        wst.wire_length = rd*1e-6; //m
        wst.nsense = number_htree_nodes/2;
        calc_wire_stats2 (Low_swing, &wst);
        *delay += wst.wire_pda.delay;
//        power->readOp.dynamic = wst.wire_pda.power.dynamic;
//        power->readOp.leakage = wst.wire_pda.power.leakage;
        *outrisetime = signal_rise_time();
        return;
//        ls_pda.area_stats.height = rd*1e-6;//(m)
//        calc_wire_stats (Low_swing, &ls_pda, number_htree_nodes);
//        *delay += ls_pda.delay;
//        power->readOp.dynamic = ls_pda.power.dynamic;
//        power->readOp.leakage = ls_pda.power.leakage;
//        *outrisetime = signal_rise_time();
//        return;
    }
#ifdef USE_GLOBAL 
	else {// if (parameters->wire_inter_mats == Global) {
		rd = (data_array->height/2 + data_array->width);
        wst.wt = parameters->wire_inter_mats;
        wst.wire_length = rd*1e-6; //m
        calc_wire_stats2 (parameters->wire_inter_mats, &wst);
        *delay += wst.wire_pda.delay;
        
        /* delay due to nand gates in the htree */
        i = hhtree_nodes-1;
        j = 4;
        if (i != 0) {
            while (i!=1) {
                len1 = data_array->width/j;
                if (wst.repeater_spacing > len1) {
                    size1 = wst.repeater_size*len1/wst.repeater_spacing;
                }
                else {
                    size1 = wst.repeater_size;
                }
                tc = 2*transreson((size1/4)*minimum_width_nmos, NCH, 1) *
                    draincap((size1/4)*minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL)*2;
                nand_d = horowitz (signal_rise_time(), tc, 0.5, 0.5, RISE);
                i--;
                j*=2;
            }
            *delay += nand_d;
        }
            i = number_htree_nodes-1;
            j=2;
        if (i != 0) {
            while (i!=1) {
                len1 =  data_array->height/j;
                if (wst.repeater_spacing > len1) {
                    size1 = wst.repeater_size*len1/wst.repeater_spacing;
                }
                else {
                    size1 = wst.repeater_size;
                }
                tc = 2*transreson((size1/4)*minimum_width_nmos, NCH, 1) *
                    draincap((size1/4)*minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL)*2;
                nand_d = horowitz (signal_rise_time(), tc, 0.5, 0.5, RISE);
                i--;
                j*=2;
            }
            *delay += nand_d;
        }

//        power->readOp.dynamic = wst.wire_pda.power.dynamic;
//        power->readOp.leakage = wst.wire_pda.power.leakage;
//        *outrisetime = signal_rise_time();
//        ls_pda.area_stats.height = rd*1e-6;//(m)
//        calc_wire_stats (Global, &ls_pda, 1);
//        *delay += (ls_pda.delay);
//        power->readOp.dynamic = ls_pda.power.dynamic;
//        power->readOp.leakage = ls_pda.power.leakage;
        *outrisetime = signal_rise_time();
        return;
	}
#endif
	for(i = 0; i < number_htree_nodes ; ++i){
        for(j = 0; j < ptr_htree_node[i].number_gates; ++j){
			if(j == 0){//NAND2 gate
                rd = transreson(ptr_htree_node[i].width_n[j], NCH, 2);
                c_load = gatecap(ptr_htree_node[i].width_n[j+1] + ptr_htree_node[i].width_p[j+1], 0.0);
                c_intrinsic = 2 * draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_htree_node[i].width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
            }
            else if(j == ptr_htree_node[i].number_gates - 1){
				rd = transreson(ptr_htree_node[i].width_n[j], NCH, 1);
				c_load = ptr_htree_node[i].c_gate_load + ptr_htree_node[i].c_wire_load;
				c_intrinsic = draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
					draincap(ptr_htree_node[i].width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
				tf = rd * (c_intrinsic + c_load) + 
					ptr_htree_node[i].r_wire_load * (ptr_htree_node[i].c_wire_load / 2 +
					ptr_htree_node[i].c_gate_load);
			}
			else{//inverter
				rd = transreson(ptr_htree_node[i].width_n[j], NCH, 1);
				c_load = gatecap(ptr_htree_node[i].width_n[j+1] + ptr_htree_node[i].width_p[j+1], 0.0);
				c_intrinsic = draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
					draincap(ptr_htree_node[i].width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
				tf = rd * ( c_intrinsic + c_load);
			}
			this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
            *delay += this_delay;
            inrisetime = this_delay/(1.0 - 0.5);
		}
	}
    *outrisetime = inrisetime;
}

void power_addr_datain_htree_node(addr_datain_htree_node *ptr_htree_node, powerDef *power, 
    input_params_t *parameters, area_type *data_array, int no_nodes, int hhtree_nodes)
{
    int i,j;
    double c_load, c_intrinsic, rd;
    double len1, size1, tc, nand_d = 0;
    wire_stats_t wst;
    
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;
    
    if (parameters->wire_inter_mats == Low_swing) {
        rd = (data_array->height/2 + data_array->width);
        wst.wt = Low_swing;
        wst.wire_length = rd*1e-6; //m
        wst.nsense = no_nodes/2;
        calc_wire_stats2 (Low_swing, &wst);
        power->readOp.dynamic = wst.wire_pda.power.dynamic;
        power->readOp.leakage = wst.wire_pda.power.leakage;
        return;
    }
#ifdef USE_GLOBAL 
	else {
		rd = (data_array->height/2 + data_array->width*3/2);
        wst.wt = parameters->wire_inter_mats;
        wst.wire_length = rd*1e-6; //m
        calc_wire_stats2 (parameters->wire_inter_mats, &wst);
        power->readOp.dynamic = wst.wire_pda.power.dynamic;
        power->readOp.leakage = wst.wire_pda.power.leakage;
        /* power due to the nand gates in the htree */
        i = hhtree_nodes-1;
        j = 4;
        if (i != 0) {
            while (i!=1) {
                len1 = data_array->width/j;
                if (wst.repeater_spacing > len1) {
                    size1 = wst.repeater_size*len1/wst.repeater_spacing;
                }
                else {
                    size1 = wst.repeater_size;
                }
                power->readOp.dynamic += 0.5 * draincap((size1/4)*minimum_width_nmos, 
                        NCH, 1, 1, DEFAULTHEIGHTCELL)*2*vdd_periph_global;
                i--;
                j*=2;
            }
        }
            i = no_nodes-1;
            j=2;
        if (i != 0) {
            while (i!=1) {
                len1 =  data_array->height/j;
                if (wst.repeater_spacing > len1) {
                    size1 = wst.repeater_size*len1/wst.repeater_spacing;
                }
                else {
                    size1 = wst.repeater_size;
                }
                power->readOp.dynamic += 0.5 * draincap((size1/4)*minimum_width_nmos, 
                        NCH, 1, 1, DEFAULTHEIGHTCELL)*2*vdd_periph_global;
                i--;
                j*=2;
            }
        }

        return;
	}
#endif
    for(j = 0; j < ptr_htree_node->number_gates; ++j){
		if(j == 0){//NAND2 gate
			c_load = gatecap(ptr_htree_node->width_n[j+1] + ptr_htree_node->width_p[j+1], 0.0);
			c_intrinsic = 2 * draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_n[j], 
                        ptr_htree_node->width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else if(j == ptr_htree_node->number_gates - 1){
			c_load = ptr_htree_node->c_gate_load + ptr_htree_node->c_wire_load;
			c_intrinsic = draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_n[j], 
				ptr_htree_node->width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else{//inverter
			c_load = gatecap(ptr_htree_node->width_n[j+1] + ptr_htree_node->width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_n[j], 
				ptr_htree_node->width_p[j], temper) * 0.5 * vdd_periph_global;
		}
	}
}


void delay_addr_datain_htree_with_nodes_at_mat_interval(addr_datain_htree_at_mat_interval 
													   *ptr_htree_at_mat_interval, 
													   int number_htree_nodes, 
													   double inrisetime, double *outrisetime, 
													   double *delay)
{
    int j, i, number_mats_to_cover_in_this_htree_segment;
    double rd, c_load, c_intrinsic, tf, this_delay;

    this_delay = 0;
    *delay = 0;
    
    //Find delay of NAND2 driving inverter
	inrisetime = 0;
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_nand_driving_inv; ++j){
		if(j == 0){//NAND2 gate
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 2);
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = 2 * draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load);
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_nand_driving_inv - 1){
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1);
			c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv + 
				ptr_htree_at_mat_interval->c_wire_load;
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load) + ptr_htree_at_mat_interval->r_wire_load * 
				ptr_htree_at_mat_interval->c_wire_load / 2;
		}
		else{//inverter
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1);
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL)
				+ draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * ( c_intrinsic + c_load);
		}
		this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
		ptr_htree_at_mat_interval->delay_nand_driving_inv += this_delay;
		inrisetime = this_delay/(1.0 - 0.5);
	}

	//Find delay of inverter driving inverter
	inrisetime = 0;
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_inv_driving_inv; ++j){
		if(j == ptr_htree_at_mat_interval->number_gates_inv_driving_inv - 1){
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1);
			c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_inv + 
				ptr_htree_at_mat_interval->c_wire_load;
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load) + ptr_htree_at_mat_interval->r_wire_load * 
				ptr_htree_at_mat_interval->c_wire_load / 2;
		}
		else{//inverter
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1);
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * ( c_intrinsic + c_load);
		}
		this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
		ptr_htree_at_mat_interval->delay_inv_driving_inv += this_delay;
		inrisetime = this_delay/(1.0 - 0.5);
	}

	//Find delay of inverter driving NAND2
	inrisetime = 0;
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_inv_driving_nand; ++j){
		if(j == ptr_htree_at_mat_interval->number_gates_inv_driving_nand - 1){//NAND2 gate
			rd = transreson(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], NCH, 1);
			c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_nand + 
				ptr_htree_at_mat_interval->c_wire_load;
			c_intrinsic = 2 * draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load) + ptr_htree_at_mat_interval->r_wire_load * 
				ptr_htree_at_mat_interval->c_wire_load / 2;
		}
		else{//inverter
			rd = transreson(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], NCH, 1);
			c_load = gatecap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j+1] + 
				ptr_htree_at_mat_interval->inv_driving_nand_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * ( c_intrinsic + c_load);
		}
	}
	this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
	ptr_htree_at_mat_interval->delay_inv_driving_nand += this_delay;
	inrisetime = this_delay/(1.0 - 0.5);

	//Find delay of NAND2 driving inverter for final segment of H-tree (over half a mat height/width)
	inrisetime = 0;
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg; ++j){
		if(j == 0){//NAND2 gate
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 2);
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j+1], 0.0);
			c_intrinsic = 2 * draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load);
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg - 1){
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 1);
			c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv + 
				ptr_htree_at_mat_interval->c_wire_load_final_seg;
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load) + ptr_htree_at_mat_interval->r_wire_load_final_seg * 
				ptr_htree_at_mat_interval->c_wire_load_final_seg / 2;
		}
		else{//inverter
			rd = transreson(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 1);
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j+1] + ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * ( c_intrinsic + c_load);
		}
		this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
		ptr_htree_at_mat_interval->delay_nand_driving_inv_final_seg += this_delay;
		inrisetime = this_delay/(1.0 - 0.5);
	}

	//Find delay of the H-tree
	for(j = 0; j < number_htree_nodes - 1; ++j){
		number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_htree_nodes - j) + 0.5)) / 2, 1);
		for(i = 0; i < number_mats_to_cover_in_this_htree_segment; ++i){
			if(i==0){//Driver structure is NAND2 driving inverter load
				*delay += ptr_htree_at_mat_interval->delay_nand_driving_inv;
				}
			else if(i==number_mats_to_cover_in_this_htree_segment - 1){//Driver structure is inverter
					//driving NAND2 load
				*delay += ptr_htree_at_mat_interval->delay_inv_driving_nand;
				}
			else{//Driver structure is inverter driving inverter load
				*delay += ptr_htree_at_mat_interval->delay_inv_driving_inv;
			}
		}
	}
	*delay += ptr_htree_at_mat_interval->delay_nand_driving_inv_final_seg;
}
	

void calculate_power_htree_with_nodes_at_mat_interval(addr_datain_htree_at_mat_interval 
													   *ptr_htree_at_mat_interval)
{
	int j;
    double c_load, c_intrinsic;
   
    //Find power of NAND2 driving inverter
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_nand_driving_inv; ++j){
		if(j == 0){//NAND2 gate
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = 2 * draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.dynamic += 
				(c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_nand_driving_inv - 1){
			c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv + 
				ptr_htree_at_mat_interval->c_wire_load;
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.leakage += cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else{//inverter
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL)
				+ draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.dynamic += (c_intrinsic +
				c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_nand_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j],
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
	}

	//Find power of inverter driving inverter
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_inv_driving_inv; ++j){
		if(j == ptr_htree_at_mat_interval->number_gates_inv_driving_inv - 1){
			c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_inv + 
				ptr_htree_at_mat_interval->c_wire_load;
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.dynamic += (c_intrinsic + 
				c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], 
					ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], temper) * 
					0.5 * vdd_periph_global;
		}
		else{//inverter
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.dynamic += (c_intrinsic +
				c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_width_n[j], 
                            ptr_htree_at_mat_interval->nand_driving_inv_width_p[j], temper) *
							0.5 * vdd_periph_global;
		}
	}

	//Find power of inverter driving NAND2
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_inv_driving_nand; ++j){
		if(j == ptr_htree_at_mat_interval->number_gates_inv_driving_nand - 1){//NAND2 gate
			c_load = ptr_htree_at_mat_interval->c_gate_load_inv_driving_nand + 
				ptr_htree_at_mat_interval->c_wire_load;
			c_intrinsic = 2 * draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_inv_driving_nand.readOp.dynamic += (c_intrinsic + 
				c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_inv_driving_nand.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], 
				ptr_htree_at_mat_interval->inv_driving_nand_width_p[j], temper) * 0.5 *
 vdd_periph_global;
		}
		else{//inverter
		    c_load = gatecap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j+1] + 
				ptr_htree_at_mat_interval->inv_driving_nand_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_inv_driving_nand.readOp.dynamic += (c_intrinsic + 
				c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_inv_driving_nand.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->inv_driving_nand_width_n[j], 
				ptr_htree_at_mat_interval->inv_driving_nand_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
	}
	

	//Find power of NAND2 driving inverter for final segment of H-tree (over half a mat height/width)
    for(j = 0; j < ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg; ++j){
		if(j == 0){//NAND2 gate
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j+1] + 
				ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j+1], 0.0);
			c_intrinsic = 2 * draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic += 
				(c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], 
				ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_nand_driving_inv_final_seg - 1){
			c_load = ptr_htree_at_mat_interval->c_gate_load_nand_driving_inv + 
				ptr_htree_at_mat_interval->c_wire_load_final_seg;
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage += cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], 
				ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else{//inverter
			c_load = gatecap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j+1] + ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage += cmos_ileakage(ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_n[j], 
				ptr_htree_at_mat_interval->nand_driving_inv_final_seg_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
	}

	//Find power of tristate inverter driving inverter load.
	for(j = 0; j < ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv; ++j){
		if(j == ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 2){//NAND2 gate
			c_load = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL) + 
				2 * draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic += 
				(c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], temper) * 0.5 *
				vdd_periph_global;
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n, 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p, temper) * 0.5 *
				vdd_periph_global;
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 1){//PMOS 
			c_load = ptr_htree_at_mat_interval->c_wire_load + ptr_htree_at_mat_interval->c_gate_load_tristate_driver_driving_inv +
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic += 
				(c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 3){//inverter driving NAND2 and NOR2
			c_load = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0.0) +
				gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n + ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p, 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic += 
				(c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
		else{//inverter
			c_load = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic += 
				(c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage += 
				cmos_ileakage(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], temper) * 0.5 * vdd_periph_global;
		}
	}
}



void power_addr_datain_htree_node_at_mat_interval(addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval, 
													  int number_mats_to_cover_in_this_htree_segment,
													  powerDef *power)
{
	int i;

	power->readOp.dynamic = 0;
	power->readOp.leakage = 0;
	power->writeOp.dynamic = 0;
	power->writeOp.leakage = 0;

	if(number_mats_to_cover_in_this_htree_segment == 1){//final segment in H-tree - going over half a 
		//subarray
		power->readOp.dynamic += ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic;
		power->readOp.leakage += ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage;
    }
    else{
	    for(i = 0; i < number_mats_to_cover_in_this_htree_segment; ++i){
		    if(i==0){//Driver structure is NAND2 driving inverter load
				power->readOp.dynamic += ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic;
				power->readOp.leakage += ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage;
		    }
		    else if(i==number_mats_to_cover_in_this_htree_segment - 1){//Driver structure is inverter
			    //driving NAND2 load
				power->readOp.dynamic += ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.dynamic;
				power->readOp.leakage += ptr_htree_at_mat_interval->power_nand_driving_inv_final_seg.readOp.leakage;
		    }
		    else{//Driver structure is inverter driving inverter load
			    power->readOp.dynamic += ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.dynamic;
				power->readOp.leakage += ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.leakage;
		    }
	    }
    }
}


void power_dataout_htree_node_at_mat_interval(addr_datain_htree_at_mat_interval *ptr_htree_at_mat_interval,
											  int number_mats_to_cover_in_this_htree_segment,
											  powerDef *power)
{
	int i;

	power->readOp.dynamic = 0;
	power->readOp.leakage = 0;

	for(i = 0; i < number_mats_to_cover_in_this_htree_segment; ++i){
		if(i==0){//Driver structure is  tristate driver driving inverter load
			power->readOp.dynamic += ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic;
			power->readOp.leakage += ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage;
		}
		else{//Driver structure is inverter driving inverter load
			power->readOp.dynamic += ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.dynamic;
			power->readOp.leakage += ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.leakage;
		}
	}
}


void delay_dataout_vertical_htree(dataout_htree_node *ptr_htree_node, int number_htree_nodes, 
        int hhtree_nodes, double inrisetime, double *outrisetime, double *delay, 
        powerDef *power, input_params_t *parameters, area_type *data_array)
{
    int i, j, htree_node_multiplier;
    double c_load, rd, tf, this_delay, c_intrinsic;
    double len1, size1, tc, nand_d = 0;
//	pda_res_t ls_pda;
    wire_stats_t wst;
    this_delay = 0;
    *delay = 0;
    *outrisetime = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;
    if (parameters->wire_inter_mats == Low_swing) { /* if we are using low swing instead of htree */
        rd = (data_array->height/2 + data_array->width);
        wst.wt = Low_swing;
        wst.wire_length = rd*1e-6; //m
        wst.nsense = number_htree_nodes/2;
        calc_wire_stats2 (Low_swing, &wst);
        *delay += wst.wire_pda.delay;
        power->readOp.dynamic = wst.wire_pda.power.dynamic;
        power->readOp.leakage = wst.wire_pda.power.leakage;
        *outrisetime = signal_rise_time();
        return;
    }
#ifdef USE_GLOBAL 
	else {
		rd = data_array->height/2 + data_array->width;
        wst.wt = parameters->wire_inter_mats;
        wst.wire_length = rd*1e-6; //m
        calc_wire_stats2 (parameters->wire_inter_mats, &wst);
        *delay += wst.wire_pda.delay;
		rd = data_array->height/2 + data_array->width*3/2;
        wst.wire_length = rd*1e-6; //m
        calc_wire_stats2 (parameters->wire_inter_mats, &wst);
        power->readOp.dynamic = wst.wire_pda.power.dynamic;
        power->readOp.leakage = wst.wire_pda.power.leakage;
        /* delay due to the nand gates in the htree */
        i = hhtree_nodes-1;
        assert (i>=0);
        j = 4;
        if (i != 0) {
            while (i!=0) {
                len1 = data_array->width/j;
                if (wst.repeater_spacing > len1) {
                    size1 = wst.repeater_size*len1/wst.repeater_spacing;
                }
                else {
                    size1 = wst.repeater_size;
                }
                tc = 2*transreson((size1/4)*minimum_width_nmos, NCH, 1) *
                    draincap((size1/4)*minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL)*2;
                nand_d += horowitz (signal_rise_time(), tc, 0.5, 0.5, RISE);
                power->readOp.dynamic += 0.5*draincap((size1/4)*minimum_width_nmos, 
                        NCH, 1, 1, DEFAULTHEIGHTCELL)*2*vdd_periph_global;
                i--;
                j*=2;
            }
        }
        if (i != 0) {
            i = number_htree_nodes-1;
            j=2;
            while (i!=1) {
                len1 =  data_array->height/j;
                if (wst.repeater_spacing > len1) {
                    size1 = wst.repeater_size*len1/wst.repeater_spacing;
                }
                else {
                    size1 = wst.repeater_size;
                }
                tc = 2*transreson((size1/4)*minimum_width_nmos, NCH, 1) *
                    draincap((size1/4)*minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL)*2;
                nand_d += horowitz (signal_rise_time(), tc, 0.5, 0.5, RISE);
                power->readOp.dynamic += 0.5*draincap((size1/4)*minimum_width_nmos, 
                        NCH, 1, 1, DEFAULTHEIGHTCELL)*2*vdd_periph_global;
                i--;
                j*=2;
            }
        }
        *delay += nand_d;
        *outrisetime = signal_rise_time();
//        ls_pda.area_stats.height = rd*1e-6;//(m)
//        calc_wire_stats (Global, &ls_pda, 1);
//        calc_wire_stats2 (Global, &ls_pda, 1);
//        *delay += ls_pda.delay;
//        power->readOp.dynamic = ls_pda.power.dynamic;
//        power->readOp.leakage = ls_pda.power.leakage;
//        *outrisetime = signal_rise_time();
        return;
	}
#endif


	for(i = number_htree_nodes - 2; i >= 0; --i){
        for(j = 0; j < ptr_htree_node[i].number_gates; ++j){
            if(j == 0){//NAND2 gate
                rd = transreson(ptr_htree_node[i].width_n[j], NCH, 2);
				if(ptr_htree_node[i].number_gates == 2){
					c_load = gatecap(ptr_htree_node[i].width_p[j+1], 0.0);//NAND2 driving PMOS output stage
				}
				else{
					c_load = gatecap(ptr_htree_node[i].width_p[j+1], 0.0);//NAND2 driving inverter
				}
                c_intrinsic = draincap(ptr_htree_node[i].width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL) + 
                    2 * draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }
            else if(j == ptr_htree_node[i].number_gates - 1){//PMOS 
				rd = transreson(ptr_htree_node[i].width_p[j], PCH, 1);
				c_load = ptr_htree_node[i].c_wire_load + ptr_htree_node[i].c_gate_load +
					draincap(ptr_htree_node[i].width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL) + 
					draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
				c_intrinsic = draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
					draincap(ptr_htree_node[i].width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
				tf = rd * (c_intrinsic + c_load) + ptr_htree_node[i].r_wire_load * 
					(ptr_htree_node[i].c_wire_load / 2 + ptr_htree_node[i].c_gate_load +
					draincap(ptr_htree_node[i].width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL) + 
					draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL));
				power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			}
			else{//inverter
				rd = transreson(ptr_htree_node[i].width_n[j], NCH, 1);
				if(j == ptr_htree_node[i].number_gates - 2){//inverter driving PMOS of output stage
					c_load = gatecap(ptr_htree_node[i].width_p[j+1], 0.0);
				}
				else{//inverter driving inverter
					c_load = gatecap(ptr_htree_node[i].width_n[j+1] + ptr_htree_node[i].width_p[j+1], 0.0);
				}				
				c_intrinsic = draincap(ptr_htree_node[i].width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) +
					draincap(ptr_htree_node[i].width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
				tf = rd * (c_intrinsic + c_load);
				power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			}

			this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
            *delay += this_delay;
            inrisetime = this_delay/(1.0 - 0.5);
        }
    }
    *outrisetime = inrisetime;

    htree_node_multiplier = 1;
    for(i = 0; i < number_htree_nodes - 1; ++i){
        for(j = 0; j < ptr_htree_node[i].number_gates; ++j){
            if(j == 0){//NAND gate, NOR gate
                power->readOp.leakage += cmos_ileakage(ptr_htree_node[i].width_n[j], 
                        ptr_htree_node[i].width_p[j], temper) * 0.5 * vdd_periph_global * 
						htree_node_multiplier;
				power->readOp.leakage += cmos_ileakage(ptr_htree_node[i].width_nor2_n, 
					ptr_htree_node[i].width_nor2_p, temper) * 0.5 * vdd_periph_global * htree_node_multiplier;
            }
            else if(j == ptr_htree_node[i].number_gates - 1){//PMOS and NMOS output stage
				power->readOp.leakage += cmos_ileakage(ptr_htree_node[i].width_n[j], 
					ptr_htree_node[i].width_p[j], temper) * 0.5 * vdd_periph_global * htree_node_multiplier;
			}
			else{//inverters in NAND2 and NOR2 paths. Assume that inverters in the NOR2 path 
				//half the width of corresp inverters in NAND2 path (because inverters in NOR2 path
				//drive NMOS of output stage compared to inverters in NAND2 path that drive PMOS of 
				//output stage)
				power->readOp.leakage += cmos_ileakage(ptr_htree_node[i].width_n[j] + ptr_htree_node[i].width_n[j] / 2,
					ptr_htree_node[i].width_p[j] + ptr_htree_node[i].width_p[j] / 2, temper) * 
					0.5 * vdd_periph_global * htree_node_multiplier;
			}
			htree_node_multiplier *= 2;
		}
	}
}

void delay_datout_vertical_htree_with_nodes_at_mat_interval(addr_datain_htree_at_mat_interval 
													   *ptr_htree_at_mat_interval, 
													   int number_htree_nodes, 
													   double inrisetime, double *outrisetime, 
													   double *delay, powerDef *power)
{
	int j, i, number_mats_to_cover_in_this_htree_segment, htree_seg_multiplier;
    double rd, c_load, c_intrinsic, tf, this_delay;

	this_delay = 0;
    *delay = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;
    
	for(j = 0; j < ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv; ++j){
		if(j == ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 2){//NAND2 gate
			rd = transreson(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 2);
			c_load = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + 
				ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL) + 
				2 * draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load);
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 1){//PMOS 
			rd = transreson(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1);
			c_load = ptr_htree_at_mat_interval->c_wire_load + ptr_htree_at_mat_interval->c_gate_load_tristate_driver_driving_inv +
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) +
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load) + ptr_htree_at_mat_interval->r_wire_load * 
				ptr_htree_at_mat_interval->c_wire_load / 2;
		}
		else if(j == ptr_htree_at_mat_interval->number_gates_tristate_driver_driving_inv - 3){//inverter driving NAND2 and NOR2
			rd = transreson(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1);
			c_load = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0.0) +
				gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_n + ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_nor2_p, 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) +
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load);
		}
		else{//inverter
			rd = transreson(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1);
			c_load = gatecap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j+1] + ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j+1], 0.0);
			c_intrinsic = draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_at_mat_interval->tristate_driver_driving_inv_width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load);
		}
		this_delay += horowitz(inrisetime, tf, 0.5, 0.5, RISE);
		inrisetime = this_delay/(1.0 - 0.5);
	}
 	ptr_htree_at_mat_interval->delay_tristate_driver_driving_inv = this_delay;

	htree_seg_multiplier = 1;
	for(i = 0; i < number_htree_nodes - 1; ++i){
		number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_htree_nodes - 1 - i) + 0.5)) / 2, 1);
		for(j = 0; j < number_mats_to_cover_in_this_htree_segment; ++j){
			if(j==0){//Driver structure is  tristate driver driving inverter load
				*delay += ptr_htree_at_mat_interval->delay_tristate_driver_driving_inv;
				power->readOp.dynamic += ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.dynamic;
				power->readOp.leakage += ptr_htree_at_mat_interval->power_tristate_driver_driving_inv.readOp.leakage *
					htree_seg_multiplier;
			}
			else{//Driver structure is inverter driving inverter load
				*delay += ptr_htree_at_mat_interval->delay_inv_driving_inv;
				power->readOp.dynamic += ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.dynamic;
				power->readOp.leakage += ptr_htree_at_mat_interval->power_inv_driving_inv.readOp.leakage *
					htree_seg_multiplier;
			}
		}
		htree_seg_multiplier *= 2;
	}
}

void initialize_point_to_point_interconnect_segment(
        point_to_point_interconnect_segment *ptr_intcnt_seg, double length, double c_gate_load)
{
    ptr_intcnt_seg->length = length;
    ptr_intcnt_seg->optimal_repeater_size = 0;
    ptr_intcnt_seg->c_gate_load = c_gate_load;
    ptr_intcnt_seg->minimum_number_gates = 1;
    ptr_intcnt_seg->power.readOp.dynamic = 0;
    ptr_intcnt_seg->power.readOp.leakage = 0;
    ptr_intcnt_seg->width_n[0] = minimum_width_nmos;
    ptr_intcnt_seg->width_p[0] = 2 * minimum_width_nmos;
    ptr_intcnt_seg->delay = 0;
}

void initialize_decoder(int num_decode_signals, int flag_way_select, double c_load_decoder_output, 
        double r_wire_decoder_output, decoder *ptr_dec)
{
    int number_addr_bits_decode, i;
    ptr_dec->flag_decoder_exists = 0;
    ptr_dec->min_number_gates = 2;
    ptr_dec->delay = 0;
    ptr_dec->power.readOp.dynamic = 0;
    ptr_dec->power.readOp.leakage = 0;
    ptr_dec->power.writeOp.dynamic = 0;
    ptr_dec->power.writeOp.leakage = 0;
    for(i = 0; i < MAX_NUMBER_GATES_STAGE; ++i){
        ptr_dec->width_decoder_n[i] = 0;
        ptr_dec->width_decoder_p[i] = 0;
    }
    number_addr_bits_decode = (int) (logbasetwo((double) (num_decode_signals)));	
	if(number_addr_bits_decode < 4){
		if(flag_way_select)
		{
			ptr_dec->flag_decoder_exists = 1;
			ptr_dec->number_input_signals = 2;
		}
		else{
			ptr_dec->number_input_signals = 0;
		}
    }
    else{
		if(flag_way_select)
		{
			ptr_dec->flag_decoder_exists = 1;
			ptr_dec->number_input_signals = 3;
		}
		else{
			ptr_dec->flag_decoder_exists = 1;
			ptr_dec->number_input_signals = 2;
		}
    }	
    ptr_dec->r_wire_decoder_output = r_wire_decoder_output;
    ptr_dec->c_load_decoder_output = c_load_decoder_output;
}

void initialize_predecoder_blocks(int num_decode_signals, predecoder_block *ptr_predec_blk_1,
        predecoder_block *ptr_predec_blk_2, decoder *ptr_dec,
        double c_wire_predecoder_block_output, 
        double r_wire_predecoder_block_output,
        int number_decode_gates_driven_per_predecode_output)
{
    int number_addr_bits_decode, branch_effort_predecoder_block_1_output,
    branch_effort_predecoder_block_2_output, i;
    double c_load_decoder_gate;
    number_addr_bits_decode = (int) (logbasetwo((double) (num_decode_signals)));
    ptr_predec_blk_1->flag_block_exists = 0;
    ptr_predec_blk_2->flag_block_exists = 0;
    ptr_predec_blk_1->number_input_addr_bits = 0;
    ptr_predec_blk_2->number_input_addr_bits = 0;
    ptr_predec_blk_1->number_gates_first_level_nand2_path = 0;
    ptr_predec_blk_1->number_gates_first_level_nand3_path = 0;
    ptr_predec_blk_1->number_gates_second_level = 0;
    ptr_predec_blk_2->number_gates_first_level_nand2_path = 0;
    ptr_predec_blk_2->number_gates_first_level_nand3_path = 0;
    ptr_predec_blk_2->number_gates_second_level = 0;
    ptr_predec_blk_1->min_number_gates_first_level = 2;
    ptr_predec_blk_1->min_number_gates_second_level = 2;
    ptr_predec_blk_2->min_number_gates_first_level = 2;
    ptr_predec_blk_2->min_number_gates_second_level = 2;
    ptr_predec_blk_1->delay_nand2_path = 0;
    ptr_predec_blk_1->delay_nand3_path = 0;
    ptr_predec_blk_2->delay_nand2_path = 0;
    ptr_predec_blk_2->delay_nand3_path = 0;
    ptr_predec_blk_1->power_nand2_path.readOp.dynamic = 0;
    ptr_predec_blk_1->power_nand3_path.readOp.dynamic = 0;
    ptr_predec_blk_1->power_second_level.readOp.dynamic = 0;
    ptr_predec_blk_1->power_nand2_path.readOp.leakage = 0;
    ptr_predec_blk_1->power_nand3_path.readOp.leakage = 0;
    ptr_predec_blk_1->power_second_level.readOp.leakage = 0;
    ptr_predec_blk_2->power_nand2_path.readOp.dynamic = 0;
    ptr_predec_blk_2->power_nand3_path.readOp.dynamic = 0;
    ptr_predec_blk_2->power_second_level.readOp.dynamic = 0;
    ptr_predec_blk_2->power_nand2_path.readOp.leakage = 0;
    ptr_predec_blk_2->power_nand3_path.readOp.leakage = 0;
    ptr_predec_blk_2->power_second_level.readOp.leakage = 0;
    ptr_predec_blk_1->power_nand2_path.writeOp.dynamic = 0;
    ptr_predec_blk_1->power_nand3_path.writeOp.dynamic = 0;
    ptr_predec_blk_1->power_second_level.writeOp.dynamic = 0;
    ptr_predec_blk_1->power_nand2_path.writeOp.leakage = 0;
    ptr_predec_blk_1->power_nand3_path.writeOp.leakage = 0;
    ptr_predec_blk_1->power_second_level.writeOp.leakage = 0;
    ptr_predec_blk_2->power_nand2_path.writeOp.dynamic = 0;
    ptr_predec_blk_2->power_nand3_path.writeOp.dynamic = 0;
    ptr_predec_blk_2->power_second_level.writeOp.dynamic = 0;
    ptr_predec_blk_2->power_nand2_path.writeOp.leakage = 0;
    ptr_predec_blk_2->power_nand3_path.writeOp.leakage = 0;
    ptr_predec_blk_2->power_second_level.writeOp.leakage = 0;
    for(i = 0; i < MAX_NUMBER_GATES_STAGE; ++i){
        ptr_predec_blk_1->width_first_level_nand2_path_n[i] = 0;
        ptr_predec_blk_1->width_first_level_nand2_path_p[i] = 0;
        ptr_predec_blk_1->width_first_level_nand3_path_n[i] = 0;
        ptr_predec_blk_1->width_first_level_nand3_path_p[i] = 0;
        ptr_predec_blk_2->width_first_level_nand2_path_n[i] = 0;
        ptr_predec_blk_2->width_first_level_nand2_path_p[i] = 0;
        ptr_predec_blk_2->width_first_level_nand3_path_n[i] = 0;
        ptr_predec_blk_2->width_first_level_nand3_path_p[i] = 0;
        ptr_predec_blk_2->width_second_level_n[i] = 0;
        ptr_predec_blk_2->width_second_level_p[i] = 0;
    }

    if((number_addr_bits_decode != 0)){
        if((number_addr_bits_decode == 1)||(number_addr_bits_decode == 2)){//Just one predecoder block is
            //required with NAND2 gates. No decoder required. The first level of predecoding directly drives 
            //the decoder output load
            ptr_predec_blk_1->flag_block_exists = 1;
            ptr_predec_blk_1->number_input_addr_bits = number_addr_bits_decode;
            ptr_predec_blk_2->number_input_addr_bits = 0;
            ptr_predec_blk_1->r_wire_predecoder_block_output = ptr_dec->r_wire_decoder_output;
            ptr_predec_blk_2->r_wire_predecoder_block_output = 0;
            ptr_predec_blk_1->c_load_predecoder_block_output = ptr_dec->c_load_decoder_output ;	
            ptr_predec_blk_2->c_load_predecoder_block_output = 0;	
        }
        else{
            if(number_addr_bits_decode == 3){//Just one predecoder block is required with NAND3 gates. No decoder required. The first level 
                    //of predecoding durectly drives the decoder output load.
				ptr_predec_blk_1->flag_block_exists = 1;
				ptr_predec_blk_1->number_input_addr_bits = 3;
				ptr_predec_blk_2->number_input_addr_bits = 0;
				ptr_predec_blk_1->r_wire_predecoder_block_output = ptr_dec->r_wire_decoder_output;
				ptr_predec_blk_2->r_wire_predecoder_block_output = 0;
                    ptr_predec_blk_1->c_load_predecoder_block_output = ptr_dec->c_load_decoder_output;	
                    ptr_predec_blk_2->c_load_predecoder_block_output = 0;	
            }
            else{
                ptr_predec_blk_1->flag_block_exists = 1;
                ptr_predec_blk_2->flag_block_exists = 1;
                if(number_addr_bits_decode%2 == 0){
                    ptr_predec_blk_1->number_input_addr_bits = number_addr_bits_decode / 2;
                }
                else{
                    ptr_predec_blk_1->number_input_addr_bits = (int) (number_addr_bits_decode / 2) + 1;
                }
                ptr_predec_blk_2->number_input_addr_bits = number_addr_bits_decode - ptr_predec_blk_1->number_input_addr_bits;
                branch_effort_predecoder_block_1_output = (int) (pow(2, ptr_predec_blk_2->number_input_addr_bits));
                branch_effort_predecoder_block_2_output = (int) (pow(2, ptr_predec_blk_1->number_input_addr_bits));
                c_load_decoder_gate = number_decode_gates_driven_per_predecode_output *
                    gatecap(ptr_dec->width_decoder_n[0] + ptr_dec->width_decoder_p[0], 0);
                ptr_predec_blk_1->r_wire_predecoder_block_output = r_wire_predecoder_block_output;
                ptr_predec_blk_2->r_wire_predecoder_block_output = r_wire_predecoder_block_output;
                ptr_predec_blk_1->c_load_predecoder_block_output = 
                    branch_effort_predecoder_block_1_output * c_load_decoder_gate + c_wire_predecoder_block_output;	
                ptr_predec_blk_2->c_load_predecoder_block_output = 
                    branch_effort_predecoder_block_2_output * c_load_decoder_gate + c_wire_predecoder_block_output;	
            }
        }
    }
}

void initialize_predecoder_block_drivers(int number_decode_signals, int flag_way_select, int way_select,
										 predecoder_block_driver *ptr_predec_blk_driver_1,
										 predecoder_block_driver *ptr_predec_blk_driver_2,
										 predecoder_block *ptr_predec_blk_1,
										 predecoder_block *ptr_predec_blk_2, decoder *ptr_dec)
{
    int number_addr_bits_decode, i;
    number_addr_bits_decode = (int) (logbasetwo((double) (number_decode_signals)));
    ptr_predec_blk_driver_1->flag_driver_exists = 0;
    ptr_predec_blk_driver_2->flag_driver_exists = 0;
    ptr_predec_blk_driver_1->flag_driving_decoder_output = 0;
    ptr_predec_blk_driver_2->flag_driving_decoder_output = 0;
    ptr_predec_blk_driver_1->number_input_addr_bits = ptr_predec_blk_1->number_input_addr_bits;
	ptr_predec_blk_driver_2->number_input_addr_bits = ptr_predec_blk_2->number_input_addr_bits;
    ptr_predec_blk_driver_1->number_gates_nand2_path = 0;
    ptr_predec_blk_driver_2->number_gates_nand2_path = 0;
    ptr_predec_blk_driver_1->number_gates_nand3_path = 0;
    ptr_predec_blk_driver_2->number_gates_nand3_path = 0;
    ptr_predec_blk_driver_1->number_parallel_instances_driving_1_nand2_load = 0;
    ptr_predec_blk_driver_1->number_parallel_instances_driving_2_nand2_load = 0;
    ptr_predec_blk_driver_1->number_parallel_instances_driving_4_nand2_load = 0;
    ptr_predec_blk_driver_1->number_parallel_instances_driving_2_nand3_load = 0;
    ptr_predec_blk_driver_1->number_parallel_instances_driving_8_nand3_load = 0;
    ptr_predec_blk_driver_2->number_parallel_instances_driving_1_nand2_load = 0;
    ptr_predec_blk_driver_2->number_parallel_instances_driving_2_nand2_load = 0;
    ptr_predec_blk_driver_2->number_parallel_instances_driving_4_nand2_load = 0;
    ptr_predec_blk_driver_2->number_parallel_instances_driving_2_nand3_load = 0;
    ptr_predec_blk_driver_2->number_parallel_instances_driving_8_nand3_load = 0;
    ptr_predec_blk_driver_1->min_number_gates = 2;
    ptr_predec_blk_driver_2->min_number_gates = 2;
    ptr_predec_blk_driver_1->delay_nand2_path = 0;
    ptr_predec_blk_driver_1->delay_nand3_path = 0;
    ptr_predec_blk_driver_2->delay_nand2_path = 0;
    ptr_predec_blk_driver_2->delay_nand3_path = 0;
    ptr_predec_blk_driver_1->power_nand2_path.readOp.dynamic = 0;
    ptr_predec_blk_driver_1->power_nand3_path.readOp.dynamic = 0;
    ptr_predec_blk_driver_1->power_nand2_path.readOp.leakage = 0;
    ptr_predec_blk_driver_1->power_nand3_path.readOp.leakage = 0;
    ptr_predec_blk_driver_1->power_nand2_path.writeOp.dynamic = 0;
    ptr_predec_blk_driver_1->power_nand3_path.writeOp.dynamic = 0;
    ptr_predec_blk_driver_1->power_nand2_path.writeOp.leakage = 0;
    ptr_predec_blk_driver_1->power_nand3_path.writeOp.leakage = 0;
    ptr_predec_blk_driver_2->power_nand2_path.readOp.dynamic = 0;
    ptr_predec_blk_driver_2->power_nand3_path.readOp.dynamic = 0;
    ptr_predec_blk_driver_2->power_nand2_path.readOp.leakage = 0;
    ptr_predec_blk_driver_2->power_nand3_path.readOp.leakage = 0;
    ptr_predec_blk_driver_2->power_nand2_path.writeOp.dynamic = 0;
    ptr_predec_blk_driver_2->power_nand3_path.writeOp.dynamic = 0;
    ptr_predec_blk_driver_2->power_nand2_path.writeOp.leakage = 0;
    ptr_predec_blk_driver_2->power_nand3_path.writeOp.leakage = 0;
    for(i = 0; i < MAX_NUMBER_GATES_STAGE; ++i){
        ptr_predec_blk_driver_1->width_nand2_path_n[i] = 0;
        ptr_predec_blk_driver_1->width_nand2_path_p[i] = 0;
        ptr_predec_blk_driver_1->width_nand3_path_n[i] = 0;
        ptr_predec_blk_driver_1->width_nand3_path_p[i] = 0;
        ptr_predec_blk_driver_2->width_nand2_path_n[i] = 0;
        ptr_predec_blk_driver_2->width_nand2_path_p[i] = 0;
        ptr_predec_blk_driver_2->width_nand3_path_n[i] = 0;
        ptr_predec_blk_driver_2->width_nand3_path_p[i] = 0;
    }
    ptr_predec_blk_driver_1->r_load_nand2_path_predecode_block_driver_output = 0;
    ptr_predec_blk_driver_1->r_load_nand3_path_predecode_block_driver_output = 0;
    ptr_predec_blk_driver_2->r_load_nand2_path_predecode_block_driver_output = 0;
    ptr_predec_blk_driver_2->r_load_nand3_path_predecode_block_driver_output = 0;

	if(flag_way_select && (way_select > 1)){
		ptr_predec_blk_driver_1->flag_driver_exists = 1;
		ptr_predec_blk_driver_1->number_input_addr_bits = way_select; 
		if(ptr_dec->number_input_signals == 2){
			 ptr_predec_blk_driver_1->c_load_nand2_path_predecode_block_driver_output = 
				 gatecap(ptr_dec->width_decoder_n[0] + ptr_dec->width_decoder_p[0], 0);
			 ptr_predec_blk_driver_1->number_parallel_instances_driving_2_nand2_load =
				 ptr_predec_blk_driver_1->number_input_addr_bits;
		}
		else{
			if(ptr_dec->number_input_signals == 3){
				ptr_predec_blk_driver_1->c_load_nand3_path_predecode_block_driver_output = 
					gatecap(ptr_dec->width_decoder_n[0] + ptr_dec->width_decoder_p[0], 0);
				ptr_predec_blk_driver_1->number_parallel_instances_driving_2_nand3_load =
					ptr_predec_blk_driver_1->number_input_addr_bits;
			}
		}
        ptr_predec_blk_driver_1->r_load_nand2_path_predecode_block_driver_output = 0;
		ptr_predec_blk_driver_1->c_load_nand3_path_predecode_block_driver_output = 0;
        ptr_predec_blk_driver_1->r_load_nand3_path_predecode_block_driver_output = 0;
        ptr_predec_blk_driver_2->c_load_nand2_path_predecode_block_driver_output = 0;
        ptr_predec_blk_driver_2->r_load_nand2_path_predecode_block_driver_output = 0;
        ptr_predec_blk_driver_2->c_load_nand3_path_predecode_block_driver_output = 0;
        ptr_predec_blk_driver_2->r_load_nand3_path_predecode_block_driver_output = 0;
	}
	if(!flag_way_select){
		/*if(number_addr_bits_decode == 0){//means that there is no decoding required, however the 
			//load still needs to be driven by a chain of inverters. Use only one predecode block
			//driver
			ptr_predec_blk_driver_1->flag_driver_exists = 1;
			ptr_predec_blk_driver_1->flag_driving_decoder_output = 1;
			ptr_predec_blk_driver_1->c_load_nand2_path_predecode_block_driver_output = 
				ptr_dec->c_load_decoder_output;
			ptr_predec_blk_driver_1->r_load_nand2_path_predecode_block_driver_output =
				ptr_dec->r_wire_decoder_output;
			ptr_predec_blk_driver_1->c_load_nand3_path_predecode_block_driver_output = 0;
			ptr_predec_blk_driver_1->r_load_nand3_path_predecode_block_driver_output = 0;
			ptr_predec_blk_driver_2->c_load_nand2_path_predecode_block_driver_output = 0;
			ptr_predec_blk_driver_2->r_load_nand2_path_predecode_block_driver_output = 0;
			ptr_predec_blk_driver_2->c_load_nand3_path_predecode_block_driver_output = 0;
			ptr_predec_blk_driver_2->r_load_nand3_path_predecode_block_driver_output = 0;
		}
		else{*/
			if(ptr_predec_blk_1->flag_block_exists){
				ptr_predec_blk_driver_1->flag_driver_exists = 1;
			}
			if(ptr_predec_blk_2->flag_block_exists){
				ptr_predec_blk_driver_2->flag_driver_exists = 1;
			}
		//}
	}
}

void
compute_widths_predecoder_block_driver(predecoder_block_driver *ptr_predec_blk_driver, 
        predecoder_block *ptr_predec_blk, decoder *ptr_dec, int way_select)
{
    //The predecode block driver accepts as input the address bits from the h-tree network. For 
    //each addr bit it then generates addr and addrbar as outputs. For now ignore the effect of
    //inversion to generate addrbar and simply treat addrbar as addr.

    int  i;
    double F, f;
    if(ptr_predec_blk_driver->flag_driver_exists){	
		if(ptr_predec_blk_driver->flag_driving_decoder_output){
			ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 1;
            ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
            ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 0;
            ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
            ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 0;
            ptr_predec_blk_driver->width_nand2_path_n[0] = minimum_width_nmos;
            ptr_predec_blk_driver->width_nand2_path_p[0] = 2 * ptr_predec_blk_driver->width_nand2_path_n[0];
            F = ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output / gatecap(ptr_predec_blk_driver->width_nand2_path_n[0] +
                    ptr_predec_blk_driver->width_nand2_path_p[0], 0);
            ptr_predec_blk_driver->number_gates_nand2_path = (int) (log(F) / log(fopt) + 0.5);
            //Check if number_gates_second_level_predecoder is odd. If it is, add 1 to make it even.
            if(ptr_predec_blk_driver->number_gates_nand2_path%2 != 0) {
                --ptr_predec_blk_driver->number_gates_nand2_path;
            }
            if(ptr_predec_blk_driver->number_gates_nand2_path < ptr_predec_blk_driver->min_number_gates){
                ptr_predec_blk_driver->number_gates_nand2_path = ptr_predec_blk_driver->min_number_gates;
            }
            //Recalculate the effective fanout of each stage for the above number_gates_last_level_predecoder
            f = pow(F, 1.0 / ptr_predec_blk_driver->number_gates_nand2_path);
            i = ptr_predec_blk_driver->number_gates_nand2_path - 1;
            ptr_predec_blk_driver->width_nand2_path_n[i] = ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output / (gatecap(1, 0) * f);
            ptr_predec_blk_driver->width_nand2_path_p[i]  = 2 * ptr_predec_blk_driver->width_nand2_path_n[i];
            for(i = ptr_predec_blk_driver->number_gates_nand2_path - 2; i >= 1; --i){
                ptr_predec_blk_driver->width_nand2_path_n[i] = ptr_predec_blk_driver->width_nand2_path_n[i+1] / f;
                ptr_predec_blk_driver->width_nand2_path_p[i] = 2 * ptr_predec_blk_driver->width_nand2_path_n[i];
            }
        }
        else{
			if(way_select == 0){
				if(ptr_predec_blk->number_input_addr_bits == 1){//2 NAND2 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load =	1;
					ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output = 
						2 * gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0], 0);

				}
				if(ptr_predec_blk->number_input_addr_bits == 2){//4 NAND2 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 2;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 0;
					ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output = 
						4 * gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 3){//8 NAND3 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 3;
					ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output = 
						8 * gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 4){//4 + 4 NAND2 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 4;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 0;
					ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output = 
						4 * gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 5){//4 NAND2 gates, 8 NAND3 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 2;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 3;
					ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output =
						4 * gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0], 0);
					ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output = 
						8 * gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 6){//8 + 8 NAND3 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 6;
					ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output = 
						8 * gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 7){//4 + 4 NAND2 gates, 8 NAND3 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 4;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 3;
					ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output =
						4 * gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0], 0);
					ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output = 
						8 * gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 8){//4 NAND2 gates, 8 + 8 NAND3 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 2;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 6;
					ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output =
						4 * gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0], 0);
					ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output = 
						8 * gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0], 0);
				}
				if(ptr_predec_blk->number_input_addr_bits == 9){//8 + 8 + 8 NAND3 gates
					ptr_predec_blk_driver->number_parallel_instances_driving_1_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_4_nand2_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_2_nand3_load = 0;
					ptr_predec_blk_driver->number_parallel_instances_driving_8_nand3_load = 9;
					ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output = 
						8 * gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0], 0);
				}
			}

			if((ptr_predec_blk->flag_two_unique_paths)||(ptr_predec_blk->number_inputs_first_level_gate == 2)
				||(ptr_predec_blk_driver->number_input_addr_bits == 0)|| 
					((way_select)&&(ptr_dec->number_input_signals == 2))){ //this means that way_select is driving NAND2 in decoder. 
						ptr_predec_blk_driver->width_nand2_path_n[0] = minimum_width_nmos;
						ptr_predec_blk_driver->width_nand2_path_p[0] = 2 * ptr_predec_blk_driver->width_nand2_path_n[0];
						F = ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output / gatecap(ptr_predec_blk_driver->width_nand2_path_n[0] +
							ptr_predec_blk_driver->width_nand2_path_p[0], 0);
						ptr_predec_blk_driver->number_gates_nand2_path = (int) (log(F) / log(fopt) + 0.5);
						if(ptr_predec_blk_driver->number_gates_nand2_path%2 != 0){
							--ptr_predec_blk_driver->number_gates_nand2_path;
						}
						if(ptr_predec_blk_driver->number_gates_nand2_path < ptr_predec_blk_driver->min_number_gates){
							ptr_predec_blk_driver->number_gates_nand2_path = ptr_predec_blk_driver->min_number_gates;
						}
						//Recalculate the effective fanout of each stage for the above number_gates_last_level_predecoder
						f = pow(F, 1.0 / ptr_predec_blk_driver->number_gates_nand2_path);
						i = ptr_predec_blk_driver->number_gates_nand2_path - 1;
						ptr_predec_blk_driver->width_nand2_path_n[i] = ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output / (gatecap(1, 0) * f);
						ptr_predec_blk_driver->width_nand2_path_p[i]  = 2 * ptr_predec_blk_driver->width_nand2_path_n[i];
						for(i = ptr_predec_blk_driver->number_gates_nand2_path - 2; i >= 1; --i){
							ptr_predec_blk_driver->width_nand2_path_n[i] = ptr_predec_blk_driver->width_nand2_path_n[i+1] / f;
							ptr_predec_blk_driver->width_nand2_path_p[i] = 2 * ptr_predec_blk_driver->width_nand2_path_n[i];
						}
			}
			if((ptr_predec_blk->flag_two_unique_paths)||(ptr_predec_blk->number_inputs_first_level_gate == 3)||
					((way_select)&&(ptr_dec->number_input_signals == 3))){ //this means that way_select is driving NAND3 in decoder. 
						ptr_predec_blk_driver->width_nand3_path_n[0] = minimum_width_nmos;
						ptr_predec_blk_driver->width_nand3_path_p[0] = 2 * ptr_predec_blk_driver->width_nand3_path_n[0];
						F = ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output / gatecap(ptr_predec_blk_driver->width_nand3_path_n[0] +
							ptr_predec_blk_driver->width_nand3_path_p[0], 0);
						ptr_predec_blk_driver->number_gates_nand3_path = (int) (log(F) / log(fopt) + 0.5);
						if(ptr_predec_blk_driver->number_gates_nand3_path%2 != 0){
							--ptr_predec_blk_driver->number_gates_nand3_path;
						}
						if(ptr_predec_blk_driver->number_gates_nand3_path < ptr_predec_blk_driver->min_number_gates){
							ptr_predec_blk_driver->number_gates_nand3_path = ptr_predec_blk_driver->min_number_gates;
						}
						//Recalculate the effective fanout of each stage for the above number_gates_last_level_predecoder
						f = pow(F, 1.0 / ptr_predec_blk_driver->number_gates_nand3_path);
						i = ptr_predec_blk_driver->number_gates_nand3_path - 1;
						ptr_predec_blk_driver->width_nand3_path_n[i] =  ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output / (gatecap(1, 0) * f);
						ptr_predec_blk_driver->width_nand3_path_p[i]  = 2 * ptr_predec_blk_driver->width_nand3_path_n[i];
						for(i = ptr_predec_blk_driver->number_gates_nand3_path - 2; i >= 1; --i){
							ptr_predec_blk_driver->width_nand3_path_n[i] = ptr_predec_blk_driver->width_nand3_path_n[i+1] / f;
							ptr_predec_blk_driver->width_nand3_path_p[i] = 2 * ptr_predec_blk_driver->width_nand3_path_n[i];
						}
			}
		}
    }
}

void
delay_predecoder_block_driver(predecoder_block_driver *ptr_predec_blk_driver, double inrisetime_nand2_path,
        double inrisetime_nand3_path, double *outrisetime_nand2_path, double *outrisetime_nand3_path)
{
    int i;
    double rd, c_gate_load, c_load, c_intrinsic, tf, this_delay;
    *outrisetime_nand2_path = 0;
    *outrisetime_nand3_path = 0;
    if(ptr_predec_blk_driver->flag_driver_exists){
        for(i = 0; i < ptr_predec_blk_driver->number_gates_nand2_path - 1; ++i){
            rd = transreson(ptr_predec_blk_driver->width_nand2_path_n[i], NCH, 1);
            c_gate_load = gatecap(ptr_predec_blk_driver->width_nand2_path_p[i+1] + 
                    ptr_predec_blk_driver->width_nand2_path_n[i+1], 0.0);
            c_intrinsic = draincap(ptr_predec_blk_driver->width_nand2_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk_driver->width_nand2_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
            tf = rd * (c_intrinsic + c_gate_load);
            this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk_driver->delay_nand2_path += this_delay;
            inrisetime_nand2_path = this_delay / (1.0 - 0.5);
            ptr_predec_blk_driver->power_nand2_path.readOp.dynamic += 
                (c_gate_load + c_intrinsic) * 0.5 * vdd_periph_global * vdd_periph_global;
        }

        //Final inverter drives the predecoder block or the decoder output load 
        if(ptr_predec_blk_driver->number_gates_nand2_path != 0){
            i = ptr_predec_blk_driver->number_gates_nand2_path - 1;
            rd = transreson(ptr_predec_blk_driver->width_nand2_path_n[i], NCH, 1);
            c_intrinsic = draincap(ptr_predec_blk_driver->width_nand2_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk_driver->width_nand2_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
            c_load = ptr_predec_blk_driver->c_load_nand2_path_predecode_block_driver_output;
            tf = rd * (c_intrinsic + c_load) + 
                ptr_predec_blk_driver->r_load_nand2_path_predecode_block_driver_output * 
                c_load/ 2;
            this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk_driver->delay_nand2_path += this_delay;
            *outrisetime_nand2_path = this_delay / (1.0 - 0.5);	
            ptr_predec_blk_driver->power_nand2_path.readOp.dynamic += 
                (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
        }


        for(i = 0; i < ptr_predec_blk_driver->number_gates_nand3_path - 1; ++i){
            rd = transreson(ptr_predec_blk_driver->width_nand3_path_n[i], NCH, 1);
            c_gate_load = gatecap(ptr_predec_blk_driver->width_nand3_path_p[i+1] + 
                    ptr_predec_blk_driver->width_nand3_path_n[i+1], 0.0);
            c_intrinsic = draincap(ptr_predec_blk_driver->width_nand3_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk_driver->width_nand3_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
            tf = rd * (c_intrinsic + c_gate_load);
            this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk_driver->delay_nand3_path += this_delay;
            inrisetime_nand3_path = this_delay / (1.0 - 0.5);
            ptr_predec_blk_driver->power_nand3_path.readOp.dynamic += 
                (c_gate_load + c_intrinsic) * 0.5 * vdd_periph_global * vdd_periph_global;
        }

        //Final inverter drives the predecoder block or the decoder output load 
        if(ptr_predec_blk_driver->number_gates_nand3_path != 0){
            i = ptr_predec_blk_driver->number_gates_nand3_path - 1;
            rd = transreson(ptr_predec_blk_driver->width_nand3_path_n[i], NCH, 1);
            c_intrinsic = draincap(ptr_predec_blk_driver->width_nand3_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk_driver->width_nand3_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
            c_load = ptr_predec_blk_driver->c_load_nand3_path_predecode_block_driver_output;
            tf = rd * (c_intrinsic + c_load) + 
                ptr_predec_blk_driver->r_load_nand3_path_predecode_block_driver_output * 
                c_load / 2;
            this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk_driver->delay_nand3_path += this_delay;
            *outrisetime_nand3_path = this_delay / (1.0 - 0.5);	
            ptr_predec_blk_driver->power_nand3_path.readOp.dynamic += 
                (c_intrinsic + c_load) * 0.5 *	vdd_periph_global * vdd_periph_global;
        }
    }
}

void
compute_widths_predecoder_block(predecoder_block *ptr_predec_blk)
{
    int i;
    double F, f, c_load_nand3_path, c_load_nand2_path, c_input, c_load;
    if(ptr_predec_blk->flag_block_exists){ 
        ptr_predec_blk->branch_effort_nand2_gate_output = 1;
        ptr_predec_blk->branch_effort_nand3_gate_output = 1;
        ptr_predec_blk->number_inputs_first_level_gate = 0;

        if(ptr_predec_blk->number_input_addr_bits == 1){
            ptr_predec_blk->flag_two_unique_paths = 0;
            ptr_predec_blk->number_inputs_first_level_gate = 2;
            ptr_predec_blk->flag_second_level_gate = 0;
        } 
        if(ptr_predec_blk->number_input_addr_bits == 2){
            ptr_predec_blk->flag_two_unique_paths = 0;
            ptr_predec_blk->number_inputs_first_level_gate = 2;
            ptr_predec_blk->flag_second_level_gate = 0;
        } 
        if(ptr_predec_blk->number_input_addr_bits == 3){
            ptr_predec_blk->flag_two_unique_paths = 0;
            ptr_predec_blk->flag_second_level_gate = 0;		
            ptr_predec_blk->number_inputs_first_level_gate = 3;
        }
        if(ptr_predec_blk->number_input_addr_bits == 4){
            ptr_predec_blk->flag_two_unique_paths = 0;
            ptr_predec_blk->flag_second_level_gate = 2;	
            ptr_predec_blk->number_inputs_first_level_gate = 2;
            ptr_predec_blk->branch_effort_nand2_gate_output = 4;
        }
        if(ptr_predec_blk->number_input_addr_bits == 5){
            ptr_predec_blk->flag_two_unique_paths = 1;
            ptr_predec_blk->flag_second_level_gate = 2;	
            ptr_predec_blk->branch_effort_nand2_gate_output = 8;
			ptr_predec_blk->branch_effort_nand3_gate_output = 4;
        }
        if(ptr_predec_blk->number_input_addr_bits == 6){
            ptr_predec_blk->flag_two_unique_paths = 0;
            ptr_predec_blk->number_inputs_first_level_gate = 3;
            ptr_predec_blk->flag_second_level_gate = 2;	
            ptr_predec_blk->branch_effort_nand3_gate_output = 8;
        }
        if(ptr_predec_blk->number_input_addr_bits == 7){
            ptr_predec_blk->flag_two_unique_paths = 1;
            ptr_predec_blk->flag_second_level_gate = 3;	
            ptr_predec_blk->branch_effort_nand2_gate_output = 32;
            ptr_predec_blk->branch_effort_nand3_gate_output = 16;
        }
        if(ptr_predec_blk->number_input_addr_bits == 8){
            ptr_predec_blk->flag_two_unique_paths = 1;
            ptr_predec_blk->flag_second_level_gate = 3;	
            ptr_predec_blk->branch_effort_nand2_gate_output = 64;
           ptr_predec_blk->branch_effort_nand3_gate_output = 32;
		}           
        if(ptr_predec_blk->number_input_addr_bits == 9){
            ptr_predec_blk->flag_two_unique_paths = 0;
			ptr_predec_blk->number_inputs_first_level_gate = 3;
            ptr_predec_blk->flag_second_level_gate = 3;
            ptr_predec_blk->branch_effort_nand3_gate_output = 64;
        }

        //Find number of gates and sizing in second level of predecoder (if there is a second level)
		if(ptr_predec_blk->flag_second_level_gate){
            if(ptr_predec_blk->flag_second_level_gate == 2){//2nd level is a NAND2 gate
                ptr_predec_blk->width_second_level_n[0] = 2 * minimum_width_nmos;
                ptr_predec_blk->width_second_level_p[0] = 2 * minimum_width_nmos;
                F = gnand2 * ptr_predec_blk->c_load_predecoder_block_output / (gatecap(ptr_predec_blk->width_second_level_n[0], 0) +
                        gatecap(ptr_predec_blk->width_second_level_p[0], 0));
            }
            else{//2nd level is a NAND3 gate
                ptr_predec_blk->width_second_level_n[0] = 3 * minimum_width_nmos;
                ptr_predec_blk->width_second_level_p[0] = 3 * minimum_width_nmos;
                F = gnand3 * ptr_predec_blk->c_load_predecoder_block_output / (gatecap(ptr_predec_blk->width_second_level_n[0], 0) +
                        gatecap(ptr_predec_blk->width_second_level_p[0], 0));
            }
            ptr_predec_blk->number_gates_second_level = (int) (log(F) / log(fopt) + 0.5);
            //Check if number_gates_second_level_predecoder is odd. If it is, add 1 to make it even.
            if(ptr_predec_blk->number_gates_second_level%2 != 0) {
                --ptr_predec_blk->number_gates_second_level;
            }
            if(ptr_predec_blk->number_gates_second_level < ptr_predec_blk->min_number_gates_second_level){
                ptr_predec_blk->number_gates_second_level = ptr_predec_blk->min_number_gates_second_level;
            }
            //Recalculate the effective fanout of each stage for the above number_gates_last_level_predecoder
            f = pow(F, 1.0 / ptr_predec_blk->number_gates_second_level);
            i = ptr_predec_blk->number_gates_second_level - 1;
			c_input = ptr_predec_blk->c_load_predecoder_block_output / f;
            ptr_predec_blk->width_second_level_n[i]  = (1.0 / 3.0) * c_input / gatecap(1, 0);
            ptr_predec_blk->width_second_level_p[i] = 2 * ptr_predec_blk->width_second_level_n[i];

			if( ptr_predec_blk->width_second_level_n[i] > MAX_NMOS_WIDTH){
				 c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
				 if(ptr_predec_blk->flag_second_level_gate == 2){
					 F =  gnand2 * c_load / gatecap(ptr_predec_blk->width_second_level_n[0] +
						 ptr_predec_blk->width_second_level_p[0] , 0);
				 }
				 else{
					 if(ptr_predec_blk->flag_second_level_gate == 3){
						F =  gnand3 * c_load / gatecap(ptr_predec_blk->width_second_level_n[0] +
							ptr_predec_blk->width_second_level_p[0], 0);
					 }
				 }
				 ptr_predec_blk->number_gates_second_level = (int) (log(F) / log(fopt) + 0.5) + 1;
				 if(ptr_predec_blk->number_gates_second_level%2 != 0){
					 --ptr_predec_blk->number_gates_second_level;
				 }
				 if(ptr_predec_blk->number_gates_second_level < ptr_predec_blk->min_number_gates_second_level){
					 ptr_predec_blk->number_gates_second_level = ptr_predec_blk->min_number_gates_second_level;
				 }
				 f = pow(F, 1.0 / ptr_predec_blk->number_gates_second_level);
				 i = ptr_predec_blk->number_gates_second_level - 1;
				 ptr_predec_blk->width_second_level_n[i] = MAX_NMOS_WIDTH;
				 ptr_predec_blk->width_second_level_p[i] = 2 *  ptr_predec_blk->width_second_level_n[i];
			}

        for(i = ptr_predec_blk->number_gates_second_level - 2; i >= 1; --i){
			ptr_predec_blk->width_second_level_n[i] = ptr_predec_blk->width_second_level_n[i + 1] / f;
            ptr_predec_blk->width_second_level_p[i] = 2 * ptr_predec_blk->width_second_level_n[i];
		}


            //Now find number of gates and widths in first level of predecoder
		if((ptr_predec_blk->flag_two_unique_paths)||(ptr_predec_blk->number_inputs_first_level_gate == 2)){//Whenever flag_two_unique_paths is 
                //TRUE, it means first level of decoder employs both NAND2 and NAND3 gates. Or when 
                //number_inputs_first_level_gate is 2, it means a NAND2 gate is used in the first level of
                //the predecoder
                c_load_nand2_path = ptr_predec_blk->branch_effort_nand2_gate_output * 
                    (gatecap(ptr_predec_blk->width_second_level_n[0], 0) + 
                     gatecap(ptr_predec_blk->width_second_level_p[0], 0));
                ptr_predec_blk->width_first_level_nand2_path_n[0] = 2 * minimum_width_nmos;
                ptr_predec_blk->width_first_level_nand2_path_p[0] = 2 * minimum_width_nmos;
                F = gnand2 * c_load_nand2_path / (gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0], 0) +
                        gatecap(ptr_predec_blk->width_first_level_nand2_path_p[0], 0));
                ptr_predec_blk->number_gates_first_level_nand2_path = (int) (log(F) / log(fopt) + 0.5);
                //Check if number_gates_second_level_predecoder is odd. If it is, add 1 to make it even.
                if(ptr_predec_blk->number_gates_first_level_nand2_path % 2 != 0) {
                    --ptr_predec_blk->number_gates_first_level_nand2_path;
                }
                if(ptr_predec_blk->number_gates_first_level_nand2_path < ptr_predec_blk->min_number_gates_first_level){
                    ptr_predec_blk->number_gates_first_level_nand2_path = ptr_predec_blk->min_number_gates_first_level;
                }
                //Recalculate the effective fanout of each stage 
                f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand2_path);
                i = ptr_predec_blk->number_gates_first_level_nand2_path - 1;
				c_input = c_load_nand2_path / f;
                ptr_predec_blk->width_first_level_nand2_path_n[i]  = (1.0 / 3.0) * c_input / gatecap(1, 0);
                ptr_predec_blk->width_first_level_nand2_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand2_path_n[i];
					
				if( ptr_predec_blk->width_first_level_nand2_path_n[i] > MAX_NMOS_WIDTH){
					c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
					F =  gnand2 * c_load / gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0] , 0);
					ptr_predec_blk->number_gates_first_level_nand2_path = (int) (log(F) / log(fopt) + 0.5) + 1;
					if(ptr_predec_blk->number_gates_first_level_nand2_path%2 != 0){
						--ptr_predec_blk->number_gates_first_level_nand2_path;
				 }
					if(ptr_predec_blk->number_gates_first_level_nand2_path < ptr_predec_blk->min_number_gates_first_level){
						ptr_predec_blk->number_gates_first_level_nand2_path = ptr_predec_blk->min_number_gates_first_level;
				 }
					f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand2_path);
					i = ptr_predec_blk->number_gates_first_level_nand2_path - 1;
					ptr_predec_blk->width_first_level_nand2_path_n[i] = MAX_NMOS_WIDTH;
					ptr_predec_blk->width_first_level_nand2_path_p[i] = 2 *  ptr_predec_blk->width_first_level_nand2_path_n[i];
				}

                for(i = ptr_predec_blk->number_gates_first_level_nand2_path - 2; i >= 1; --i){
                    ptr_predec_blk->width_first_level_nand2_path_n[i] = ptr_predec_blk->width_first_level_nand2_path_n[i + 1] / f;
                    ptr_predec_blk->width_first_level_nand2_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand2_path_n[i];
                }
            }

            //Now find widths of gates along path in which first gate is a NAND3
            if((ptr_predec_blk->flag_two_unique_paths)||(ptr_predec_blk->number_inputs_first_level_gate == 3)){//Whenever flag_two_unique_paths is 
                //TRUE, it means first level of decoder employs both NAND2 and NAND3 gates. Or when 
                //number_inputs_first_level_gate is 3, it means a NAND3 gate is used in the first level of
                //the predecoder
                c_load_nand3_path = ptr_predec_blk->branch_effort_nand3_gate_output * 
                    (gatecap(ptr_predec_blk->width_second_level_n[0], 0) + 
                     gatecap(ptr_predec_blk->width_second_level_p[0], 0));
                ptr_predec_blk->width_first_level_nand3_path_n[0] = 3 * minimum_width_nmos;
                ptr_predec_blk->width_first_level_nand3_path_p[0] = 3 * minimum_width_nmos;
                F = gnand3 * c_load_nand3_path / (gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0], 0) +
                        gatecap(ptr_predec_blk->width_first_level_nand3_path_p[0], 0));
                ptr_predec_blk->number_gates_first_level_nand3_path = (int) (log(F) / log(fopt) + 0.5);
                //Check if number_gates_second_level_predecoder is odd. If it is, add 1 to make it even.
                if(ptr_predec_blk->number_gates_first_level_nand3_path % 2 != 0) {
                    --ptr_predec_blk->number_gates_first_level_nand3_path;
                }
                if(ptr_predec_blk->number_gates_first_level_nand3_path < ptr_predec_blk->min_number_gates_first_level){
                    ptr_predec_blk->number_gates_first_level_nand3_path = ptr_predec_blk->min_number_gates_first_level;
                }
                //Recalculate the effective fanout of each stage 
                f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand3_path);
                i = ptr_predec_blk->number_gates_first_level_nand3_path - 1;
				c_input = c_load_nand3_path / f;
                ptr_predec_blk->width_first_level_nand3_path_n[i]  = (1.0 / 3.0) * c_input / gatecap(1, 0);
                ptr_predec_blk->width_first_level_nand3_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand3_path_n[i];
				
				if( ptr_predec_blk->width_first_level_nand3_path_n[i] > MAX_NMOS_WIDTH){
					c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
					F =  gnand3 * c_load / gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0] , 0);
					ptr_predec_blk->number_gates_first_level_nand3_path = (int) (log(F) / log(fopt) + 0.5) + 1;
					if(ptr_predec_blk->number_gates_first_level_nand3_path%2 != 0){
						--ptr_predec_blk->number_gates_first_level_nand3_path;
				 }
					if(ptr_predec_blk->number_gates_first_level_nand3_path < ptr_predec_blk->min_number_gates_first_level){
						ptr_predec_blk->number_gates_first_level_nand3_path = ptr_predec_blk->min_number_gates_first_level;
				 }
					f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand3_path);
					i = ptr_predec_blk->number_gates_first_level_nand3_path - 1;
					ptr_predec_blk->width_first_level_nand3_path_n[i] = MAX_NMOS_WIDTH;
					ptr_predec_blk->width_first_level_nand3_path_p[i] = 2 *  ptr_predec_blk->width_first_level_nand3_path_n[i];
				}

               for(i = ptr_predec_blk->number_gates_first_level_nand3_path - 2; i >= 1; --i){
                    ptr_predec_blk->width_first_level_nand3_path_n[i] = ptr_predec_blk->width_first_level_nand3_path_n[i + 1] / f;
                    ptr_predec_blk->width_first_level_nand3_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand3_path_n[i];
                }
            }	
        }
		else{//Find number of gates and widths in first level of predecoder block when there is no second level 
            if(ptr_predec_blk->number_inputs_first_level_gate == 2){
                ptr_predec_blk->width_first_level_nand2_path_n[0] = 2 * minimum_width_nmos;
                ptr_predec_blk->width_first_level_nand2_path_p[0] = 2 * minimum_width_nmos;
                F = ptr_predec_blk->c_load_predecoder_block_output / (gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0], 0) +
                        gatecap(ptr_predec_blk->width_first_level_nand2_path_p[0], 0));
                ptr_predec_blk->number_gates_first_level_nand2_path = (int) (log(F) / log(fopt) + 0.5);
                //Check if number_gates_second_level_predecoder is odd. If it is, add 1 to make it even.
                if(ptr_predec_blk->number_gates_first_level_nand2_path % 2 != 0) {
                    --ptr_predec_blk->number_gates_first_level_nand2_path;
                }
                if(ptr_predec_blk->number_gates_first_level_nand2_path < ptr_predec_blk->min_number_gates_first_level){
                    ptr_predec_blk->number_gates_first_level_nand2_path = ptr_predec_blk->min_number_gates_first_level;
                }
                //Recalculate the effective fanout of each stage 
                f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand2_path);
                i = ptr_predec_blk->number_gates_first_level_nand2_path - 1;
                
				c_input = ptr_predec_blk->c_load_predecoder_block_output / f;
                ptr_predec_blk->width_first_level_nand2_path_n[i]  = (1.0 / 3.0) * c_input / gatecap(1, 0);
                ptr_predec_blk->width_first_level_nand2_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand2_path_n[i];
					
				if( ptr_predec_blk->width_first_level_nand2_path_n[i] > MAX_NMOS_WIDTH){
					c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
					F =  gnand2 * c_load / gatecap(ptr_predec_blk->width_first_level_nand2_path_n[0] +
						ptr_predec_blk->width_first_level_nand2_path_p[0] , 0);
					ptr_predec_blk->number_gates_first_level_nand2_path = (int) (log(F) / log(fopt) + 0.5) + 1;
					if(ptr_predec_blk->number_gates_first_level_nand2_path%2 != 0){
						--ptr_predec_blk->number_gates_first_level_nand2_path;
				 }
					if(ptr_predec_blk->number_gates_first_level_nand2_path < ptr_predec_blk->min_number_gates_first_level){
						ptr_predec_blk->number_gates_first_level_nand2_path = ptr_predec_blk->min_number_gates_first_level;
				 }
					f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand2_path);
					i = ptr_predec_blk->number_gates_first_level_nand2_path - 1;
					ptr_predec_blk->width_first_level_nand2_path_n[i] = MAX_NMOS_WIDTH;
					ptr_predec_blk->width_first_level_nand2_path_p[i] = 2 *  ptr_predec_blk->width_first_level_nand2_path_n[i];
				}

                for(i = ptr_predec_blk->number_gates_first_level_nand2_path - 2; i >= 1; --i){
                    ptr_predec_blk->width_first_level_nand2_path_n[i] = ptr_predec_blk->width_first_level_nand2_path_n[i + 1] / f;
                    ptr_predec_blk->width_first_level_nand2_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand2_path_n[i];
                }
            }

            if(ptr_predec_blk->number_inputs_first_level_gate == 3){
                ptr_predec_blk->width_first_level_nand3_path_n[0] = 3 * minimum_width_nmos;
                ptr_predec_blk->width_first_level_nand3_path_p[0] = 3 * minimum_width_nmos;
                F = ptr_predec_blk->c_load_predecoder_block_output / (gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0], 0) +
                        gatecap(ptr_predec_blk->width_first_level_nand3_path_p[0], 0));
                ptr_predec_blk->number_gates_first_level_nand3_path = (int) (log(F) / log(fopt) + 0.5);
                //Check if number_gates_second_level_predecoder is odd. If it is, add 1 to make it even.
                if(ptr_predec_blk->number_gates_first_level_nand3_path % 2 != 0) {
                    --ptr_predec_blk->number_gates_first_level_nand3_path;
                }
                if(ptr_predec_blk->number_gates_first_level_nand3_path < ptr_predec_blk->min_number_gates_first_level){
                    ptr_predec_blk->number_gates_first_level_nand3_path = ptr_predec_blk->min_number_gates_first_level;
                }
                //Recalculate the effective fanout of each stage 
                f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand3_path);
                i = ptr_predec_blk->number_gates_first_level_nand3_path - 1;
				c_input = ptr_predec_blk->c_load_predecoder_block_output / f;
                ptr_predec_blk->width_first_level_nand3_path_n[i]  = (1.0 / 3.0) * c_input / gatecap(1, 0);
                ptr_predec_blk->width_first_level_nand3_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand3_path_n[i];
				
				if( ptr_predec_blk->width_first_level_nand3_path_n[i] > MAX_NMOS_WIDTH){
					c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
					F =  gnand3 * c_load / gatecap(ptr_predec_blk->width_first_level_nand3_path_n[0] +
						ptr_predec_blk->width_first_level_nand3_path_p[0] , 0);
					ptr_predec_blk->number_gates_first_level_nand3_path = (int) (log(F) / log(fopt) + 0.5) + 1;
					if(ptr_predec_blk->number_gates_first_level_nand3_path%2 != 0){
						--ptr_predec_blk->number_gates_first_level_nand3_path;
				 }
					if(ptr_predec_blk->number_gates_first_level_nand3_path < ptr_predec_blk->min_number_gates_first_level){
						ptr_predec_blk->number_gates_first_level_nand3_path = ptr_predec_blk->min_number_gates_first_level;
				 }
					f = pow(F, 1.0 / ptr_predec_blk->number_gates_first_level_nand3_path);
					i = ptr_predec_blk->number_gates_first_level_nand3_path - 1;
					ptr_predec_blk->width_first_level_nand3_path_n[i] = MAX_NMOS_WIDTH;
					ptr_predec_blk->width_first_level_nand3_path_p[i] = 2 *  ptr_predec_blk->width_first_level_nand3_path_n[i];
				}

                for(i = ptr_predec_blk->number_gates_first_level_nand3_path - 2; i >= 1; --i){
                    ptr_predec_blk->width_first_level_nand3_path_n[i] = ptr_predec_blk->width_first_level_nand3_path_n[i + 1] / f;
                    ptr_predec_blk->width_first_level_nand3_path_p[i] = 2 * ptr_predec_blk->width_first_level_nand3_path_n[i];
                }
            }
        }
    }
}

void
delay_predecoder_block(predecoder_block *ptr_predec_blk, double inrisetime_nand2_path,
        double inrisetime_nand3_path, double *outrisetime_nand2_path,
        double *outrisetime_nand3_path)
{
    int i;
    double   rd, c_load, c_intrinsic, tf, this_delay;
    *outrisetime_nand2_path = 0;
    *outrisetime_nand3_path = 0;
    //First check whether a predecoder block is required
    if(ptr_predec_blk->flag_block_exists){
        //Find delay in first level of predecoder block
        //First find delay in path 
        if((ptr_predec_blk->flag_two_unique_paths)||(ptr_predec_blk->number_inputs_first_level_gate == 2)){
            //First gate is a NAND2 gate
            rd = transreson(ptr_predec_blk->width_first_level_nand2_path_n[0], NCH, 2);
            c_load = gatecap(ptr_predec_blk->width_first_level_nand2_path_n[1] +
                    ptr_predec_blk->width_first_level_nand2_path_p[1], 0.0);
            c_intrinsic = 2 * draincap(ptr_predec_blk->width_first_level_nand2_path_p[0], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk->width_first_level_nand2_path_n[0], NCH, 2, 1, DEFAULTHEIGHTCELL);
            tf = rd * (c_intrinsic + c_load);
            this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk->delay_nand2_path += this_delay;
            inrisetime_nand2_path = this_delay / (1.0 - 0.5);
            ptr_predec_blk->power_nand2_path.readOp.dynamic += 
                (c_load + c_intrinsic) * 0.5 * vdd_periph_global * vdd_periph_global;

            //Add delays of all but the last inverter in the chain
            for(i = 1; i < ptr_predec_blk->number_gates_first_level_nand2_path - 1; ++i){
                rd = transreson(ptr_predec_blk->width_first_level_nand2_path_n[i], NCH, 1);
                c_load = gatecap(ptr_predec_blk->width_first_level_nand2_path_n[i+1] +
                        ptr_predec_blk->width_first_level_nand2_path_p[i+1], 0.0);
                c_intrinsic = draincap(ptr_predec_blk->width_first_level_nand2_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_first_level_nand2_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand2_path += this_delay;
                inrisetime_nand2_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_nand2_path.readOp.dynamic += (c_intrinsic + c_load) * 0.5 *
                    vdd_periph_global * vdd_periph_global;
            }

            //Add delay of the last inverter
            i = ptr_predec_blk->number_gates_first_level_nand2_path - 1;
            rd = transreson(ptr_predec_blk->width_first_level_nand2_path_n[i], NCH, 1);
            if(ptr_predec_blk->flag_second_level_gate){
                c_load = ptr_predec_blk->branch_effort_nand2_gate_output * 
                    (gatecap(ptr_predec_blk->width_second_level_n[0], 0) + 
                     gatecap(ptr_predec_blk->width_second_level_p[0], 0));
                c_intrinsic = draincap(ptr_predec_blk->width_first_level_nand2_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_first_level_nand2_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand2_path += this_delay;
                inrisetime_nand2_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_nand2_path.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }
            else{//First level directly drives decoder output load
                c_load = ptr_predec_blk->c_load_predecoder_block_output;
                c_intrinsic = draincap(ptr_predec_blk->width_first_level_nand2_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_first_level_nand2_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load) + ptr_predec_blk->r_wire_predecoder_block_output *
                    c_load / 2; 
                this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand2_path += this_delay;
                *outrisetime_nand2_path = this_delay / (1.0 - 0.5);		
                ptr_predec_blk->power_nand2_path.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }
        }

        if((ptr_predec_blk->flag_two_unique_paths)||(ptr_predec_blk->number_inputs_first_level_gate == 3)){
            //Check if the number of gates in the first level is more than 1. 
            //First gate is a NAND3 gate
            rd = transreson(ptr_predec_blk->width_first_level_nand3_path_n[0], NCH, 3);
            c_load = gatecap(ptr_predec_blk->width_first_level_nand3_path_n[1] +
                    ptr_predec_blk->width_first_level_nand3_path_p[1], 0.0);
            c_intrinsic = 3 * draincap(ptr_predec_blk->width_first_level_nand3_path_p[0], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk->width_first_level_nand3_path_n[0], NCH, 3, 1, DEFAULTHEIGHTCELL);
            tf = rd * (c_intrinsic + c_load);
            this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk->delay_nand3_path += this_delay;
            inrisetime_nand3_path = this_delay / (1.0 - 0.5);
            ptr_predec_blk->power_nand3_path.readOp.dynamic += 
                (c_intrinsic + c_load) * 0.5 *	vdd_periph_global * vdd_periph_global;

            //Add delays of all but the last inverter in the chain
            for(i = 1; i < ptr_predec_blk->number_gates_first_level_nand3_path - 1; ++i){
                rd = transreson(ptr_predec_blk->width_first_level_nand3_path_n[i], NCH, 3);
                c_load = gatecap(ptr_predec_blk->width_first_level_nand3_path_n[i+1] +
                        ptr_predec_blk->width_first_level_nand3_path_p[i+1], 0.0);
                c_intrinsic = draincap(ptr_predec_blk->width_first_level_nand3_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_first_level_nand3_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand3_path += this_delay;
                inrisetime_nand3_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_nand3_path.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }

            //Add delay of the last inverter
            i = ptr_predec_blk->number_gates_first_level_nand3_path - 1;
            rd = transreson(ptr_predec_blk->width_first_level_nand3_path_n[i], NCH, 1);
            if(ptr_predec_blk->flag_second_level_gate){
                c_load = ptr_predec_blk->branch_effort_nand3_gate_output * 
                    (gatecap(ptr_predec_blk->width_second_level_n[0], 0) + 
                     gatecap(ptr_predec_blk->width_second_level_p[0], 0));
                c_intrinsic = draincap(ptr_predec_blk->width_first_level_nand3_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_first_level_nand3_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand3_path += this_delay;
                inrisetime_nand3_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_nand3_path.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 *	vdd_periph_global * vdd_periph_global;
            }
            else{//First level directly drives decoder output load
                c_load = ptr_predec_blk->c_load_predecoder_block_output;
                c_intrinsic = draincap(ptr_predec_blk->width_first_level_nand3_path_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_first_level_nand3_path_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load) + ptr_predec_blk->r_wire_predecoder_block_output *
                    c_load / 2; 
                this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand3_path += this_delay;
                *outrisetime_nand3_path = this_delay / (1.0 - 0.5);	
                ptr_predec_blk->power_nand3_path.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 *	vdd_periph_global * vdd_periph_global;
            }
        }	

        //Find delay through second level 
        if(ptr_predec_blk->flag_second_level_gate){
            if(ptr_predec_blk->flag_second_level_gate == 2){
                rd = transreson(ptr_predec_blk->width_second_level_n[0], NCH, 2);
                c_load = gatecap(ptr_predec_blk->width_second_level_n[1] +
                        ptr_predec_blk->width_second_level_p[1], 0.0);
                c_intrinsic = 2 * draincap(ptr_predec_blk->width_second_level_p[0], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_second_level_n[0], NCH, 2, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand2_path += this_delay;
                inrisetime_nand2_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_second_level.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }
            else{//flag_second_level_gate = 3)
                rd = transreson(ptr_predec_blk->width_second_level_n[0], NCH, 3);
                c_load = gatecap(ptr_predec_blk->width_second_level_n[1] +
                        ptr_predec_blk->width_second_level_p[1], 0.0);
                c_intrinsic = 3 * draincap(ptr_predec_blk->width_second_level_p[0], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_second_level_n[0], NCH, 3, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand3_path += this_delay;
                inrisetime_nand3_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_second_level.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }


            for(i = 1; i < ptr_predec_blk->number_gates_second_level - 1; ++i){
                rd = transreson(ptr_predec_blk->width_second_level_n[i], NCH, 1);
                c_load = gatecap(ptr_predec_blk->width_second_level_n[i+1] +
                        ptr_predec_blk->width_second_level_p[i+1], 0.0);
                c_intrinsic = draincap(ptr_predec_blk->width_second_level_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                    draincap(ptr_predec_blk->width_second_level_p[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand2_path += this_delay;
                inrisetime_nand2_path = this_delay / (1.0 - 0.5);
                this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
                ptr_predec_blk->delay_nand3_path += this_delay;
                inrisetime_nand3_path = this_delay / (1.0 - 0.5);
                ptr_predec_blk->power_second_level.readOp.dynamic += 
                    (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            }
            //Add delay of final inverter that drives the wordline decoders
            i = ptr_predec_blk->number_gates_second_level - 1;
            c_load = ptr_predec_blk->c_load_predecoder_block_output;
            rd = transreson(ptr_predec_blk->width_second_level_n[i], NCH, 1);
            c_intrinsic = draincap(ptr_predec_blk->width_second_level_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(ptr_predec_blk->width_second_level_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
            tf = rd * (c_intrinsic + c_load) + ptr_predec_blk->r_wire_predecoder_block_output *
                c_load / 2;
            this_delay = horowitz(inrisetime_nand2_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk->delay_nand2_path += this_delay;
            *outrisetime_nand2_path = this_delay / (1.0 - 0.5);	
            this_delay = horowitz(inrisetime_nand3_path, tf, 0.5, 0.5, RISE);
            ptr_predec_blk->delay_nand3_path += this_delay;
            *outrisetime_nand3_path = this_delay / (1.0 - 0.5);	
            ptr_predec_blk->power_second_level.readOp.dynamic += 
                (c_intrinsic + c_load) * 0.5 *	vdd_periph_global * vdd_periph_global;
        }
    }
}


void 
compute_widths_decoder(decoder *ptr_dec)
{
    int i;
    double F, f, c_input, c_load;
    if(ptr_dec->flag_decoder_exists){
        if(ptr_dec->number_input_signals == 2){
            ptr_dec->width_decoder_n[0] = 2 * minimum_width_nmos;
            ptr_dec->width_decoder_p[0] = 2 * minimum_width_nmos;
			F = gnand2 * ptr_dec->c_load_decoder_output / (gatecap(ptr_dec->width_decoder_n[0], 0) +
                    gatecap(ptr_dec->width_decoder_p[0], 0));
        }
        else{
            if(ptr_dec->number_input_signals == 3){
                ptr_dec->width_decoder_n[0] = 3 * minimum_width_nmos;
                ptr_dec->width_decoder_p[0] = 3 * minimum_width_nmos;
                F = gnand3 * ptr_dec->c_load_decoder_output / (gatecap(ptr_dec->width_decoder_n[0], 0) +
                        gatecap(ptr_dec->width_decoder_p[0], 0));
            }
        }
        ptr_dec->number_gates = (int) (log(F) / log(fopt) + 0.5);
        //Check if ptr_dec->number_gates  is odd. If it is, add 1 to make it even.
        if(ptr_dec->number_gates % 2 != 0){
            --ptr_dec->number_gates;
        }
        if(ptr_dec->number_gates < ptr_dec->min_number_gates){
            ptr_dec->number_gates = ptr_dec->min_number_gates;
        }

        //Recalculate the effective fanout of each stage 
        f = pow(F, 1.0 / ptr_dec->number_gates);
        i = ptr_dec->number_gates - 1;
        c_input = ptr_dec->c_load_decoder_output /  f;
		ptr_dec->width_decoder_n[i] = (1.0 / 3.0) * c_input / gatecap(1, 0);
		ptr_dec->width_decoder_p[i] = 2 * ptr_dec->width_decoder_n[i];

		if(ptr_dec->width_decoder_n[i] > MAX_NMOS_WIDTH){
			c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
			if(ptr_dec->number_input_signals == 2){
				F =  gnand2 * c_load / gatecap(ptr_dec->width_decoder_n[0] +
					ptr_dec->width_decoder_p[0] , 0);
			}
			else{
				if(ptr_dec->number_input_signals == 3){
					F =  gnand3 * c_load / gatecap(ptr_dec->width_decoder_n[0] +
						ptr_dec->width_decoder_p[0], 0);
				}
			}
			ptr_dec->number_gates = (int) (log(F) / log(fopt) + 0.5) + 1;
			if(ptr_dec->number_gates%2 != 0){
				--ptr_dec->number_gates;
			}
			if(ptr_dec->number_gates < ptr_dec->min_number_gates){
				ptr_dec->number_gates = ptr_dec->min_number_gates;
			}
			f = pow(F, 1.0 / ptr_dec->number_gates);
			i = ptr_dec->number_gates - 1;
			ptr_dec->width_decoder_n[i] = MAX_NMOS_WIDTH;
			ptr_dec->width_decoder_p[i] = 2 * ptr_dec->width_decoder_n[i];
		}
        
        for(i = ptr_dec->number_gates - 2; i >= 1; --i){
            ptr_dec->width_decoder_n[i] = ptr_dec->width_decoder_n[i+1] / f;
            ptr_dec->width_decoder_p[i] = 2 * ptr_dec->width_decoder_n[i];
        }
    }
}


void
delay_decoder(decoder *ptr_dec, double inrisetime, double *outrisetime)
{
    int i;
    double rd, tf, this_delay, c_load, c_intrinsic, vdd_wordline;
    *outrisetime = 0;
	if((is_wordline_transistor)&&(is_dram)){
		vdd_wordline = vpp;
	}
	else if(is_wordline_transistor){
		vdd_wordline = vdd_sram_cell;
	}
	else{
		vdd_wordline = vdd_periph_global;
	}
    if(ptr_dec->flag_decoder_exists){//First check whether a decoder is required at all
        if(ptr_dec->number_input_signals == 2){
            rd = transreson(ptr_dec->width_decoder_n[0], NCH, 2);
            c_load = gatecap(ptr_dec->width_decoder_n[1] +
                    ptr_dec->width_decoder_p[1], 0.0);
            c_intrinsic = 2 * draincap(ptr_dec->width_decoder_p[1], PCH, 1, 1, 4 * height_cell) + 
                draincap(ptr_dec->width_decoder_n[1], NCH, 2, 1, 4 * height_cell);
            tf = rd * (c_intrinsic + c_load);
            this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
            ptr_dec->delay += this_delay;
            inrisetime = this_delay / (1.0 - 0.5);
            ptr_dec->power.readOp.dynamic += (c_load + c_intrinsic) * 0.5 * vdd_wordline * vdd_wordline;
        }
        else{
            if(ptr_dec->number_input_signals == 3){
                rd = transreson(ptr_dec->width_decoder_n[0], NCH, 3);
                c_load = gatecap(ptr_dec->width_decoder_n[1] +
                        ptr_dec->width_decoder_p[1], 0.0);
                c_intrinsic = 2 * draincap(ptr_dec->width_decoder_p[0], PCH, 1, 1, 4 * height_cell) + 
                    draincap(ptr_dec->width_decoder_p[0], NCH, 2, 1, 4 * height_cell);
                tf = rd * (c_intrinsic + c_load);
                this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
                ptr_dec->delay += this_delay;
                inrisetime = this_delay / (1.0 - 0.5);	
                ptr_dec->power.readOp.dynamic += (c_load + c_intrinsic) * 0.5 * vdd_wordline *
                    vdd_wordline;
            }
        }
        for(i = 1; i < ptr_dec->number_gates - 1; ++i){
            rd = transreson(ptr_dec->width_decoder_n[i], NCH, 3);
            c_load = gatecap(ptr_dec->width_decoder_p[i+1] +
                    ptr_dec->width_decoder_n[i+1], 0.0);
            c_intrinsic = draincap(ptr_dec->width_decoder_p[i], PCH, 1, 1, 4 * height_cell) + 
                draincap(ptr_dec->width_decoder_n[i], NCH, 1, 1, 4 * height_cell);
            tf = rd * (c_intrinsic + c_load);
            this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
            ptr_dec->delay += this_delay;
            inrisetime = this_delay / (1.0 - 0.5);
            ptr_dec->power.readOp.dynamic += (c_load + c_intrinsic) * 0.5 * vdd_wordline * vdd_wordline;
        }
        //Add delay of final inverter that drives the wordline 
        i = ptr_dec->number_gates - 1;
        c_load = ptr_dec->c_load_decoder_output;
        rd = transreson(ptr_dec->width_decoder_n[i], NCH, 1);
        c_intrinsic = draincap(ptr_dec->width_decoder_p[i], PCH, 1, 1, 4 * height_cell) + 
            draincap(ptr_dec->width_decoder_n[i], NCH, 1, 1, 4 * height_cell);
        tf = rd * (c_intrinsic + c_load) + ptr_dec->r_wire_decoder_output * c_load / 2;
        this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
        ptr_dec->delay += this_delay;
        *outrisetime = this_delay / (1.0 - 0.5);
        ptr_dec->power.readOp.dynamic += (c_load + c_intrinsic) * 0.5 * vdd_wordline * vdd_wordline;
    }
}


void
delay_comparator(int tagbits, int A, double inputtime, double  *outputtime, 
        double *delay, powerDef *power)

{
    double Req, Ceq, tf, st1del, st2del, st3del, nextinputtime, m;
    double c1, c2, r1, r2, tstep, a, b, c, lkgCurrent;
    double Tcomparatorni;

    *delay = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;

    tagbits = tagbits / 4;//Assuming there are 4 quarter comparators. input tagbits is already
    //a multiple of 4.

    /* First Inverter */
    Ceq = gatecap(Wcompinvn2+Wcompinvp2, 0) +
        draincap(Wcompinvp1,PCH,1, 1, DEFAULTHEIGHTCELL) + draincap(Wcompinvn1,NCH,1, 1, DEFAULTHEIGHTCELL);
    Req = transreson(Wcompinvp1,PCH,1);
    tf = Req*Ceq;
    st1del = horowitz(inputtime,tf,VTHCOMPINV,VTHCOMPINV,FALL);
    nextinputtime = st1del/VTHCOMPINV;
    power->readOp.dynamic += 0.5 * Ceq * vdd_periph_global * vdd_periph_global * 4 * A;//For each degree of associativity 
    //there are 4 such quarter comparators
    lkgCurrent = 0.5 * cmos_ileakage(Wcompinvn1, Wcompinvp1, temper) * 4 * A;

    /* Second Inverter */
    Ceq = gatecap(Wcompinvn3+Wcompinvp3, 0) +
        draincap(Wcompinvp2,PCH,1, 1, DEFAULTHEIGHTCELL) + draincap(Wcompinvn2,NCH,1, 1, DEFAULTHEIGHTCELL);
    Req = transreson(Wcompinvn2,NCH,1);
    tf = Req*Ceq;
    st2del = horowitz(nextinputtime,tf,VTHCOMPINV,VTHCOMPINV,RISE);
    nextinputtime = st2del/(1.0-VTHCOMPINV);
    power->readOp.dynamic += 0.5 * Ceq * vdd_periph_global * vdd_periph_global * 4 * A;
    lkgCurrent += 0.5 * cmos_ileakage(Wcompinvn2, Wcompinvp2, temper) * 4 * A;

    /* Third Inverter */
    Ceq = gatecap(Wevalinvn+Wevalinvp, 0) +
        draincap(Wcompinvp3,PCH,1, 1, DEFAULTHEIGHTCELL) + draincap(Wcompinvn3,NCH,1, 1, DEFAULTHEIGHTCELL);
    Req = transreson(Wcompinvp3,PCH,1);
    tf = Req*Ceq;
    st3del = horowitz(nextinputtime,tf,VTHCOMPINV,VTHEVALINV,FALL);
    nextinputtime = st3del/(VTHEVALINV);
    power->readOp.dynamic += 0.5 * Ceq * vdd_periph_global * vdd_periph_global * 4 * A;
    lkgCurrent += 0.5 * cmos_ileakage(Wcompinvn3,Wcompinvp3,  temper) * 4 * A;

    /* Final Inverter (virtual ground driver) discharging compare part */
    r1 = transreson(Wcompn,NCH,2);
    r2 = transreson(Wevalinvn,NCH,1); /* was switch */
    c2 = (tagbits)*(draincap(Wcompn,NCH,1, 1, DEFAULTHEIGHTCELL)+draincap(Wcompn,NCH,2, 1, DEFAULTHEIGHTCELL))+
        draincap(Wevalinvp,PCH,1, 1, DEFAULTHEIGHTCELL) + draincap(Wevalinvn,NCH,1, 1, DEFAULTHEIGHTCELL);
    c1 = (tagbits)*(draincap(Wcompn,NCH,1, 1, DEFAULTHEIGHTCELL)+draincap(Wcompn,NCH,2, 1, DEFAULTHEIGHTCELL))
        +draincap(Wcompp,PCH,1, 1, DEFAULTHEIGHTCELL) + gatecap(WmuxdrvNANDn+WmuxdrvNANDp,0);
    power->readOp.dynamic += 0.5 * c2 * vdd_periph_global * vdd_periph_global * 4 * A;
    power->readOp.dynamic += c1 * vdd_periph_global * vdd_periph_global *  (A - 1);
    lkgCurrent += 0.5 * cmos_ileakage(Wevalinvn,Wevalinvp, temper) * 4 * A;
    lkgCurrent += 0.2 * 0.5 * cmos_ileakage(Wcompn,Wcompp, temper) * 4 * A;//stack factor of 0.2

    /* time to go to threshold of mux driver */
    tstep = (r2*c2+(r1+r2)*c1)*log(1.0/VTHMUXNAND);
    /* take into account non-zero input rise time */
    m = vdd_periph_global/nextinputtime;

    if((tstep) <= (0.5*(vdd_periph_global-v_th_periph_global)/m)) {
        a = m;
        b = 2*((vdd_periph_global*VTHEVALINV)-v_th_periph_global);
        c = -2*(tstep)*(vdd_periph_global-v_th_periph_global)+1/m*((vdd_periph_global*VTHEVALINV)-v_th_periph_global)*((vdd_periph_global*VTHEVALINV)-v_th_periph_global);
        Tcomparatorni = (-b+sqrt(b*b-4*a*c))/(2*a);
    }
    else{
        Tcomparatorni = (tstep) + (vdd_periph_global+v_th_periph_global)/(2*m) - (vdd_periph_global*VTHEVALINV)/m;
    }
    *outputtime = Tcomparatorni/(1.0-VTHMUXNAND);
    *delay = Tcomparatorni+st1del+st2del+st3del;
    power->readOp.leakage = lkgCurrent * vdd_periph_global;
}



void 
reset_powerDef(powerDef *power) {
    power->readOp.dynamic = 0.0;
    power->readOp.leakage = 0.0;

    power->writeOp.dynamic = 0.0;
    power->writeOp.leakage = 0.0;
}

void copy_powerDef(powerDef *dest, powerDef source) {
    dest->readOp.dynamic = source.readOp.dynamic;
    dest->readOp.leakage = source.readOp.leakage;

    dest->writeOp.dynamic = source.writeOp.dynamic;
    dest->writeOp.leakage = source.writeOp.leakage;
}

double objective_function(int flag_opt_for_delay, int flag_opt_for_dynamic_power, int flag_opt_for_leakage_power,
        int flag_opt_for_cycle_time, double delay_weight, double dynamic_power_weight, 
        double leakage_power_weight, double cycle_time_weight, double delay_wrt_min_delay, 
        double dyn_power_wrt_min_dyn_power, double leak_power_wrt_min_leak_power, 
        double cycle_time_wrt_min_cycle_time)
{
    double t = (double)(delay_wrt_min_delay * flag_opt_for_delay * delay_weight +
                 dyn_power_wrt_min_dyn_power * flag_opt_for_dynamic_power * dynamic_power_weight + 
                 leak_power_wrt_min_leak_power * flag_opt_for_leakage_power * leakage_power_weight + 
                 cycle_time_wrt_min_cycle_time * flag_opt_for_cycle_time * cycle_time_weight);
    return t;
}

/* 
 * Check whether the current organization
 * meets the power/delay/area/cycle time 
 * requirements
 */
int
check_org (uca_org_t *u, double  min_delay, double min_dyn,
        double min_leak, double min_cyc, double min_area)
{
    if (((u->access_time - min_delay)/min_delay)*100 > u->params->delay_dev) {
        return 0;
    }
    if (((u->power.readOp.dynamic - min_dyn)/min_dyn)*100 > 
                        u->params->dynamic_power_dev) {
        return 0;
    }
    if (((u->power.readOp.leakage - min_leak)/min_leak)*100 > 
                        u->params->leakage_power_dev) {
        return 0;
    }
    if (((u->cycle_time - min_cyc)/min_cyc)*100 > 
                         u->params->cycle_time_dev) {
        return 0;
    }
    if (((u->area - min_area)/min_area)*100 > u->params->area_dev) {
        return 0;
    }
    return 1;
}

int
check_nuca_org (nuca_org_t *n, double  min_delay, double min_dyn,
        double min_leak, double min_cyc, double min_area)
{
    if (((n->nuca_pda.delay - min_delay)*100/min_delay) > n->params->delay_dev_nuca) {
        return 0;
    }
    if (((n->nuca_pda.power.dynamic - min_dyn)/min_dyn)*100 > 
                        n->params->dynamic_power_dev_nuca) {
        return 0;
    }
    if (((n->nuca_pda.power.leakage - min_leak)/min_leak)*100 > 
                        n->params->leakage_power_dev_nuca) {
        return 0;
    }
    if (((n->nuca_pda.cycle_time - min_cyc)/min_cyc)*100 > 
                         n->params->cycle_time_dev_nuca) {
        return 0;
    }
    if (((n->nuca_pda.area_stats.area - min_area)/min_area)*100 > 
                         n->params->area_dev_nuca) {
        return 0;
    }
    return 1;
}

void
filter_arr (results_mem_array **a, double  min_delay, double min_dyn, 
        int s /* size */, int data)
{
    int i=0;
    int track = 0;
    int th1, th2;
	PRINTD(fprintf(stderr, "\nmin delay %g, dyn %g\n", min_delay, min_dyn));
    if (data) { // data array?
        th1 = th2 = 80;
    }
    else {
        th1 = th2 = 50;
    }
    while (i<s) {
        if ((((a[i]->access_time - min_delay)/min_delay)*100 > th1) &&
            (((a[i]->total_power.readOp.dynamic - min_dyn)/min_dyn)*100 > th2)) {
                free(a[i]);
                a[i] = NULL;
                track++;
        }
        i++;
    }
    PRINTD(fprintf(stderr, "filter_arr: Total no. of entries freed %d\n\n", track));
}

/* Objective function -- Version 6 */
uca_org_t *
find_optimal_uca (uca_res_lentry_t *ll, double min_delay, 
        double min_dyn, double min_leak, double min_cyc, 
        double min_area, input_params_t *parameters)
{
    double cost = 0;
    double min_cost = BIGNUM;
    uca_org_t *res = NULL;
    float d, a, dp, lp, c;
    int v;
    int te = 0;
    dp = parameters->dynamic_power_wt;
    lp = parameters->leakage_power_wt;
    a = parameters->area_wt;
    d = parameters->delay_wt;
    c = parameters->cycle_time_wt;
    
    while  (ll) {
        if (!(ll->fin_res.data_array)) break;
        te++;
        if (parameters->ed == 1) {
            cost = (ll->fin_res.access_time/min_delay)*
                   (ll->fin_res.power.readOp.dynamic/min_dyn);
            if (min_cost > cost) {
                min_cost = cost;
                res = &(ll->fin_res);
            }
        }
        else if (parameters->ed == 2) {
            cost = (ll->fin_res.access_time/min_delay)*
                    (ll->fin_res.access_time/min_delay)*
                    (ll->fin_res.power.readOp.dynamic/min_dyn);
            if (min_cost > cost) {
                min_cost = cost;
                res = &(ll->fin_res);
            }
        }
        else {

            /* 
             * check whether the current organization
             * meets the input deviation constraints
             */
            v = check_org(&(ll->fin_res), min_delay, min_dyn, 
                    min_leak, min_cyc, min_area);
            if (v) {
                cost = (d  * (ll->fin_res.access_time/min_delay) +
                        c  * (ll->fin_res.cycle_time/min_cyc) +
                        dp * (ll->fin_res.power.readOp.dynamic/min_dyn) +
                        lp * (ll->fin_res.power.readOp.leakage/min_leak) +
                        a  * (ll->fin_res.area/min_area));

                if (min_cost > cost) {
                    min_cost = cost;
                    res = &(ll->fin_res);
                }
            }
        }
        ll = ll->next_lentry;
    }
    return res;
}

nuca_org_t *
find_optimal_nuca (nuca_res_lentry_t *ll, double min_delay,
        double min_dyn, double min_leak, double min_cyc, 
        double min_area)
{
    double cost = 0;
    double min_cost = BIGNUM;
    nuca_org_t *res = NULL;
    float d, a, dp, lp, c;
    int v;
    input_params_t *parameters = ll->nres.params;
    dp = parameters->dynamic_power_wt_nuca;
    lp = parameters->leakage_power_wt_nuca;
    a = parameters->area_wt_nuca;
    d = parameters->delay_wt_nuca;
    c = parameters->cycle_time_wt_nuca;

//    fprintf(stderr, "min %g %g %g %g %g\n", min_delay, min_dyn,
//                    min_leak, min_cyc, min_area);
    
    while  (ll) {
//        PRINTD(fprintf(stderr, "\n-----------------------------"
//                       "---------------\n"));
//
//
//        fprintf(stderr, "NUCA___stats %d \tbankcount: lat = %g \tdynP = %g \twt = %d\t "
//                        "bank_dpower = %g \tleak = %g \tcycle = %g\n",
//                ll->nres.bank_count,
//                ll->nres.nuca_pda.delay, 
//                ll->nres.nuca_pda.power.dynamic,
//                ll->nres.h_wire.wt,
//                ll->nres.bank_pda.power.dynamic,
//                ll->nres.nuca_pda.power.leakage,
//                ll->nres.nuca_pda.cycle_time);
//        if (!(ll->nres.data_array)) break;


        if (parameters->ed == 1) {
            cost = (ll->nres.nuca_pda.delay/min_delay)*
                   (ll->nres.nuca_pda.power.dynamic/min_dyn);
            if (min_cost > cost) {
                min_cost = cost;
                res = &(ll->nres);
            }
        }
        else if (parameters->ed == 2) {
            cost = (ll->nres.nuca_pda.delay/min_delay)*
                   (ll->nres.nuca_pda.delay/min_delay)*
                   (ll->nres.nuca_pda.power.dynamic/min_dyn);
            if (min_cost > cost) {
                min_cost = cost;
                res = &(ll->nres);
            }
        }
        else {
            /* 
             * check whether the current organization
             * meets the input deviation constraints
             */
            v = check_nuca_org(&(ll->nres), min_delay, min_dyn, 
                    min_leak, min_cyc, min_area);
            if (min_leak == 0) min_leak = 0.1; //remove this after leakage modeling

            if (v) {
                cost = (d  * (ll->nres.nuca_pda.delay/min_delay) +
                        c  * (ll->nres.nuca_pda.cycle_time/min_cyc) +
                        dp * (ll->nres.nuca_pda.power.dynamic/min_dyn) +
                        lp * (ll->nres.nuca_pda.power.leakage/min_leak) +
                        a  * (ll->nres.nuca_pda.area_stats.area/min_area));
                //            fprintf(stderr, "cost = %g\n", cost);

                if (min_cost > cost) {
                    min_cost = cost;
                    res = &(ll->nres);
                }
            }
        }
        ll = ll->next_lentry;
    }
    return res;
}



double 
bitline_delay(int number_rows_subarray, int number_cols_subarray, int number_subarrays, 
					 double inrisetime, double *outrisetime,powerDef *power, double Tpre, 
					 int deg_bitline_muxing, int number_activated_mats_horizontal_direction, 
					 double *writeback_delay, int RWP, int ERP, int EWP)
{
 
    double Tbit, tau, v_bitline_precharge, v_th_mem_cell, v_wordline;
    double m, tstep;
	double dynRdEnergy = 0.0;
    double Icell, Iport;
    double Cbitrow_drain_capacitance, R_cell_pull_down, R_cell_acc, Cbitline, 
		Rbitline, C_drain_bit_mux, R_bit_mux, C_drain_sense_amp_iso, R_sense_amp_iso, 
		C_sense_amp_latch, C_drain_sense_amp_mux, r_dev;
	double leakage_power_cc_inverters_sram_cell, 
		leakage_power_access_transistors_read_write_or_write_only_port_sram_cell,
		leakage_power_read_only_port_sram_cell;

	*writeback_delay  = 0;
	leakage_power_cc_inverters_sram_cell = 0;
	leakage_power_access_transistors_read_write_or_write_only_port_sram_cell = 0;
	leakage_power_read_only_port_sram_cell = 0;

	if(is_dram){
		v_bitline_precharge = Vbitpre_dram;
		v_th_mem_cell = v_th_dram_access_transistor;
		v_wordline = vpp;
		is_access_transistor = 1;
		Cbitrow_drain_capacitance = draincap(Wmemcella_dram, NCH, 1, 0, width_cell) / 2.0;	/*due to shared contact*/ 
		//The access transistor is not folded. So we just need to specify a theshold value for the
		//folding width that is equal to or greater than Wmemcella. 
		R_cell_acc = transreson(Wmemcella_dram, NCH, 1);
		is_access_transistor = 0;
		Cbitline = number_rows_subarray / 2 * Cbitrow_drain_capacitance + number_rows_subarray * 
			Cbitmetal ;
		r_dev = vdd_dram_cell / I_on_dram_cell;
	}
	else{//SRAM
		v_bitline_precharge = Vbitpre_sram;
		v_th_mem_cell = v_th_sram_cell_transistor;
		v_wordline = vdd_sram_cell;
		 is_sram_cell = 1; 
		 Cbitrow_drain_capacitance = draincap(Wmemcella_sram, NCH, 1, 0, width_cell) / 2.0;	/* due to shared contact */
		 R_cell_pull_down = transreson(Wmemcellnmos_sram, NCH, 1);
		 R_cell_acc = transreson(Wmemcella_sram, NCH, 1);
		 Cbitline = number_rows_subarray * (Cbitrow_drain_capacitance + Cbitmetal);
	     //Leakage current of an SRAM cell
		 Iport = cmos_ileakage(Wmemcella_sram, 0,  temper); 
		 Icell = cmos_ileakage(Wmemcellnmos_sram, Wmemcellpmos_sram, temper);
		 is_sram_cell = 0; 
		 leakage_power_cc_inverters_sram_cell = Icell * vdd_sram_cell;
		 leakage_power_access_transistors_read_write_or_write_only_port_sram_cell = Iport * vdd_sram_cell;
		 leakage_power_read_only_port_sram_cell = 
			 leakage_power_access_transistors_read_write_or_write_only_port_sram_cell *
				 NAND2_LEAK_STACK_FACTOR;
	}

  
	Rbitline = number_rows_subarray * Rbitmetal;
	C_drain_bit_mux = draincap(width_nmos_bit_mux, NCH, 1, 0, width_cell / (2 *(RWP + ERP + RWP)));
	R_bit_mux = transreson(width_nmos_bit_mux, NCH, 1);
	C_drain_sense_amp_iso = draincap(Wiso, PCH, 1, 0, width_cell * deg_bitline_muxing / (RWP + ERP));
	R_sense_amp_iso = transreson(Wiso, PCH, 1);
	C_sense_amp_latch = gatecap(WsenseP + WsenseN, 0) + draincap(WsenseN, NCH, 1, 0, width_cell * deg_bitline_muxing / (RWP + ERP)) + 
		draincap(WsenseP, PCH, 1, 0, width_cell * deg_bitline_muxing / (RWP + ERP));
	C_drain_sense_amp_mux = draincap(width_nmos_sense_amp_mux, NCH, 1, 0, width_cell * deg_bitline_muxing / (RWP + ERP));
	if(is_dram){
		tstep = 2.3 * r_dev * (c_dram_cell * (Cbitline + 2 * C_drain_sense_amp_iso + 
			C_sense_amp_latch + C_drain_sense_amp_mux)) / (c_dram_cell + (Cbitline + 2 * 
			C_drain_sense_amp_iso + C_sense_amp_latch + C_drain_sense_amp_mux));
		*writeback_delay = tstep;
		dynRdEnergy +=  (Cbitline + 2 * C_drain_sense_amp_iso + 
			C_sense_amp_latch + C_drain_sense_amp_mux) * 
            vdd_dram_cell * vdd_dram_cell * number_cols_subarray * 4 *  number_activated_mats_horizontal_direction;
	}
	else{
		Vbitsense = (0.05 * vdd_sram_cell > VBITSENSEMIN) ? 0.05 * vdd_sram_cell : VBITSENSEMIN;
		if(deg_bitline_muxing > 1){
			tau = (R_cell_pull_down + R_cell_acc) * (Cbitline + 2 * C_drain_bit_mux + 
				2 * C_drain_sense_amp_iso + C_sense_amp_latch + C_drain_sense_amp_mux) + Rbitline * 
				(Cbitline / 2 + 2 * C_drain_bit_mux + 	2 * C_drain_sense_amp_iso + C_sense_amp_latch +
				C_drain_sense_amp_mux)+ R_bit_mux * (C_drain_bit_mux + 2 * C_drain_sense_amp_iso +
				C_sense_amp_latch + C_drain_sense_amp_mux) + R_sense_amp_iso * (C_drain_sense_amp_iso + C_sense_amp_latch + 
				C_drain_sense_amp_mux);
			dynRdEnergy += (Cbitline + 2 * C_drain_bit_mux) * Vbitsense * vdd_sram_cell * 
				number_cols_subarray * 4 * number_activated_mats_horizontal_direction;
			dynRdEnergy += (C_drain_sense_amp_iso +	C_sense_amp_latch +	C_drain_sense_amp_mux) * 
				Vbitsense * vdd_sram_cell * (number_cols_subarray * 4 / deg_bitline_muxing) *
				number_activated_mats_horizontal_direction;
		}
		else{//deg_bitline_muxing == 1
			tau = (R_cell_pull_down + R_cell_acc) * (Cbitline + C_drain_sense_amp_iso +
				C_sense_amp_latch + C_drain_sense_amp_mux) + Rbitline * Cbitline / 2 + 
				R_sense_amp_iso * (C_drain_sense_amp_iso + C_sense_amp_latch + 	C_drain_sense_amp_mux);
			dynRdEnergy += (Cbitline + C_drain_sense_amp_iso + C_sense_amp_latch + 
				C_drain_sense_amp_mux) * Vbitsense * vdd_sram_cell * number_cols_subarray * 4 * 
				number_activated_mats_horizontal_direction;
		}
		tstep = tau * log(v_bitline_precharge / (v_bitline_precharge - Vbitsense));
		power->readOp.leakage = (Icell + Iport) * vdd_sram_cell;
		power->readOp.leakage = leakage_power_cc_inverters_sram_cell + 
			leakage_power_access_transistors_read_write_or_write_only_port_sram_cell + 
			leakage_power_access_transistors_read_write_or_write_only_port_sram_cell * 
			(RWP + EWP -1) + leakage_power_read_only_port_sram_cell * ERP;
	}
  
    /* take input rise time into account */
    m = v_wordline / inrisetime;
    if (tstep <= (0.5 * (v_wordline - v_th_mem_cell) / m))
	{
		Tbit = sqrt(2 * tstep * (v_wordline - v_th_mem_cell)/ m);
    }
    else
    {
        Tbit = tstep + (v_wordline - v_th_mem_cell) / (2 * m);
    }

	power->readOp.dynamic = dynRdEnergy;
    *outrisetime = 0;
	return(Tbit);
}


double
sense_amp_input_cap()
{
    return draincap(Wiso, PCH, 1, 1, DEFAULTHEIGHTCELL) +  
        gatecap(WsenseP + WsenseN, 0) + 
        draincap(WsenseN, NCH, 1, 1, DEFAULTHEIGHTCELL) +  
        draincap(WsenseP, PCH, 1, 1, DEFAULTHEIGHTCELL);
}


/* Time taken by the pull-down transistor in the memory
   cell to discharge the bitline by Vbitsense. 
   Employs distributed RC wire model */
double 
bitline_delay_sram(int number_rows_subarray, 
    int number_cols_subarray, int number_subarrays, 
	double inrisetime, double *outrisetime,powerDef *power, 
    double Tpre, int deg_bitline_muxing, 
    int number_activated_mats_horizontal_direction, 
	double *writeback_delay, int RWP, int ERP, int EWP)
{
    double Rpd, Cpd, /* pulldown trans. res and cap */
           Rpass, Cpass, 
           r, c, /* res and cap /len of wire */
           Rbmux, Cbmux, /* bit line mux */
           Csmux, /* sense amp mux */
           Csensein;

    double tau, delay, m, Tbit;

    double Vth = v_th_sram_cell_transistor;
    double Vdd = vdd_sram_cell;
    /* length of bitline in terms of cell count */
    double len = number_rows_subarray;

    /* for leakage calc. */
    double Iport, Icell, cell_leakage,
        rw_w_port_leakge, r_port_leakage;

    *writeback_delay = 0;
    Vbitsense = (0.05 * vdd_sram_cell > VBITSENSEMIN) ? 0.05 * vdd_sram_cell : VBITSENSEMIN;

    Rpd = transreson (minimum_width_nmos, NCH, 1);
    Cpd = draincap(minimum_width_nmos, NCH, 1, 1, width_cell)/2; 
                            /*shared contact in the layout */
                            
    Rpass = transreson(minimum_width_nmos, NCH, 1);
    Cpass = draincap(minimum_width_nmos, NCH, 1, 1, 
                              DEFAULTHEIGHTCELL)/2; 
                  /*shared contact in the layout */

    r = wire_res (height_cell*1e-6 /* in m */, 1, 1);
    c = Cbitmetal + Cpass;

    if (deg_bitline_muxing > 1) {
        Rbmux = transreson(minimum_width_nmos, NCH, 1);
        Cbmux = draincap(minimum_width_nmos, NCH, 1, 1, 
                                DEFAULTHEIGHTCELL);
    }
    else {
        Rbmux = 0;
        Cbmux = 0;
    }
    Csmux = Cbmux;

    Csensein = sense_amp_input_cap();

    /* time const */
    tau = (Rpass + Rpd)*Cpass + 
            (Rpd + Rpass + r*len + Rbmux)*Cbmux +
            (Rpd + Rpass)*(c*len) + len*(len+1)*r*c/2;
    
    delay = -tau*log((Vbitpre_sram - Vbitsense)/Vbitpre_sram);

    power->readOp.dynamic = (Vbitsense*vdd_sram_cell*
                                  (Cpd + Cpass + c*len + Csensein)) *
                            number_cols_subarray * 4 * 
                            number_activated_mats_horizontal_direction;

    /* Leakage power of an SRAM cell
       NOTE: Total leakage value is calculated in calculate_time */
    Iport = cmos_ileakage(Wmemcella_sram, 0,  temper); 
    Icell = cmos_ileakage(Wmemcellnmos_sram, 
                           Wmemcellpmos_sram, temper);
    
    /* cross coupled inv. */
    cell_leakage = 2*Icell * vdd_sram_cell;
    /* rw port or w port */
    rw_w_port_leakge = Iport * vdd_sram_cell;
    /* read port */
    r_port_leakage = rw_w_port_leakge * NAND2_LEAK_STACK_FACTOR;

    power->readOp.leakage = cell_leakage + 
        rw_w_port_leakge + 
        rw_w_port_leakge * 
        (RWP + EWP -1) + r_port_leakage * ERP;

    /* if the wordline signal has a finite risetime */
    if (inrisetime != 0) {
        m = Vdd/inrisetime;
        if (delay <= (0.5*Vdd-Vth)/m) {
            Tbit = sqrt(2*delay*(Vdd-Vth)/m);
        }
        else {
            Tbit = delay + (Vdd - Vth)/(2*m);
        }
    }
    else 
        Tbit = delay;

    *outrisetime = 0;

    return Tbit;
}


void delay_sense_amplifier(int number_cols_subarray, int RWP, int ERP, double inrisetime,double *outrisetime, powerDef *power, 
        dataout_htree_node *ptr_htree_node, int deg_bitline_muxing, int number_mats,
        int number_activated_mats_horizontal_direction,  double *delay)
{
    double c_load;
    int  number_sense_amps_subarray;
    double IsenseEn, IsenseN, IsenseP, Iiso;
    double lkgIdlePh, lkgReadPh, lkgWritePh;
    double lkgRead, lkgIdle;

    *delay = 0;
    *outrisetime = 0;
    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;
    //v5.0
//    double  tau;
    number_sense_amps_subarray = number_cols_subarray / deg_bitline_muxing;//in a subarray

    //Sense-amp leakage 
    Iiso = simplified_pmos_leakage(Wiso, temper);
    IsenseEn = simplified_nmos_leakage(WsenseEn, temper);
    IsenseN  = simplified_nmos_leakage(WsenseN, temper);
    IsenseP  = simplified_pmos_leakage(WsenseP, temper);


    lkgIdlePh = IsenseEn;
    lkgWritePh = Iiso + IsenseEn;
    lkgReadPh = Iiso + IsenseN + IsenseP;
    lkgRead = lkgReadPh * number_sense_amps_subarray * 4 * number_activated_mats_horizontal_direction + 
		lkgIdlePh * number_sense_amps_subarray * 4 *
		(number_mats - number_activated_mats_horizontal_direction);
	lkgIdle = lkgIdlePh * number_sense_amps_subarray * 4 * number_mats;

    power->readOp.leakage = lkgIdle * vdd_periph_global;
    /* using SPICE generated values for senseamp dynamic power and delay */
    *delay = SENSE_AMP_D;
    power->readOp.dynamic = SENSE_AMP_P * number_sense_amps_subarray * 4;
}

#define VTHMODEL 1 /* FIXME: as of now it is for high performance */
#define ITRS_GLOBAL 1
#define ASPECT_RATIO 2.2 /* wire height/width */
double EPS0 = 8.8541878176e-12; /* in SI unit */
double EPSR = 1.5;

/* Low swing wire model */

/*
 * The fall time of an input signal to the first stage of a circuit is
 * assumed to be same as the fall time of the output signal of two 
 * inverters connected in series (refer: CACTI 1 Technical report, 
 * section 6.1.3) 
 */
double
signal_fall_time () 
{

    /* rise time of inverter 1's output */
    double rt; 
    /* fall time of inverter 2's output */
    double ft; 
    double timeconst; /* rc time constant */

    timeconst = (draincap(minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL) +
         draincap(minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) +
         gatecap(minimum_width_pmos + minimum_width_nmos, 0)) *
         transreson(minimum_width_pmos, PCH, 1);
    rt = horowitz (0, timeconst, v_th[VTHMODEL], v_th[VTHMODEL], FALL) / (1 - v_th[VTHMODEL]);
    timeconst = (draincap(minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL) +
         draincap(minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) +
         gatecap(minimum_width_pmos + minimum_width_nmos, 0)) *
         transreson(minimum_width_nmos, NCH, 1);
    ft = horowitz (rt, timeconst, v_th[0], v_th[0], RISE) / v_th[0];
    return ft;
}

double
signal_rise_time () 
{

    /* rise time of inverter 1's output */
    double ft; 
    /* fall time of inverter 2's output */
    double rt; 
    double timeconst; /* rc time constant */

    timeconst = (draincap(minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL) +
         draincap(minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) +
         gatecap(minimum_width_pmos + minimum_width_nmos, 0)) *
         transreson(minimum_width_nmos, NCH, 1);
    rt = horowitz (0, timeconst, v_th[VTHMODEL], v_th[VTHMODEL], RISE) / v_th[VTHMODEL];
    timeconst = (draincap(minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL) +
         draincap(minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) +
         gatecap(minimum_width_pmos + minimum_width_nmos, 0)) *
         transreson(minimum_width_pmos, PCH, 1);
    ft = horowitz (rt, timeconst, v_th[VTHMODEL], v_th[VTHMODEL], FALL) / (1 - v_th[VTHMODEL]);
    return ft; //sec
}

/* Wire resistance and capacitance calculations
 *   wire width
 *
 *    /__/
 *   |  |
 *   |  |  height = ASPECT_RATIO*wire width (ASPECT_RATIO = 2.2, ref: ITRS)
 *   |__|/
 *
 *   spacing between wires in same level = wire width
 *   spacing between wires in adjacent levels = wire width
 */

double
wire_cap (double len /* in m */,
          /* width and spacing scaling factor can be used
           * to model low level wires or special
           * fat wires
           */
          double w_scaling_factor, 
          double s_scaling_factor)
{
    double sidewall, adj, tot_cap;
    double wire_height;
    double wire_width, wire_spacing;
    double M = 1.4; 
    wire_width = wire_pitch[VTHMODEL][ITRS_GLOBAL];
    wire_width *= 0.5e-6; /* width is half the pitch (in m) */
    wire_height = wire_width*ASPECT_RATIO*1.2;
//    wire_height = wire_width*3;
    wire_spacing = wire_width;

    /* for non-standard wires */
    wire_width *= w_scaling_factor;
    wire_spacing *= s_scaling_factor;


    /* capacitance between wires in the same level */
    sidewall = M * EPS0 * EPSR * (wire_height/wire_spacing);

    /* capacitance between wires in adjacent levels */
    adj = M * EPS0 * EPSR * (w_scaling_factor);

    tot_cap =  (sidewall + adj + (c_fringe[VTHMODEL] * 1e6)); //F/m

    PRINTDW(fprintf(stderr, "wire_cap:len %f, adj %g, side %g, tot %g, w_scaling %f,"
                         " s_scaling %f\n", len, adj, sidewall, tot_cap*len,
                         w_scaling_factor, s_scaling_factor));
    
    return (tot_cap*len); // (F)
}

double
wire_res (double len /*(in m)*/, 
          /* width and spacing scaling factor can be used
           * to model low level wires or special
           * fat wires
           */
          double w_scaling_factor,
          double s_scaling_factor)
{
    double wire_width;
#ifdef SEMI_GLOBAL_WIRES
    wire_width = semi_global_interconnect_pitch[VTHMODEL];
#else
    wire_width = wire_pitch[VTHMODEL][ITRS_GLOBAL];
#endif

#define dish .9
    wire_width *= 0.5e-6; /* width is half the pitch (in m) */

    PRINTDW(fprintf(stderr, "wire_res: len %f, w_scaling %f,"
                    " s_scaling %f, res %g\n", len, w_scaling_factor,
                    s_scaling_factor, 
                    (dish*Cu_resistivity * len/(ASPECT_RATIO*wire_width*
                    wire_width*w_scaling_factor))));

    return (dish * Cu_resistivity * len/(ASPECT_RATIO*wire_width*
            wire_width*w_scaling_factor)); 
}

/*
 * Calculates the delay, power and area of the transmitter circuit.
 *
 * The transmitter delay is the sum of nand gate delay, inverter delay
 * low swing nmos delay, and the wire delay
 * (ref: Technical report 6)
 */
void
low_swing_model(wire_stats_t *wstat) 
{
    double delay;
    double power;
    double req_cin, st_eff, inv_size;
    double inputrise, timeconst, res_eq, cap_eq;
    double driver_res;
    double wire_width, cwire, rwire;
    double nsize;
    double len = wstat->wire_length;
    double rt, sense_delay;
    double c_load, tau;
    double sense_power;

    PRINTDW(fprintf(stderr, "low_swing_model: length of the wire "
                    " %g\n", len));
    
    inputrise = signal_rise_time();
#define HPerf 0

    /* Final nmos low swing driver size calculation:
     * Try to size the driver such that the delay
     * is less than 8FO4. 
     * If the driver size is greater than
     * the max allowable size, assume max size for the driver.
     * In either case, recalculate the delay using
     * the final driver size assuming slow input with
     * finite rise time instead of ideal step input
     * 
     * (ref: Technical report 6)
     */
    cwire = wire_cap(len, 1 , 
                      1); /* load capacitance */
    rwire = wire_res(len, 1, 1);
#define RES_ADJ (8.6) // Increase in resistance due to low driving vol.
    driver_res = (-8*FO4*1e-12/(log(0.5) * cwire))/RES_ADJ; 
    nsize = restowidth(driver_res, NCH);
    if (nsize > 100 * FEATURESIZE) {
        PRINTDW(fprintf(stderr, "transmitter_delay: max nsize "
            "@ len %f\n", len));
        nsize = 100 * FEATURESIZE;
    }
    else if (nsize < 3 * FEATURESIZE/2) {
        nsize = minimum_width_nmos;
    }
    if(rwire*cwire > 8*FO4*1e-12) {
        nsize = 100 * FEATURESIZE;
    }

    /* size the inverter appropriately to minimize the transmitter delay */
    st_eff = sqrt((4/3)*gatecap(nsize, 0)/(gatecap(2*minimum_width_nmos, 0)
            + gatecap(2*minimum_width_pmos, 0)));
    req_cin = ((4/3)*gatecap(nsize, 0))/st_eff;
    inv_size = req_cin/(gatecap(minimum_width_pmos, 0) +
                               gatecap(minimum_width_nmos, 0));

    if(inv_size < 1) inv_size = 1;

    /* nand gate delay */
    res_eq = (2 * transreson(minimum_width_nmos, NCH, 1));
    cap_eq = 2 * draincap(minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) + 
             draincap(2*minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL) +
             gatecap(inv_size*minimum_width_nmos, 0) +
             gatecap(inv_size*minimum_width_pmos, 0);

    timeconst = res_eq * cap_eq;
    
    delay = horowitz(inputrise, timeconst, v_th[HPerf]/vdd_periph_global,  
            v_th[HPerf]/vdd_periph_global, RISE);
    power = cap_eq*vdd_periph_global*vdd_periph_global;

    PRINTDW(fprintf(stderr, "transmitter_delay: nand gate delay, "
                    "power %g, %g\n", delay, power*2));
    inputrise = delay / (1 - v_th[HPerf]); /* for the next stage */

    /* Inverter delay:
     * The load capacitance of this inv depends on
     * the gate capacitace of the final stage nmos
     * transistor which in turn depends on nsize
     */
    res_eq = transreson(inv_size*minimum_width_pmos, PCH, 1);
    cap_eq = draincap(inv_size*minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) + 
             draincap(inv_size*minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL) +
             gatecap(nsize, 0);
    timeconst = res_eq * cap_eq;

    delay += horowitz(inputrise, timeconst, v_th[HPerf]/vdd_periph_global,
             v_th[HPerf]/vdd_periph_global, FALL);
    power += cap_eq*vdd_periph_global*vdd_periph_global;

    PRINTDW(fprintf(stderr, "transmitter_delay: nand gate + inv delay, "
                    "power %g, %g\n", delay, power*2));

    wstat->transmitter.delay = delay;
    wstat->transmitter.power.dynamic = power*2; /* since it is a diff. model*/
    wstat->transmitter.power.leakage = 0.5 * vdd_periph_global * 
	    (4 * cmos_ileakage(minimum_width_nmos, minimum_width_pmos, temper) * 
				NAND2_LEAK_STACK_FACTOR +
        4 * cmos_ileakage(minimum_width_nmos, minimum_width_nmos, temper));
    
    inputrise = delay / v_th[HPerf];

    /* nmos delay + wire delay */
    wire_width = wire_pitch[HPerf][ITRS_GLOBAL];
    wire_width *= 1e-6; /* low swing wire width = 2* global 
                         * wire width 
						 * and spacing = width. 
						 */
    wstat->wire_width = wire_width;
    wstat->wire_spacing = wire_width/2;
    cap_eq = cwire + draincap(nsize, NCH, 1, 1, DEFAULTHEIGHTCELL)*2 + 
        wstat->nsense * sense_amp_input_cap(); //+receiver cap
    /* 
     * NOTE: nmos is used as both pull up and pull down transistor
     * in the transmitter. This is because for low voltage swing, drive 
     * resistance of nmos is less than pmos 
     * (for a detailed graph ref: On-Chip Wires: Scaling and Efficiency)
     */
    timeconst = (transreson(nsize, NCH, 1)*RES_ADJ) * (cwire + 
                draincap(nsize, NCH, 1, 1, DEFAULTHEIGHTCELL)*2) +
                rwire*cwire/2 + 
                (transreson(nsize, NCH, 1) + rwire) * 
                wstat->nsense * sense_amp_input_cap();

    /* 
     * since we are pre-equalizing and overdriving the low 
     * swing wires, the net time constant is less 
     * than the actual value
     */
    delay += horowitz(inputrise, timeconst, v_th[HPerf]/vdd_periph_global, 
             .25, 0);
    power += cap_eq*VOL_SWING*.400; /* .4v is the over drive voltage */
    power *= 2; /* differential wire */

    wstat->l_wire.delay = delay - wstat->transmitter.delay;
    wstat->l_wire.power.dynamic = power - wstat->transmitter.power.dynamic;
    wstat->l_wire.power.leakage = 0.5 * vdd_periph_global * 
        (4* simplified_nmos_leakage (nsize, temper));

    PRINTDW(fprintf(stderr, "transmitter_delay: nand gate + inv "
                    "+ nmos + wire delay, power %g, %g\n", delay, power));

    rt = horowitz(inputrise, timeconst, v_th[HPerf]/vdd_periph_global,  
                    v_th[HPerf]/vdd_periph_global, RISE)/v_th[HPerf];

    /* sense amp delay */
    c_load = gatecap(WsenseP + WsenseN, 0) + 
             draincap(WsenseN, NCH, 1, 1, DEFAULTHEIGHTCELL) + 
             draincap(WsenseP, PCH, 1, 1, DEFAULTHEIGHTCELL) +
		     draincap(Wiso, PCH, 1, 1, DEFAULTHEIGHTCELL);
    tau = c_load / Gm_sense_amp_latch;
    sense_delay = SENSE_AMP_D; // in s
    sense_power = SENSE_AMP_P; // J 
    delay += sense_delay;

    wstat->sense_amp.delay = sense_delay;
    wstat->sense_amp.power.dynamic = sense_power;
    wstat->sense_amp.power.leakage = 0; //FIXME

    PRINTDW(fprintf(stderr, "sense_amp: power %g delay %g, "
                "input cap %g \n", sense_delay, sense_power, 
                sense_amp_input_cap()));

    wstat->wire_pda.delay = delay;
    wstat->wire_pda.power.dynamic = power;
    wstat->wire_pda.power.leakage = wstat->transmitter.power.leakage +
                                     wstat->l_wire.power.leakage +
                                     wstat->sense_amp.power.leakage;
}

void
wire_energy2 (wire_stats_t *w)
{
    double len = w->wire_length;
    pda_res_t *wire_stats = &(w->wire_pda);
     /* width is half the pitch (in m) */
    double min_wire_width = wire_pitch[VTHMODEL][ITRS_GLOBAL] * .5e-6;
    double w_scaling_factor = w->wire_width/min_wire_width;
    double s_scaling_factor = w->wire_spacing/min_wire_width;
    /* switching energy */
    double switching = 0;
    /* short-circuit energy */
    double short_ckt = 0;
    /* time constant */
    double tc = 0;
#define CO_CAP 1.25e-15
#define CP_CAP 1.25e-15
    /* input cap of min sized driver */
    double input_cap = gatecap (minimum_width_nmos + 
                                minimum_width_pmos, 0);
    //double input_cap = CO_CAP;
    /* 
     * output parasitic capacitance of
     * the min. sized driver
     */
    double out_cap = draincap(minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) +
            draincap(minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL);
    //double out_cap = CP_CAP;
    /* drive resistance */
    double out_res = (transreson(minimum_width_nmos, NCH, 1) +
                transreson(minimum_width_pmos, PCH, 1))/2;
    /* wire res */
    double wr = wire_res(len, w_scaling_factor, 
                    s_scaling_factor);
    /* wire cap /m */
    double wc = wire_cap(len, w_scaling_factor,
                    s_scaling_factor);
    /* 
     * size the repeater such that the delay of the wire
     * is minimum
     */
    double repeater_scaling = sqrt(out_res*wc/(wr*input_cap)); // len will cancel
    /* 
     * calc the optimum spacing between the repeaters (m)
     */
    double repeater_spacing = sqrt(2 * out_res * (out_cap + input_cap)/
                                    ((wr/len)*(wc/len)));
    w->repeater_size = repeater_scaling;
    w->repeater_spacing = repeater_spacing;
    PRINTDW(fprintf(stderr, "wire_energy: repeater spacing (m) %g, "
                                          "repeater sizing %g\n",
                            repeater_spacing, repeater_scaling));

    switching = (repeater_scaling * (input_cap + out_cap) +
         repeater_spacing * (wc/len)) * vdd_periph_global * vdd_periph_global;

    tc = out_res * (input_cap + out_cap) + 
         out_res * (wc/len) * repeater_spacing/repeater_scaling +
         (wr/len) * repeater_spacing * input_cap * repeater_scaling +
         0.5 * (wr/len) * (wc/len)* repeater_spacing * repeater_spacing;
#define Ishort_ckt 65e-6 /* across all tech Ref:Banerjee et al. {IEEE TED} */
    short_ckt = vdd_periph_global * minimum_width_nmos * Ishort_ckt * 1.0986 * 
            repeater_scaling * tc;

    PRINTDW(fprintf(stderr, "wire_energy: (energy of a repater + stretch"
         " of wire between repeaters) len %g, switching energy %g,"
            " time constant %g, short circuit energy %g\n", len,
            switching, tc, short_ckt));

    wire_stats->power.dynamic = ((len/repeater_spacing)*(switching + short_ckt));
    wire_stats->power.leakage = ((len/repeater_spacing)*
                    1.5*vdd_periph_global*
                    I_off_n[VTHMODEL][temper-300]*
//                    I_off_n[0][temper-300]*
                    minimum_width_nmos*repeater_scaling);
}

#ifndef __linux__
void bcopy (char *a, char *b, int c) {
	while (c>0) {
		*b++ = *a++;
		c--;
	}
}
#endif

void
calc_wire_stats2 (enum wire_type wire_model, wire_stats_t *wire_st) 
{
    double sf; /* scale factor */
    /* 
     * make sure that the length of the wire is already 
     * calculated and is in reasonable range
     */
    if (wire_model != Low_swing) {
        sf = ((wire_model == Semi_global) ? .5:1.0);
        /* delay increases linearly with length for wires with 
         * repeaters */
        wire_st->wire_pda.delay /* delay in s*/ = 
            2.13 * sqrt(wire_res(wire_st->wire_length , sf, sf) * 
            wire_cap(wire_st->wire_length, sf, sf) * FO4*1e-12/3); 
            /* FO1 = FO4/3 Ref: Future of wires FIXME*/
        wire_st->wire_spacing =
                sf * wire_pitch[ITRS_GLOBAL][VTHMODEL]*1e-6/2 /* (m) */;
        wire_st->wire_width = wire_st->wire_spacing; 
        wire_energy2(wire_st);
        if (wire_model == Global_10) {
            wire_st->wire_pda.delay += wire_st->wire_pda.delay * .1;
            wire_st->wire_pda.power.dynamic /= 2;
            wire_st->wire_pda.power.leakage /= 2;
        }
        else if (wire_model == Global_20) {
            wire_st->wire_pda.delay += wire_st->wire_pda.delay * .2;
            wire_st->wire_pda.power.dynamic *= .4;
            wire_st->wire_pda.power.leakage *= .4;
        }
        else if (wire_model == Global_30) {
            wire_st->wire_pda.delay += wire_st->wire_pda.delay * .3;
            wire_st->wire_pda.power.dynamic *= .35;
            wire_st->wire_pda.power.leakage *= .35;
        }
    }
    else if (wire_model == Low_swing) {
        /* 
         * make sure that the length of the wire is already 
         * calculated and is in reasonable range
         */
        low_swing_model (wire_st);
    }
    else {
        assert(0); 
    }
}


void 
delay_output_driver_at_subarray(int deg_out_driver_muxing, double inrisetime, double *outrisetime, powerDef *power,
        dataout_htree_node *ptr_htree_node, int number_mats,
        double *delay, 
        double *delay_final_stage_subarray_output_driver)
{
    int j, flag_final_stage_subarray_output_driver;
    double c_load, rd, tf, this_delay, c_intrinsic;

    power->readOp.dynamic = 0;
    power->readOp.leakage = 0;
    power->writeOp.dynamic = 0;
    power->writeOp.leakage = 0;
    *delay = 0;
    *delay_final_stage_subarray_output_driver = 0;
    flag_final_stage_subarray_output_driver = 0;

   
	//delay of signal through pass-transistor to input of subarray output driver.
	rd = transreson(width_nmos_sense_amp_mux, NCH, 1);
	c_load = deg_out_driver_muxing * draincap(width_nmos_sense_amp_mux, NCH, 1, 1, DEFAULTHEIGHTCELL) + gatecap(ptr_htree_node->width_n[0] + 
		ptr_htree_node->width_p[0] + ptr_htree_node->width_nor2_n + ptr_htree_node->width_nor2_p, 0.0);
	tf = rd * c_load;
	this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
	*delay += this_delay;
	inrisetime = this_delay/(1.0 - 0.5);
	power->readOp.dynamic += c_load * 0.5 * vdd_periph_global * vdd_periph_global;
	power->readOp.leakage += 0;//for now, let leakage of the pass transistor be 0

    for(j = 0; j < ptr_htree_node->number_gates; ++j){
        if(j == 0){//NAND2 gate
			rd = transreson(ptr_htree_node->width_n[j], NCH, 2);
			if(ptr_htree_node->number_gates ==2){
				c_load = gatecap(ptr_htree_node->width_p[j+1], 0.0);//NAND2 drives PMOS of output stage
			}
			else{
				c_load = gatecap(ptr_htree_node->width_n[j+1] + ptr_htree_node->width_p[j+1], 0.0);//NAND2 drives inverter
			}
            c_intrinsic = draincap(ptr_htree_node->width_n[j], NCH, 2, 1, DEFAULTHEIGHTCELL) + 
                2 * draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
            tf = rd * (c_intrinsic + c_load);
            power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
            power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_n[j], 
                    ptr_htree_node->width_p[j],  temper) * 0.5 * vdd_periph_global;
			power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_nor2_n, 
					ptr_htree_node->width_nor2_p, temper) * 0.5 * vdd_periph_global;
        }
        else if(j == ptr_htree_node->number_gates - 1){//PMOS
			flag_final_stage_subarray_output_driver = 1;
			rd = transreson(ptr_htree_node->width_p[j], PCH, 1);
			c_load = ptr_htree_node->c_wire_load + ptr_htree_node->c_gate_load +
				draincap(ptr_htree_node->width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL);
			c_intrinsic = draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load) + ptr_htree_node->r_wire_load * 
				(ptr_htree_node->c_wire_load / 2 + ptr_htree_node->c_gate_load +
				draincap(ptr_htree_node->width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL));
			power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_n[j], 
				ptr_htree_node->width_p[j], temper) * 0.5 * vdd_periph_global;
			flag_final_stage_subarray_output_driver = 1;
		}
		else{//inverter
			rd = transreson(ptr_htree_node->width_n[j], NCH, 1);
			if(j == ptr_htree_node->number_gates - 2){//inverter driving PMOS of output stage
				c_load = gatecap(ptr_htree_node->width_p[j+1], 0.0);
			}
			else{//inverter driving inverter
				c_load = gatecap(ptr_htree_node->width_n[j+1] + ptr_htree_node->width_p[j+1], 0.0);
			}
			c_intrinsic = draincap(ptr_htree_node->width_p[j], PCH, 1, 1, DEFAULTHEIGHTCELL) + 
				draincap(ptr_htree_node->width_n[j], NCH, 1, 1, DEFAULTHEIGHTCELL);
			tf = rd * (c_intrinsic + c_load);
			power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
			power->readOp.leakage += cmos_ileakage(ptr_htree_node->width_n[j] + ptr_htree_node->width_n[j] / 2, 
				ptr_htree_node->width_p[j] + ptr_htree_node->width_p[j] /2, temper) * 0.5 * vdd_periph_global;
		}
        this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
        *delay += this_delay;
        if(flag_final_stage_subarray_output_driver){
            *delay_final_stage_subarray_output_driver = this_delay;
        }
        inrisetime = this_delay/(1.0 - 0.5);
    }

    *outrisetime = inrisetime;

}

void
delay_routing_to_bank(int number_repeaters_htree_route_to_bank, double length_htree_route_to_bank,
        double sizing_repeater_htree_route_to_bank, double inrisetime,double *outrisetime,  
		double *delay, powerDef *power)
{
    double rd, c_intrinsic, c_load, r_wire, tf, this_delay;
    int i;
    *delay = 0;
    //Add delay of cascade of inverters that drives the first repeater. FIX this.
    if(length_htree_route_to_bank > 0){
        for(i = 0; i < number_repeaters_htree_route_to_bank; ++i){
            rd = transreson(sizing_repeater_htree_route_to_bank * minimum_width_nmos, NCH, 1);
            c_intrinsic = draincap(sizing_repeater_htree_route_to_bank * minimum_width_pmos, PCH, 1, 1, DEFAULTHEIGHTCELL) + 
                draincap(sizing_repeater_htree_route_to_bank * minimum_width_nmos, NCH, 1, 1, DEFAULTHEIGHTCELL);
            c_load = sizing_repeater_htree_route_to_bank * gatecap(minimum_width_nmos +
				minimum_width_pmos, 0);
            //The last repeater actually sees the load at the input of the bank as well. FIX this.
            //if(i == number_repeaters_htree_route_to_bank - 1){
            //c_load = ptr_intcnt_seg->c_gate_load;
            //}

            r_wire = (length_htree_route_to_bank / number_repeaters_htree_route_to_bank) *
                wire_outside_mat_r_per_micron;
            tf = rd * c_intrinsic + (rd + r_wire) * c_load;
            this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
            *delay += this_delay;
            inrisetime = this_delay / (1.0 - 0.5);	
            power->readOp.dynamic += (c_intrinsic + c_load) * 0.5 *	vdd_periph_global * vdd_periph_global;
            power->readOp.leakage += cmos_ileakage(
                    sizing_repeater_htree_route_to_bank * minimum_width_nmos,
                    sizing_repeater_htree_route_to_bank * minimum_width_pmos,  temper) * 0.5 * vdd_periph_global;
		}
		//From Ron Ho's thesis, energy-delay optimal repeaters improve dynamic energy approximately by 
		//30% at the expense of approx 10% increase in delay. We simply apply this correction at the end here. 
		//Assuming that leakage also improves by 30%. Doing it like this not right. FIX this - arrive at these numbers by deriving sizing, number of
		//repeaters....
		//*delay = 1.1 * (*delay);
		//power->readOp.dynamic = 0.7 * power->readOp.dynamic;
		//power->readOp.leakage = 0.7 * power->readOp.leakage;
    }
    *outrisetime = inrisetime;
}


void 
initialize_driver(driver *ptr_driver)
{
	int i;
	ptr_driver->delay = 0;
	ptr_driver->power.readOp.dynamic = 0;
	ptr_driver->power.readOp.leakage = 0;
	ptr_driver->power.writeOp.dynamic = 0;
	ptr_driver->power.writeOp.leakage = 0;
	ptr_driver->min_number_gates = 2;
	for(i = 0; i < MAX_NUMBER_GATES_STAGE; ++i){
		ptr_driver->width_n[i] = 0;
		ptr_driver->width_p[i] = 0;
	}
}

void 
compute_widths_driver(driver *ptr_driver)
{ 
	double c_load, f, F, c_input;
	int i;
	c_load = ptr_driver->c_gate_load + ptr_driver->c_wire_load;
	ptr_driver->width_n[0] = minimum_width_nmos;
	ptr_driver->width_p[0] = minimum_width_pmos;
	F = c_load / gatecap(ptr_driver->width_n[0] + ptr_driver->width_p[0], 0);
	ptr_driver->number_gates = (int) (log(F) / log(fopt) + 0.5);
	f = pow(F, 1.0 / ptr_driver->number_gates);
	if(ptr_driver->number_gates%2 != 0){
            --ptr_driver->number_gates;
	}
	if(ptr_driver->number_gates < ptr_driver->min_number_gates){
		ptr_driver->number_gates = ptr_driver->min_number_gates;
	}
	f = pow(F, 1.0 / ptr_driver->number_gates);

	i = ptr_driver->number_gates - 1;
	c_input = c_load /  f;
	ptr_driver->width_n[i] = (1.0 / 3.0) * c_input / gatecap(1, 0);
	ptr_driver->width_p[i] = 2 * ptr_driver->width_n[i];

	if(ptr_driver->width_n[i] > MAX_NMOS_WIDTH){
		c_load = gatecap(MAX_NMOS_WIDTH + 2 * MAX_NMOS_WIDTH, 0);
		F = c_load / gatecap(ptr_driver->width_n[0] + ptr_driver->width_p[0], 0);
		ptr_driver->number_gates = (int) (log(F) / log(fopt) + 0.5) + 1;
		if(ptr_driver->number_gates%2 != 0){
				--ptr_driver->number_gates;
		}
		if(ptr_driver->number_gates < ptr_driver->min_number_gates){
			ptr_driver->number_gates = ptr_driver->min_number_gates;
		}
		f = pow(F, 1.0 / ptr_driver->number_gates);
		i = ptr_driver->number_gates - 1;
		ptr_driver->width_n[i] = MAX_NMOS_WIDTH;
		ptr_driver->width_p[i] = 2 * ptr_driver->width_n[i];
	}

	for(i = ptr_driver->number_gates - 2; i >= 1; --i){
		ptr_driver->width_n[i] = ptr_driver->width_n[i+1] / f;
		ptr_driver->width_p[i] = 2 * ptr_driver->width_n[i];
	}
}

void 
delay_driver(driver *ptr_driver, double inrisetime, double *outrisetime)
{
	int i;
    double rd, c_load, c_intrinsic, tf, this_delay;
    this_delay = 0;
	for(i = 0; i < ptr_driver->number_gates - 1; ++i){
		rd = transreson(ptr_driver->width_n[i], NCH, 1);
		c_load = gatecap(ptr_driver->width_n[i+1] + ptr_driver->width_p[i+1], 0.0);
		c_intrinsic = draincap(ptr_driver->width_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) +  
			draincap(ptr_driver->width_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
		tf = rd * (c_intrinsic + c_load);
		this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
		ptr_driver->delay += this_delay;
		inrisetime = this_delay / (1.0 - 0.5);
		ptr_driver->power.readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
		ptr_driver->power.readOp.leakage += cmos_ileakage(ptr_driver->width_n[i], 
			ptr_driver->width_p[i], temper) * 0.5 * vdd_periph_global;
    }
	i = ptr_driver->number_gates - 1;
	c_load = ptr_driver->c_gate_load + ptr_driver->c_wire_load;
	rd = transreson(ptr_driver->width_n[i], NCH, 1);
	c_intrinsic = draincap(ptr_driver->width_p[i], PCH, 1, 1, DEFAULTHEIGHTCELL) +  
		draincap(ptr_driver->width_n[i], NCH, 1, 1, DEFAULTHEIGHTCELL);
	tf = rd * (c_intrinsic + c_load) + ptr_driver->r_wire_load * (ptr_driver->c_wire_load  / 2 +
		ptr_driver->c_gate_load);
	this_delay = horowitz(inrisetime, tf, 0.5, 0.5, RISE);
	ptr_driver->delay += this_delay;
	*outrisetime = this_delay / (1.0 - 0.5);
	ptr_driver->power.readOp.dynamic += (c_intrinsic + c_load) * 0.5 * vdd_periph_global * vdd_periph_global;
	ptr_driver->power.readOp.leakage += cmos_ileakage(ptr_driver->width_n[i], ptr_driver->width_p[i],
		temper) * 0.5 * vdd_periph_global;
}


/* check whether the input ndwl, ndbl, nspd, ndcm, and ndsam values make sense
 * for the given cache parameters
 */
int
is_partition_valid (input_params_t *parameters, int tag_array, int Ndwl, int Ndbl, double Nspd, int Ndcm, int Ndsam)
{
    int dram; 

    int sub_hor_bits /* no. of rowsin a sub-array */, 
        sub_ver_bits /* no. of columns in a sub-array */;

    int num_sets, num_sets_pbank, set_bits, tag_bits, rem;

    int bits_read_per_mat;
    
    dram = parameters->dram;

    /* DRAM's degree of bitline muxing can never be greater than 1 */
    if (dram && !tag_array && Ndcm > 1) {
        return 0;
    }
    /* Both tag and data array should have atleast 4 sub-arrays */
    if (Ndwl == 1 || Ndbl == 1) {
        return 0;
    }

    if (Nspd == 0) return 0;

    if (tag_array) {
        num_sets = (int)(parameters->cache_size/
                (parameters->block_size * 
                parameters->tag_associativity));
    }
    else {
        num_sets = (int)(parameters->cache_size/
                (parameters->block_size * 
                parameters->data_associativity));
    }

    num_sets_pbank = num_sets/parameters->uca_banks;

    if (tag_array) {
        if(parameters->force_tag == 1){
            tag_bits = parameters->tag_size;
        }
        else{
            tag_bits = ADDRESS_BITS + EXTRA_TAG_BITS - 
                (int)logtwo((double)(parameters->cache_size)) + 
                (int)(ceil(logtwo((double)(parameters->tag_associativity)))) - 
                (int)(logtwo(parameters->uca_banks));
        }
        rem = tag_bits%4;
        if(rem != 0){
            tag_bits = tag_bits + 4 - rem;
        }
        set_bits = tag_bits * parameters->tag_associativity;
    }
    else {
        set_bits = parameters->block_size * 
                parameters->data_associativity * 8;
    }

    sub_ver_bits = (int)(num_sets_pbank/(Ndbl*Nspd) + EPSILON);

    sub_hor_bits = (int)((set_bits * Nspd / Ndwl) + EPSILON);

    /* sub-arrays cant be arbitrarily large */
    if (sub_ver_bits > MAXSUBARRAYROWS || 
        sub_ver_bits < MINSUBARRAYROWS) {
        return 0;
    }
    if (sub_hor_bits > MAXSUBARRAYCOLS || 
        sub_hor_bits < MINSUBARRAYCOLS) {
        return 0;
    }

    if ((Ndwl*Ndbl) > MAXSUBARRAYS) { /* there is an upper limit 
                                      on tot. no. sub-arrays */
        return 0;
    }
    if (Ndwl*Ndbl < 4) return 0;

    bits_read_per_mat = 4/* there are 4 sub-arrays/mat*/ * 
                        sub_hor_bits/(Ndcm * Ndsam);

    if(bits_read_per_mat < 4) return 0;


    if(tag_array) {
        if (bits_read_per_mat * (Ndwl/2) !=
            (tag_bits * parameters->tag_associativity)) {
            return 0;
        }
    }
    else {
        if(parameters->fast_access == 1) {
            if(bits_read_per_mat * (Ndwl/2) != 
                    (parameters->output_width * parameters->data_associativity)) {
                return 0;
            }
        }
        else {
            if(bits_read_per_mat * (Ndwl/2) != 
                    (parameters->output_width)) {
                return 0;
            }
        }
    }
    return 1;
}


/*======================================================================*/

int calculate_time(input_params_t *parameters, int NSubbanks, int is_tag, int pure_ram, double Nspd, int Ndwl, 
        int Ndbl, int Ndcm, int Ndsam, mem_array *ptr_array, int flag_results_populate,
        results_mem_array *ptr_results, final_results *ptr_fin_res)
{
    wire_stats_t wst;
    arearesult_type arearesult_temp;

    double length_wire_htree_node[MAX_NUMBER_GATES_STAGE];

    predecoder_block row_predec_blk_1, row_predec_blk_2,  bit_mux_predec_blk_1, bit_mux_predec_blk_2,
        senseamp_mux_predec_blk_1, senseamp_mux_predec_blk_2, dummy_way_select_predec_blk_1, 
		dummy_way_select_predec_blk_2;

    decoder row_dec, bit_mux_dec, senseamp_mux_dec;

	driver bitline_precharge_eq_driver;

    predecoder_block_driver row_predec_blk_driver_1, row_predec_blk_driver_2, bit_mux_predec_blk_driver_1, 
        bit_mux_predec_blk_driver_2, senseamp_mux_predec_blk_driver_1, senseamp_mux_predec_blk_driver_2,
		way_select_driver_1, dummy_way_select_driver_2;

    addr_datain_htree_node vertical_addr_din_htree_node[MAX_NUMBER_HTREE_NODES], 
                   horizontal_addr_din_htree_node[MAX_NUMBER_HTREE_NODES];
    dataout_htree_node dout_htree_node[MAX_NUMBER_HTREE_NODES], 
                   subarray_output_htree_node;
    powerDef tot_power_addr_vertical_htree, tot_power_datain_vertical_htree, tot_power_bitlines, 
        tot_power_sense_amps, tot_power_subarray_output_drivers, 
        tot_power_dataout_vertical_htree, tot_power_comparators;

    powerDef tot_power, tot_power_row_predecode_block_drivers,
             tot_power_bit_mux_predecode_block_drivers, 
             tot_power_senseamp_mux_predecode_block_drivers,
             tot_power_row_predecode_blocks, tot_power_bit_mux_predecode_blocks,
             tot_power_senseamp_mux_predecode_blocks, tot_power_row_decoders, 
             tot_power_bit_mux_decoders, tot_power_senseamp_mux_decoders, 
			 tot_power_bitlines_precharge_eq_driver;

	
    double access_time = 0;
    double bitline_data = 0.0;
    powerDef sense_amp_data_power,  bitline_data_power;
    double cycle_time = 0.0;
    double outrisetime = 0.0, inrisetime = 0.0;
    double Tpre;
    double writeback_delay;
    int number_rows_subarray, number_cols_subarray, deg_bitline_muxing, deg_senseamp_muxing_non_associativity;


    double outrisetime_nand2_row_decode_path_1, outrisetime_nand3_row_decode_path_1,
           outrisetime_nand2_bit_mux_decode_path_1, outrisetime_nand3_bit_mux_decode_path_1,
           outrisetime_nand2_senseamp_mux_decode_path_1, outrisetime_nand3_senseamp_mux_decode_path_1,
           outrisetime_nand2_row_decode_path_2, outrisetime_nand3_row_decode_path_2,
           outrisetime_nand2_bit_mux_decode_path_2, outrisetime_nand3_bit_mux_decode_path_2,
           outrisetime_nand2_senseamp_mux_decode_path_2, outrisetime_nand3_senseamp_mux_decode_path_2;

    double delay_before_row_decoder_nand2_path_1, delay_before_row_decoder_nand3_path_1,
           delay_before_row_decoder_nand2_path_2, delay_before_row_decoder_nand3_path_2,
           delay_before_bit_mux_decoder_nand2_path_1, delay_before_bit_mux_decoder_nand3_path_1,
           delay_before_bit_mux_decoder_nand2_path_2, delay_before_bit_mux_decoder_nand3_path_2,
           delay_before_senseamp_mux_decoder_nand2_path_1, delay_before_senseamp_mux_decoder_nand3_path_1,
           delay_before_senseamp_mux_decoder_nand2_path_2, delay_before_senseamp_mux_decoder_nand3_path_2;

    double c_wordline_driver_load, c_gate_load_bit_mux_decoder_output, c_wire_load_bit_mux_decoder_output, c_load_bit_mux_decoder_output,
           c_gate_load_senseamp_mux_decoder_output, c_wire_load_senseamp_mux_decoder_output, c_load_senseamp_mux_decoder_output,
           r_wire_wordline_driver_output, r_wire_bit_mux_decoder_output, r_wire_senseamp_mux_decoder_output,
           c_wire_row_predecode_block_output, r_wire_row_predecode_block_output,
           c_wire_bit_mux_predecode_block_output, r_wire_bit_mux_predecode_block_output, 
           c_wire_senseamp_mux_predecode_block_output, r_wire_senseamp_mux_predecode_block_output;

    double row_dec_outrisetime, bit_mux_dec_outrisetime, senseamp_mux_dec_outrisetime;

    double delay_addr_din_vertical_htree, delay_dout_vertical_htree, delay_before_subarray_output_driver,
           delay_from_subarray_output_driver_to_output;


    double max_delay_before_row_decoder, max_delay_before_bit_mux_decoder, max_delay_before_senseamp_mux_decoder;

    double row_dec_inrisetime, bit_mux_dec_inrisetime, senseamp_mux_dec_inrisetime;

    double delay_data_access_row_path, delay_data_access_col_path, delay_data_access_out_drv_path,
           temp_delay_data_access_path, delay_data_access_path; 

    int number_addr_bits_nand2_row_decode_path_1, number_addr_bits_nand3_row_decode_path_1, 
        number_addr_bits_nand2_row_decode_path_2, number_addr_bits_nand3_row_decode_path_2,
        number_addr_bits_nand2_bit_mux_decode_path_1, number_addr_bits_nand3_bit_mux_decode_path_1,
        number_addr_bits_nand2_bit_mux_decode_path_2, number_addr_bits_nand3_bit_mux_decode_path_2,
        number_addr_bits_nand2_senseamp_mux_decode_path_1, number_addr_bits_nand3_senseamp_mux_decode_path_1,
        number_addr_bits_nand2_senseamp_mux_decode_path_2, number_addr_bits_nand3_senseamp_mux_decode_path_2;

    powerDef power_dout_htree, power_output_drivers_at_subarray, power_addr_datain_htree;

    int number_addr_bits_mat, number_datain_bits_mat, number_dataout_bits_mat, 
        number_subarrays, number_mats, number_sense_amps_subarray,
	    number_output_drivers_subarray,
        number_way_select_signals_mat = 0;

    area_type subarray, mat, bank, all_banks;

    point_to_point_interconnect_segment horizontal_addr_intcnt_segment_within_bank,
                                        horizontal_datain_intcnt_segment_within_bank;
    double  delay_dout_horizontal_htree;

    int number_subarrays_horizontal_direction, number_subarrays_vertical_direction, 
		number_vertical_htree_nodes, number_horizontal_htree_nodes, 
        number_mats_horizontal_direction, number_mats_vertical_direction, 
        number_activated_mats_horizontal_direction, way_select = 0;
    int number_decode_gates_driven_per_predecode_output, 
        number_addr_bits_row_decode, number_addr_bits_bit_mux_decode, 
        number_addr_bits_senseamp_mux_decode_non_associativity,  number_subbanks, 
        number_subbanks_decode, number_addr_bits_routed_to_bank,
        number_comparator_bits_routed_to_bank, number_data_bits_routed_to_bank, 
        number_bits_routed_to_bank, number_repeaters_htree_route_to_bank,
		number_dataout_bits_subbank;

    double sizing_repeater_htree_route_to_bank, length_htree_route_to_bank;

    double vertical_htree_seg_length, horizontal_htree_seg_length, delay_addr_din_horizontal_htree;
    powerDef power_addr_datain_horizontal_htree;
	powerDef power_addr_datain_horizontal_htree_node, power_addr_datain_vertical_htree_node;


    int rem, tagbits,  k, htree_seg_multiplier, number_mats_to_cover_in_this_htree_segment;
    double Cbitrow_drain_cap, Cbitline;
    double delay_route_to_bank, wordline_data, delay_subarray_output_driver, delay_final_stage_subarray_output_driver,
           delay_sense_amp,  multisubbank_interleave_cycle_time;
    powerDef tot_power_addr_horizontal_htree, tot_power_datain_horizontal_htree,
             tot_power_dataout_horizontal_htree, //tot_power_comparators,
             tot_power_routing_to_bank;
    double comparator_delay;
    powerDef comparator_power, power_routing_to_bank;

    double area_all_dataramcells;

    double c_load, rd, c_intrinsic, tf, wordline_reset_delay;
	double dummy_precharge_outrisetime;

	double r_bitline_precharge, r_bitline, bitline_restore_delay;

	double temp, v_storage_worst, dram_array_availability;

	double leakage_power_cc_inverters_sram_cell, leakage_power_access_transistors_read_write_port_sram_cell,
			leakage_power_read_only_port_sram_cell;

    if((is_dram)&&(!is_tag)&&(Ndcm > 1)){goto partition_not_valid;}
    if(Ndwl == 1) {goto partition_not_valid;}
    if(Ndbl == 1) {goto partition_not_valid;}
    tagbits = 0;//if data array, let tagbits = 0
    if(is_tag){
        if(parameters->force_tag == 1){
			tagbits = parameters->tag_size;
        }
        else{
            tagbits = ADDRESS_BITS + EXTRA_TAG_BITS - (int)logtwo((double)(parameters->cache_size)) + 
                (int)(ceil(logtwo((double)(parameters->tag_associativity)))) - (int)(logtwo(NSubbanks));
        }
        rem = tagbits%4;
        if(rem != 0){
            tagbits = tagbits + 4 - rem;
        }
        number_rows_subarray = (int)(parameters->cache_size / (parameters->uca_banks *
                    parameters->block_size * parameters->tag_associativity * Ndbl * Nspd) + EPSILON);
        number_cols_subarray = (int)((tagbits * parameters->tag_associativity * Nspd / Ndwl) + EPSILON);
    }
    else{
        number_rows_subarray = (int)(parameters->cache_size / (parameters->uca_banks * 
                    parameters->block_size * parameters->data_associativity * Ndbl * Nspd) + EPSILON);
        number_cols_subarray = (int)((8 * parameters->block_size * parameters->data_associativity * Nspd / Ndwl) + EPSILON);
    }		
    if(number_rows_subarray < MINSUBARRAYROWS) {goto partition_not_valid;}
    if(number_rows_subarray > MAXSUBARRAYROWS) {goto partition_not_valid;}
    if(number_cols_subarray < MINSUBARRAYCOLS) {goto partition_not_valid;}
    if(number_cols_subarray > MAXSUBARRAYCOLS) {goto partition_not_valid;}
    number_subarrays = Ndwl * Ndbl;	
    if(number_subarrays < 4) {goto partition_not_valid;}//number of subarrays is at least 4 i.e. at least one mat.
    if(number_subarrays > MAXSUBARRAYS) {goto partition_not_valid;}

	//Calculate wire parameters
	if(is_tag){
		height_cell = BitHeight_sram + 2 * wire_local_pitch * (parameters->rw_ports - 1 + 
			parameters->excl_read_ports);
		width_cell = BitWidth_sram + 2 * wire_local_pitch * (parameters->rw_ports - 1 + 
			(parameters->excl_read_ports - parameters->single_ended_read_ports)) + 
			wire_local_pitch * parameters->single_ended_read_ports;
	}
	else{
		if(is_dram){
			height_cell = BitHeight_dram;
			width_cell = BitWidth_dram;
		}
		else{
			height_cell = BitHeight_sram + 2 * wire_local_pitch * (parameters->excl_write_ports + 
				parameters->rw_ports - 1 + parameters->excl_read_ports);
			width_cell = BitWidth_sram + 2 * wire_local_pitch * (parameters->rw_ports - 1 + 
				(parameters->excl_read_ports - parameters->single_ended_read_ports) + 
				parameters->excl_write_ports) + wire_local_pitch * parameters->single_ended_read_ports;
		}
	}

	Cbitmetal = height_cell * wire_local_c_per_micron;
	Cwordmetal = width_cell * wire_local_c_per_micron;
	Rbitmetal = height_cell * wire_local_r_per_micron;
	Rwordmetal = width_cell * wire_local_r_per_micron;

    if(is_dram){
        deg_bitline_muxing = 1;
		is_access_transistor = 1;
        Cbitrow_drain_cap = draincap(Wmemcella_dram, NCH, 1, 0, width_cell) / 2.0;
		is_access_transistor = 0;
		Cbitline = number_rows_subarray / 2 * Cbitrow_drain_cap + number_rows_subarray * 
			Cbitmetal;
        Vbitsense = (vdd_dram_cell/2) * c_dram_cell /(c_dram_cell + Cbitline);
        if(Vbitsense < VBITSENSEMIN) goto partition_not_valid;
		Vbitsense = VBITSENSEMIN;
		v_storage_worst = vdd_dram_cell / 2 - VBITSENSEMIN * (c_dram_cell + Cbitline) / c_dram_cell;
		dram_refresh_period = 1.1 * c_dram_cell * v_storage_worst / I_off_dram_cell_worst_case_length_temp;
    }
    else{//SRAM
        deg_bitline_muxing = Ndcm;
		is_sram_cell = 1; 
		Cbitrow_drain_cap = draincap(Wmemcella_sram, NCH, 1, 0, width_cell) / 2.0;	/* due to shared contact */
		is_sram_cell = 0; 
		Cbitline = number_rows_subarray * (Cbitrow_drain_cap + Cbitmetal);
		dram_refresh_period = 0;
    }
    number_subarrays_horizontal_direction = Ndwl;
    number_subarrays_vertical_direction = Ndbl;
    number_mats_horizontal_direction = number_subarrays_horizontal_direction / 2;
    number_mats_vertical_direction = number_subarrays_vertical_direction / 2;
    number_mats = number_mats_horizontal_direction * number_mats_vertical_direction;
    number_dataout_bits_mat = MAX(4 * number_cols_subarray / (deg_bitline_muxing * Ndsam), 1);
    if(number_dataout_bits_mat < 4) {goto partition_not_valid;}
   

    if(!is_tag){
		if(parameters->fast_access == 1){
			number_dataout_bits_subbank = parameters->output_width * parameters->data_associativity;
		}
		else{
			number_dataout_bits_subbank = parameters->output_width;
		}
//        if(Ndsam < parameters->data_associativity) {goto partition_not_valid;} //FIXME
        if(number_dataout_bits_mat * number_mats_horizontal_direction !=
                number_dataout_bits_subbank) {goto partition_not_valid;}
        deg_senseamp_muxing_non_associativity = MAX(Ndsam / parameters->data_associativity, 1);
    }
    else{
        if(number_dataout_bits_mat * number_mats_horizontal_direction !=
                tagbits * parameters->tag_associativity) {goto partition_not_valid;}
        deg_senseamp_muxing_non_associativity = Ndsam;
    }

    number_activated_mats_horizontal_direction = number_mats_horizontal_direction;
    number_datain_bits_mat = number_dataout_bits_mat;
    number_addr_bits_row_decode = (int)(logbasetwo((double)(number_rows_subarray)));
    number_addr_bits_bit_mux_decode = (int)(logbasetwo((double)(deg_bitline_muxing)));
    number_addr_bits_senseamp_mux_decode_non_associativity = 
        (int)(logbasetwo((double)(deg_senseamp_muxing_non_associativity)));
    number_subbanks = number_mats / number_activated_mats_horizontal_direction;
    number_subbanks_decode = (int)(logbasetwo((double)(number_subbanks)));
    number_addr_bits_mat = number_addr_bits_row_decode + number_addr_bits_bit_mux_decode +
        number_addr_bits_senseamp_mux_decode_non_associativity;
    number_addr_bits_routed_to_bank = number_addr_bits_mat + number_subbanks_decode;
    number_comparator_bits_routed_to_bank = 0;
    if((!is_tag)&&(!pure_ram)){//data array of cache
        number_comparator_bits_routed_to_bank = parameters->data_associativity;
    }
    if(is_tag){
        number_data_bits_routed_to_bank = tagbits + parameters->data_associativity;//input to tag
        //array= tagbits, output from tag array = parameters->data_associativity number of comparator
        //bits for set-associative cache. parameters->data_associativity = 1 valid bit for direct-mapped
        //cache or sequential-access cache.
    }
    else{
        number_data_bits_routed_to_bank = parameters->output_width * 2;//datain and dataout
    }
    number_bits_routed_to_bank = number_addr_bits_routed_to_bank + number_data_bits_routed_to_bank +
        number_comparator_bits_routed_to_bank;

    number_sense_amps_subarray = number_cols_subarray / deg_bitline_muxing;
	number_output_drivers_subarray = number_sense_amps_subarray / Ndsam;
    if((!is_tag)&&(parameters->data_associativity > 1)){
        way_select = parameters->data_associativity;
        number_way_select_signals_mat = parameters->data_associativity;
    }		

    if(is_dram){
        is_access_transistor = 1;
        c_wordline_driver_load = (gatecappass(Wmemcella_dram, BitWidth_dram) + Cwordmetal) * 	number_cols_subarray;
        is_access_transistor = 0;
    }
    else{
        is_sram_cell = 1;
        c_wordline_driver_load = (gatecappass(Wmemcella_sram,(BitWidth_sram-2*Wmemcella_sram)/2.0) + gatecappass(Wmemcella_sram,
                    (BitWidth_sram-2*Wmemcella_sram)/2.0) + Cwordmetal) * number_cols_subarray;
        is_sram_cell = 0;
    }

    reset_powerDef(&tot_power);
    reset_powerDef(&tot_power_routing_to_bank);
    reset_powerDef(&tot_power_addr_horizontal_htree);
    reset_powerDef(&tot_power_datain_horizontal_htree);
    reset_powerDef(&tot_power_dataout_horizontal_htree);
    reset_powerDef(&tot_power_addr_vertical_htree);
	reset_powerDef(&tot_power_datain_vertical_htree);
    reset_powerDef(&tot_power_row_predecode_block_drivers);
    reset_powerDef(&tot_power_bit_mux_predecode_block_drivers);
    reset_powerDef(&tot_power_senseamp_mux_predecode_block_drivers);
    reset_powerDef(&tot_power_row_predecode_blocks);
    reset_powerDef(&tot_power_bit_mux_predecode_blocks);
    reset_powerDef(&tot_power_senseamp_mux_predecode_blocks);
    reset_powerDef(&tot_power_row_decoders);
    reset_powerDef(&tot_power_bit_mux_decoders);
    reset_powerDef(&tot_power_senseamp_mux_decoders);
    reset_powerDef(&tot_power_bitlines);
	reset_powerDef(&tot_power_sense_amps);
    reset_powerDef(&tot_power_subarray_output_drivers);
    reset_powerDef(&tot_power_dataout_vertical_htree);
    reset_powerDef(&tot_power_comparators);
    refresh_power = 0;

	inrisetime = 0;
	rd = transreson(minimum_width_nmos, NCH, 1);
    c_load = gatecap(minimum_width_nmos + minimum_width_pmos, 0.0);
    tf = rd * c_load;
	kinv = horowitz(inrisetime, tf, 0.5, 0.5, RISE);

    subarray = subarraymem_area(number_rows_subarray, number_cols_subarray, number_subarrays,
            parameters->rw_ports,parameters->excl_read_ports, parameters->excl_write_ports,
            parameters->single_ended_read_ports);	
    height_subarray = subarray.height;
    width_subarray = subarray.width;
    area_all_dataramcells = subarray.height * subarray.width * number_subarrays *
        parameters->uca_banks;

	//number_vertical_htree_nodes = 1 when number_mats_vertical_direction = 1 or 2
	if(number_mats_vertical_direction > 1){
		number_vertical_htree_nodes = (int) (logbasetwo((double)(number_mats_vertical_direction / 2))) + 1;
	}
	else{//number_mats_vertical_direction = 1
		number_vertical_htree_nodes = 1;
	}
    number_horizontal_htree_nodes = (int) (logbasetwo((double)(number_mats_horizontal_direction))) + 1;
    compute_widths_subarray_output_driver(&subarray_output_htree_node, number_vertical_htree_nodes);

	if(deg_bitline_muxing > 1){
		c_gate_load_bit_mux_decoder_output = (4 * number_cols_subarray / deg_bitline_muxing) * gatecap(width_nmos_bit_mux, 0);
		c_wire_load_bit_mux_decoder_output = 2 * number_cols_subarray * wire_inside_mat_c_per_micron * width_cell;
		c_load_bit_mux_decoder_output = c_gate_load_bit_mux_decoder_output + c_wire_load_bit_mux_decoder_output;
	}
	else{
		c_gate_load_bit_mux_decoder_output = 0;
		c_wire_load_bit_mux_decoder_output = 0;
		c_load_bit_mux_decoder_output = 0;
	}

	if(Ndsam > 1){
		c_gate_load_senseamp_mux_decoder_output = (4 * number_sense_amps_subarray / Ndsam) * gatecap(width_nmos_sense_amp_mux, 0);
		c_wire_load_senseamp_mux_decoder_output = 2 * number_cols_subarray * wire_inside_mat_c_per_micron * width_cell;
		c_load_senseamp_mux_decoder_output = c_gate_load_senseamp_mux_decoder_output + c_wire_load_senseamp_mux_decoder_output;
	}
	else{
		c_gate_load_senseamp_mux_decoder_output = 0;
		c_wire_load_senseamp_mux_decoder_output = 0;
		c_load_senseamp_mux_decoder_output = 0;
	}

    r_wire_wordline_driver_output = number_cols_subarray * Rbitmetal;
    r_wire_bit_mux_decoder_output = number_cols_subarray * wire_inside_mat_r_per_micron * width_cell / 2;
    r_wire_senseamp_mux_decoder_output = number_cols_subarray * wire_inside_mat_r_per_micron * width_cell / 2;

    initialize_decoder(number_rows_subarray, 0, c_wordline_driver_load, r_wire_wordline_driver_output, &row_dec);
    initialize_decoder(deg_bitline_muxing, 0, c_load_bit_mux_decoder_output, r_wire_bit_mux_decoder_output, &bit_mux_dec);
    initialize_decoder(deg_senseamp_muxing_non_associativity, way_select, c_load_senseamp_mux_decoder_output, r_wire_senseamp_mux_decoder_output, &senseamp_mux_dec);

	is_wordline_transistor = 1;
    compute_widths_decoder(&row_dec);
	is_wordline_transistor = 0;
    compute_widths_decoder(&bit_mux_dec);
    compute_widths_decoder(&senseamp_mux_dec);

    c_wire_row_predecode_block_output = 2 * number_rows_subarray * wire_inside_mat_c_per_micron * height_cell;
    r_wire_row_predecode_block_output = number_rows_subarray * wire_inside_mat_r_per_micron * height_cell / 2;
    c_wire_bit_mux_predecode_block_output = 0;
    r_wire_bit_mux_predecode_block_output = 0;
    c_wire_senseamp_mux_predecode_block_output = 0;
    r_wire_senseamp_mux_predecode_block_output = 0;

    number_decode_gates_driven_per_predecode_output = 4;//Because of wordline decoders in 4 subarrays
    initialize_predecoder_blocks(number_rows_subarray, &row_predec_blk_1, &row_predec_blk_2, &row_dec,
            c_wire_row_predecode_block_output, r_wire_row_predecode_block_output, 
            number_decode_gates_driven_per_predecode_output);
    number_decode_gates_driven_per_predecode_output = 1;
    initialize_predecoder_blocks(deg_bitline_muxing, &bit_mux_predec_blk_1, &bit_mux_predec_blk_2, &bit_mux_dec, 
            c_wire_bit_mux_predecode_block_output, r_wire_bit_mux_predecode_block_output, 
            number_decode_gates_driven_per_predecode_output);
    number_decode_gates_driven_per_predecode_output = 1;
    initialize_predecoder_blocks(deg_senseamp_muxing_non_associativity, &senseamp_mux_predec_blk_1, &senseamp_mux_predec_blk_2,
            &senseamp_mux_dec, c_wire_senseamp_mux_predecode_block_output, r_wire_senseamp_mux_predecode_block_output, 
            number_decode_gates_driven_per_predecode_output);
	initialize_predecoder_blocks(1, &dummy_way_select_predec_blk_1, &dummy_way_select_predec_blk_2,
            &senseamp_mux_dec, 0, 0, 0);

    compute_widths_predecoder_block(&row_predec_blk_1);
    compute_widths_predecoder_block(&row_predec_blk_2);
    compute_widths_predecoder_block(&bit_mux_predec_blk_1);
    compute_widths_predecoder_block(&bit_mux_predec_blk_2);
    compute_widths_predecoder_block(&senseamp_mux_predec_blk_1);
    compute_widths_predecoder_block(&senseamp_mux_predec_blk_2);

    initialize_predecoder_block_drivers(number_rows_subarray, 0, 0, &row_predec_blk_driver_1,
            &row_predec_blk_driver_2, &row_predec_blk_1, &row_predec_blk_2,	&row_dec);
    initialize_predecoder_block_drivers(deg_bitline_muxing, 0, 0, &bit_mux_predec_blk_driver_1,
            &bit_mux_predec_blk_driver_2, &bit_mux_predec_blk_1, &bit_mux_predec_blk_2,	&bit_mux_dec);
    initialize_predecoder_block_drivers(deg_senseamp_muxing_non_associativity, 0, 0, &senseamp_mux_predec_blk_driver_1,
            &senseamp_mux_predec_blk_driver_2, &senseamp_mux_predec_blk_1, &senseamp_mux_predec_blk_2,	&senseamp_mux_dec);
	initialize_predecoder_block_drivers(1, 1, way_select, &way_select_driver_1,
            &dummy_way_select_driver_2, &dummy_way_select_predec_blk_1, &dummy_way_select_predec_blk_2, 
			&senseamp_mux_dec);

	compute_widths_predecoder_block_driver(&row_predec_blk_driver_1, &row_predec_blk_1, &row_dec, 0);
    compute_widths_predecoder_block_driver(&row_predec_blk_driver_2, &row_predec_blk_2, &row_dec, 0);
    compute_widths_predecoder_block_driver(&bit_mux_predec_blk_driver_1, &bit_mux_predec_blk_1, &bit_mux_dec, 0);
    compute_widths_predecoder_block_driver(&bit_mux_predec_blk_driver_2, &bit_mux_predec_blk_2, &bit_mux_dec, 0);
    compute_widths_predecoder_block_driver(&senseamp_mux_predec_blk_driver_1, &senseamp_mux_predec_blk_1, &senseamp_mux_dec, 0);
    compute_widths_predecoder_block_driver(&senseamp_mux_predec_blk_driver_2, &senseamp_mux_predec_blk_2, &senseamp_mux_dec, 0);
	compute_widths_predecoder_block_driver(&way_select_driver_1, &dummy_way_select_predec_blk_1, &senseamp_mux_dec, way_select);


	bitline_precharge_eq_driver.c_gate_load = gatecap(2 * width_pmos_bitline_precharge +
		width_pmos_bitline_equalization, 0);
	bitline_precharge_eq_driver.c_wire_load = number_cols_subarray * width_cell * 
		wire_outside_mat_c_per_micron;
	bitline_precharge_eq_driver.r_wire_load = number_cols_subarray * width_cell * 
		wire_outside_mat_r_per_micron;
	initialize_driver(&bitline_precharge_eq_driver);
	compute_widths_driver(&bitline_precharge_eq_driver);


    mat = area_mat(is_tag, number_rows_subarray, number_cols_subarray, 
            number_subarrays, deg_bitline_muxing, deg_senseamp_muxing_non_associativity, Ndsam,
			number_addr_bits_mat, number_datain_bits_mat, number_dataout_bits_mat,
            number_way_select_signals_mat, parameters,
            &row_predec_blk_1, &row_predec_blk_2, &bit_mux_predec_blk_1, &bit_mux_predec_blk_2, 
            &senseamp_mux_predec_blk_1, &senseamp_mux_predec_blk_2, &dummy_way_select_predec_blk_1, &row_dec, &bit_mux_dec, &senseamp_mux_dec,
            &row_predec_blk_driver_1, &row_predec_blk_driver_2, &bit_mux_predec_blk_driver_1,
            &bit_mux_predec_blk_driver_2, &senseamp_mux_predec_blk_driver_1, &senseamp_mux_predec_blk_driver_2,
            &way_select_driver_1, &subarray_output_htree_node);
    height_mat = mat.height;
    width_mat = mat.width;


    horizontal_htree_seg_length = width_mat / 2;

    initialize_addr_datain_htree(horizontal_addr_din_htree_node, 
                number_horizontal_htree_nodes, horizontal_htree_seg_length, 
                &delay_addr_din_horizontal_htree, 
                &power_addr_datain_horizontal_htree, 
                length_wire_htree_node);

    compute_widths_addr_datain_htree(horizontal_addr_din_htree_node, 
                number_horizontal_htree_nodes);

	initialize_addr_datain_htree_with_nodes_at_mat_interval(
        &horizontal_addr_datain_htree_at_mat_interval, &delay_addr_din_horizontal_htree, 
		&power_addr_datain_horizontal_htree, horizontal_htree_seg_length);

	compute_widths_addr_datain_htree_with_nodes_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval);

    vertical_htree_seg_length = height_mat / 2;
    initialize_addr_datain_htree(vertical_addr_din_htree_node, number_vertical_htree_nodes, 
            vertical_htree_seg_length, &delay_addr_din_vertical_htree, 
            &power_addr_datain_htree, length_wire_htree_node);	

    compute_widths_addr_datain_htree(vertical_addr_din_htree_node, number_vertical_htree_nodes);

	initialize_addr_datain_htree_with_nodes_at_mat_interval(
        &vertical_addr_datain_htree_at_mat_interval, &delay_addr_din_vertical_htree, 
		&power_addr_datain_htree, vertical_htree_seg_length);

	compute_widths_addr_datain_htree_with_nodes_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval);

    initialize_dataout_htree(dout_htree_node, number_vertical_htree_nodes, &delay_dout_vertical_htree, 
                             &power_dout_htree, length_wire_htree_node);

    compute_widths_dataout_htree(dout_htree_node, number_vertical_htree_nodes);


    bank = area_single_bank(number_rows_subarray, is_tag, parameters, &arearesult_temp, 
            number_horizontal_htree_nodes, number_vertical_htree_nodes, 
            number_mats_horizontal_direction, 
            number_mats_vertical_direction, 
            number_addr_bits_mat, number_way_select_signals_mat, tagbits, 
            number_datain_bits_mat, number_dataout_bits_mat, 
            horizontal_addr_din_htree_node, vertical_addr_din_htree_node, 
            dout_htree_node, &horizontal_addr_intcnt_segment_within_bank, 
            &horizontal_datain_intcnt_segment_within_bank, 
            &row_predec_blk_1, &row_predec_blk_2, &bit_mux_predec_blk_1, 
            &bit_mux_predec_blk_2, &senseamp_mux_predec_blk_1, &senseamp_mux_predec_blk_2, 
            &row_dec, &bit_mux_dec, &senseamp_mux_dec, &row_predec_blk_driver_1, 
            &row_predec_blk_driver_2, &bit_mux_predec_blk_driver_1,
            &bit_mux_predec_blk_driver_2, &senseamp_mux_predec_blk_driver_1, 
            &senseamp_mux_predec_blk_driver_2, &tot_power, 
            &tot_power_row_predecode_block_drivers, 
            &tot_power_bit_mux_predecode_block_drivers, 
            &tot_power_senseamp_mux_predecode_block_drivers,
            &tot_power_row_predecode_blocks,
            &tot_power_bit_mux_predecode_blocks,
            &tot_power_senseamp_mux_predecode_blocks, &tot_power_row_decoders,
            &tot_power_bit_mux_decoders, &tot_power_senseamp_mux_decoders);

    all_banks = area_all_banks(parameters->uca_banks, bank.height, bank.width, 
            number_bits_routed_to_bank, &length_htree_route_to_bank, number_mats_vertical_direction);
    //compute_delay_optimal_repeater(length_htree_route_to_bank, &sizing_repeater_htree_route_to_bank,
            //&number_repeaters_htree_route_to_bank);
    wst.wt = parameters->wire_inter_mats;
    wst.wire_length = length_htree_route_to_bank*1e-6;//m
    wst.nsense = 1;

    reset_powerDef(&power_routing_to_bank);
    if (wst.wire_length != 0) {
        calc_wire_stats2(parameters->wire_inter_mats, &wst);
        //	compute_delay_optimal_repeater(length_htree_route_to_bank, &sizing_repeater_htree_route_to_bank,
        //            &number_repeaters_htree_route_to_bank, &delay_route_to_bank, &power_routing_to_bank);
        delay_route_to_bank = wst.wire_pda.delay;
        power_routing_to_bank.readOp.dynamic = wst.wire_pda.power.dynamic;
        power_routing_to_bank.readOp.leakage = wst.wire_pda.power.leakage;
    }
    else {
        delay_route_to_bank = 0;
        power_routing_to_bank.readOp.dynamic = 0;
        power_routing_to_bank.readOp.leakage = 0;
    }


    /* if the bank size is too small then dont use low-swing wires */
    if((bank.height/2+bank.width)<1000 /* 1mm */) {
        if (parameters->wire_inter_mats == Low_swing) {
            parameters->wire_inter_mats = Global;
        }
    }

    inrisetime = 0;
    //delay_routing_to_bank(number_repeaters_htree_route_to_bank, length_htree_route_to_bank,
            //sizing_repeater_htree_route_to_bank, inrisetime, &outrisetime,  
            //&delay_route_to_bank, &power_routing_to_bank);
    inrisetime = 0;	

    delay_addr_datain_htree(horizontal_addr_din_htree_node, number_horizontal_htree_nodes, 
            number_vertical_htree_nodes, 
            inrisetime, &outrisetime, &delay_addr_din_horizontal_htree, parameters, &bank);
    inrisetime = 0;
    /* for low_swing, the entire bus delay is calculated in the vertical htree function */
        if (parameters->wire_inter_mats == Low_swing) delay_addr_din_horizontal_htree = 0;
#ifdef USE_GLOBAL
            delay_addr_din_horizontal_htree = 0;
#endif

    delay_addr_datain_htree(vertical_addr_din_htree_node, number_vertical_htree_nodes, 
            number_horizontal_htree_nodes, inrisetime, &outrisetime, 
            &delay_addr_din_vertical_htree, parameters, &bank);
    inrisetime = 0;

	if(HTREE_NODES_AT_MAT_INTERVALS){
		delay_addr_datain_htree_with_nodes_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval,
		number_horizontal_htree_nodes, inrisetime, &outrisetime, &delay_addr_din_horizontal_htree);
		inrisetime = 0;
		delay_addr_datain_htree_with_nodes_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval,
		number_vertical_htree_nodes, inrisetime, &outrisetime, &delay_addr_din_vertical_htree);
		inrisetime = 0;
	}

	delay_dout_horizontal_htree = delay_addr_din_horizontal_htree;

    delay_predecoder_block_driver(&row_predec_blk_driver_1, inrisetime, inrisetime, 
            &outrisetime_nand2_row_decode_path_1, &outrisetime_nand3_row_decode_path_1);
    delay_predecoder_block_driver(&row_predec_blk_driver_2, inrisetime, inrisetime,
            &outrisetime_nand2_row_decode_path_2, &outrisetime_nand3_row_decode_path_2);
    delay_predecoder_block_driver(&bit_mux_predec_blk_driver_1, inrisetime, inrisetime, 
            &outrisetime_nand2_bit_mux_decode_path_1, &outrisetime_nand3_bit_mux_decode_path_1);
    delay_predecoder_block_driver(&bit_mux_predec_blk_driver_2, inrisetime, inrisetime, 
            &outrisetime_nand2_bit_mux_decode_path_2, &outrisetime_nand3_bit_mux_decode_path_2);
    delay_predecoder_block_driver(&senseamp_mux_predec_blk_driver_1, inrisetime, inrisetime, 
            &outrisetime_nand2_senseamp_mux_decode_path_1, &outrisetime_nand3_senseamp_mux_decode_path_1);
    delay_predecoder_block_driver(&senseamp_mux_predec_blk_driver_2, inrisetime, inrisetime, 
            &outrisetime_nand2_senseamp_mux_decode_path_2, &outrisetime_nand3_senseamp_mux_decode_path_2);

    delay_predecoder_block(&row_predec_blk_1, outrisetime_nand2_row_decode_path_1,
            outrisetime_nand3_row_decode_path_1, &outrisetime_nand2_row_decode_path_1,
            &outrisetime_nand3_row_decode_path_1);
    delay_predecoder_block(&row_predec_blk_2, outrisetime_nand2_row_decode_path_2,
            outrisetime_nand3_row_decode_path_2, &outrisetime_nand2_row_decode_path_2,
            &outrisetime_nand3_row_decode_path_2);
    delay_predecoder_block(&bit_mux_predec_blk_1, outrisetime_nand2_bit_mux_decode_path_1,
            outrisetime_nand3_bit_mux_decode_path_1, &outrisetime_nand2_bit_mux_decode_path_1,
            &outrisetime_nand3_bit_mux_decode_path_1);
    delay_predecoder_block(&bit_mux_predec_blk_2, outrisetime_nand2_bit_mux_decode_path_2,
            outrisetime_nand3_bit_mux_decode_path_2, &outrisetime_nand2_bit_mux_decode_path_2,
            &outrisetime_nand3_bit_mux_decode_path_2);
    delay_predecoder_block(&senseamp_mux_predec_blk_1, outrisetime_nand2_senseamp_mux_decode_path_1,
            outrisetime_nand3_senseamp_mux_decode_path_1, &outrisetime_nand2_senseamp_mux_decode_path_1,
            &outrisetime_nand3_senseamp_mux_decode_path_1);
    delay_predecoder_block(&senseamp_mux_predec_blk_2, outrisetime_nand2_senseamp_mux_decode_path_2,
            outrisetime_nand3_senseamp_mux_decode_path_2, &outrisetime_nand2_senseamp_mux_decode_path_2,
            &outrisetime_nand3_senseamp_mux_decode_path_2);

    delay_before_row_decoder_nand2_path_1 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree +
        row_predec_blk_driver_1.delay_nand2_path + row_predec_blk_1.delay_nand2_path;
    delay_before_row_decoder_nand3_path_1 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        row_predec_blk_driver_1.delay_nand3_path + row_predec_blk_1.delay_nand3_path;
    delay_before_row_decoder_nand2_path_2 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree +
        row_predec_blk_driver_1.delay_nand2_path + row_predec_blk_1.delay_nand2_path;
    delay_before_row_decoder_nand3_path_2 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        row_predec_blk_driver_1.delay_nand3_path + row_predec_blk_1.delay_nand3_path;

    max_delay_before_row_decoder = delay_before_row_decoder_nand2_path_1;
    row_dec_inrisetime = outrisetime_nand2_row_decode_path_1;
    if(delay_before_row_decoder_nand3_path_1 > max_delay_before_row_decoder){
        max_delay_before_row_decoder = delay_before_row_decoder_nand3_path_1;
        row_dec_inrisetime = outrisetime_nand3_row_decode_path_1;
    }
    if(delay_before_row_decoder_nand2_path_2 > max_delay_before_row_decoder){
        max_delay_before_row_decoder = delay_before_row_decoder_nand2_path_2;
        row_dec_inrisetime = outrisetime_nand2_row_decode_path_2;
    }
    if(delay_before_row_decoder_nand3_path_2 > max_delay_before_row_decoder){
        max_delay_before_row_decoder = delay_before_row_decoder_nand3_path_2;
        row_dec_inrisetime = outrisetime_nand3_row_decode_path_2;
    }

    delay_before_bit_mux_decoder_nand2_path_1 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        bit_mux_predec_blk_driver_1.delay_nand2_path + bit_mux_predec_blk_1.delay_nand2_path;
    delay_before_bit_mux_decoder_nand3_path_1 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree +
        bit_mux_predec_blk_driver_1.delay_nand3_path + bit_mux_predec_blk_1.delay_nand3_path;
    delay_before_bit_mux_decoder_nand2_path_2 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree +	
        bit_mux_predec_blk_driver_1.delay_nand2_path + bit_mux_predec_blk_1.delay_nand2_path;
    delay_before_bit_mux_decoder_nand3_path_2 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        bit_mux_predec_blk_driver_1.delay_nand3_path + bit_mux_predec_blk_1.delay_nand3_path;

    max_delay_before_bit_mux_decoder = delay_before_bit_mux_decoder_nand2_path_1;
    bit_mux_dec_inrisetime = outrisetime_nand2_bit_mux_decode_path_1;
    if(delay_before_bit_mux_decoder_nand3_path_1 > max_delay_before_bit_mux_decoder){
        max_delay_before_bit_mux_decoder = delay_before_bit_mux_decoder_nand3_path_1;
        bit_mux_dec_inrisetime = outrisetime_nand3_bit_mux_decode_path_1;
    }
    if(delay_before_bit_mux_decoder_nand2_path_2 > max_delay_before_bit_mux_decoder){
        max_delay_before_bit_mux_decoder = delay_before_bit_mux_decoder_nand2_path_2;
        bit_mux_dec_inrisetime = outrisetime_nand2_bit_mux_decode_path_2;
    }
    if(delay_before_bit_mux_decoder_nand3_path_2 > max_delay_before_bit_mux_decoder){
        max_delay_before_bit_mux_decoder = delay_before_bit_mux_decoder_nand3_path_2;
        bit_mux_dec_inrisetime = outrisetime_nand3_bit_mux_decode_path_2;
    }

    delay_before_senseamp_mux_decoder_nand2_path_1 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        senseamp_mux_predec_blk_driver_1.delay_nand2_path + senseamp_mux_predec_blk_1.delay_nand2_path;
    delay_before_senseamp_mux_decoder_nand3_path_1 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        senseamp_mux_predec_blk_driver_1.delay_nand3_path + senseamp_mux_predec_blk_1.delay_nand3_path;
    delay_before_senseamp_mux_decoder_nand2_path_2 = delay_route_to_bank +
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        senseamp_mux_predec_blk_driver_1.delay_nand2_path + senseamp_mux_predec_blk_1.delay_nand2_path;
    delay_before_senseamp_mux_decoder_nand3_path_2 = delay_route_to_bank + 
        delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree + 
        senseamp_mux_predec_blk_driver_1.delay_nand3_path + senseamp_mux_predec_blk_1.delay_nand3_path;

    max_delay_before_senseamp_mux_decoder = delay_before_senseamp_mux_decoder_nand2_path_1;
    senseamp_mux_dec_inrisetime = outrisetime_nand2_senseamp_mux_decode_path_1;
    if(delay_before_senseamp_mux_decoder_nand3_path_1 > max_delay_before_senseamp_mux_decoder){
        max_delay_before_senseamp_mux_decoder = delay_before_senseamp_mux_decoder_nand3_path_1;
        senseamp_mux_dec_inrisetime = outrisetime_nand3_senseamp_mux_decode_path_1;
    }
    if(delay_before_senseamp_mux_decoder_nand2_path_2 > max_delay_before_senseamp_mux_decoder){
        max_delay_before_senseamp_mux_decoder = delay_before_senseamp_mux_decoder_nand2_path_2;
        senseamp_mux_dec_inrisetime = outrisetime_nand2_senseamp_mux_decode_path_2;
    }
    if(delay_before_senseamp_mux_decoder_nand3_path_2 > max_delay_before_senseamp_mux_decoder){
        max_delay_before_senseamp_mux_decoder = delay_before_senseamp_mux_decoder_nand3_path_2;
        senseamp_mux_dec_inrisetime = outrisetime_nand3_senseamp_mux_decode_path_2;
    }

	is_wordline_transistor = 1;
    delay_decoder(&row_dec, row_dec_inrisetime, &row_dec_outrisetime);
	is_wordline_transistor = 0;
    delay_decoder(&bit_mux_dec, bit_mux_dec_inrisetime, &bit_mux_dec_outrisetime);
    delay_decoder(&senseamp_mux_dec, senseamp_mux_dec_inrisetime, &senseamp_mux_dec_outrisetime);

    Tpre = max_delay_before_row_decoder + row_dec.delay;
    reset_powerDef(&bitline_data_power);
	leakage_power_cc_inverters_sram_cell = 0;
	leakage_power_access_transistors_read_write_port_sram_cell = 0;
	leakage_power_read_only_port_sram_cell = 0;
  
    if (is_dram) {
    bitline_data = bitline_delay(number_rows_subarray, number_cols_subarray, number_subarrays, 
            row_dec_outrisetime, &outrisetime, &bitline_data_power, Tpre, deg_bitline_muxing,
            number_activated_mats_horizontal_direction, &writeback_delay,
			parameters->rw_ports, parameters->excl_read_ports, parameters->excl_write_ports);
    }
    else {
    bitline_data = bitline_delay_sram(number_rows_subarray, number_cols_subarray, number_subarrays, 
            row_dec_outrisetime, &outrisetime, &bitline_data_power, Tpre, deg_bitline_muxing,
            number_activated_mats_horizontal_direction, &writeback_delay,
			parameters->rw_ports, parameters->excl_read_ports, parameters->excl_write_ports);
//    bitline_data = bitline_delay(number_rows_subarray, number_cols_subarray, number_subarrays, 
//          row_dec_outrisetime, &outrisetime, &bitline_data_power, Tpre, deg_bitline_muxing,
//          number_activated_mats_horizontal_direction, &writeback_delay,
//			parameters->rw_ports, parameters->excl_read_ports, parameters->excl_write_ports);
    }

    wordline_data = row_dec.delay;
    inrisetime = outrisetime;

    reset_powerDef(&sense_amp_data_power);
	delay_sense_amplifier(number_cols_subarray, parameters->rw_ports, parameters->excl_read_ports, inrisetime, &outrisetime, &sense_amp_data_power, 
            &subarray_output_htree_node, deg_bitline_muxing, number_mats, number_activated_mats_horizontal_direction, 
            &delay_sense_amp);
    inrisetime = outrisetime;
    reset_powerDef(&power_output_drivers_at_subarray);
    delay_output_driver_at_subarray(Ndsam, inrisetime, &outrisetime, &power_output_drivers_at_subarray, 
            &subarray_output_htree_node, number_mats, &delay_subarray_output_driver,
            &delay_final_stage_subarray_output_driver);
    delay_dataout_vertical_htree(dout_htree_node, number_vertical_htree_nodes, number_vertical_htree_nodes,  inrisetime, &outrisetime, &delay_dout_vertical_htree,
            &power_dout_htree, parameters, &bank);
    inrisetime = outrisetime;

	if(HTREE_NODES_AT_MAT_INTERVALS){
		delay_datout_vertical_htree_with_nodes_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval,
			number_vertical_htree_nodes, inrisetime, &outrisetime, &delay_dout_vertical_htree, &power_dout_htree);
		inrisetime = 0;
	}

    delay_data_access_row_path =  max_delay_before_row_decoder + row_dec.delay + bitline_data + 
        delay_sense_amp +  delay_subarray_output_driver + delay_dout_vertical_htree + 
        delay_dout_horizontal_htree + delay_route_to_bank;
    delay_data_access_col_path = max_delay_before_bit_mux_decoder + bit_mux_dec.delay +	delay_sense_amp +  
        delay_subarray_output_driver + delay_dout_vertical_htree + delay_dout_horizontal_htree + 
        delay_route_to_bank;
    delay_data_access_out_drv_path = max_delay_before_senseamp_mux_decoder + senseamp_mux_dec.delay + 
        delay_final_stage_subarray_output_driver + 	delay_dout_vertical_htree + 	delay_dout_horizontal_htree + 
        delay_route_to_bank;
    temp_delay_data_access_path = MAX(delay_data_access_row_path, delay_data_access_col_path);
    delay_data_access_path = MAX(temp_delay_data_access_path, delay_data_access_out_drv_path);
    delay_before_subarray_output_driver = delay_data_access_path - 
        (delay_final_stage_subarray_output_driver + delay_dout_vertical_htree + delay_dout_horizontal_htree + 
         delay_route_to_bank);
    delay_from_subarray_output_driver_to_output = delay_final_stage_subarray_output_driver + delay_dout_vertical_htree +
        delay_dout_horizontal_htree + delay_route_to_bank;

    access_time = delay_data_access_path;
    reset_powerDef(&comparator_power);
    if(is_tag){
        delay_comparator(tagbits, parameters->tag_associativity, inrisetime, &outrisetime,
                &comparator_delay, &comparator_power);
        access_time += comparator_delay;
    }
    else{
        comparator_delay = 0;
    }

    delay_driver(&bitline_precharge_eq_driver, 0, &dummy_precharge_outrisetime);
    if(is_dram){
		k = row_dec.number_gates - 1;
        c_load = row_dec.c_load_decoder_output;
        rd = transreson(row_dec.width_decoder_n[k], NCH, 1);
        c_intrinsic = draincap(row_dec.width_decoder_p[k], PCH, 1, 1, 4 * height_cell) + 
            draincap(row_dec.width_decoder_n[k], NCH, 1, 1, 4 * height_cell);
        tf = rd * (c_intrinsic + c_load) + row_dec.r_wire_decoder_output * c_load / 2;
        wordline_reset_delay = 0.69 * tf;
		r_bitline_precharge = transreson(width_pmos_bitline_precharge, PCH, 1);
		r_bitline = number_rows_subarray * Rbitmetal;
		bitline_restore_delay = bitline_precharge_eq_driver.delay + 
			2.3 * (r_bitline_precharge * Cbitline + r_bitline * Cbitline / 2);
        cycle_time = wordline_data + bitline_data + delay_sense_amp + writeback_delay +
			wordline_reset_delay + bitline_restore_delay;
		dram_array_availability = (1 - number_rows_subarray * cycle_time / dram_refresh_period) * 100;
    }
    else{
        k = row_dec.number_gates - 1;
        c_load = row_dec.c_load_decoder_output;
        rd = transreson(row_dec.width_decoder_n[k], NCH, 1);
        c_intrinsic = draincap(row_dec.width_decoder_p[k], PCH, 1, 1, 4 * height_cell) + 
            draincap(row_dec.width_decoder_n[k], NCH, 1, 1, 4 * height_cell);
        tf = rd * (c_intrinsic + c_load) + row_dec.r_wire_decoder_output * c_load / 2;
        wordline_reset_delay = horowitz(0, tf, 0.5, 0.5, RISE);
		r_bitline_precharge = transreson(width_pmos_bitline_precharge, PCH, 1);
		r_bitline = number_rows_subarray * Rbitmetal;
		bitline_restore_delay = bitline_precharge_eq_driver.delay + 
			log((Vbitpre_sram - 0.1 * Vbitsense) / (Vbitpre_sram - Vbitsense)) * (r_bitline_precharge * Cbitline + 
			r_bitline * Cbitline / 2);
        cycle_time = wordline_data + bitline_data + wordline_reset_delay + bitline_restore_delay;
		dram_array_availability = 0;
    }	
	temp = MAX(max_delay_before_row_decoder, max_delay_before_bit_mux_decoder);
	temp = MAX(temp, max_delay_before_senseamp_mux_decoder);
	temp = MAX(temp, delay_subarray_output_driver + delay_dout_vertical_htree + delay_dout_horizontal_htree +
		delay_route_to_bank);
	multisubbank_interleave_cycle_time = MIN(temp, cycle_time);

    tot_power_routing_to_bank.readOp.dynamic = power_routing_to_bank.readOp.dynamic *
        (number_addr_bits_routed_to_bank + number_data_bits_routed_to_bank / 2 +
        number_comparator_bits_routed_to_bank);
    tot_power.readOp.dynamic += tot_power_routing_to_bank.readOp.dynamic;

	if(HTREE_NODES_AT_MAT_INTERVALS){
		calculate_power_htree_with_nodes_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval);
		calculate_power_htree_with_nodes_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval);
	}
    htree_seg_multiplier = 1;
    /* htree power calculation happens here */
	for(k = 0; k < number_horizontal_htree_nodes; ++k){
        /* for use_global and low_swing, power calculation is taken care at vertical htree */
        if (parameters->wire_inter_mats == Low_swing) break;
#ifdef USE_GLOBAL
            break;
#endif
		power_addr_datain_htree_node(&horizontal_addr_din_htree_node[k], 
			&power_addr_datain_horizontal_htree_node, parameters, &bank, number_horizontal_htree_nodes, number_horizontal_htree_nodes);
		if(HTREE_NODES_AT_MAT_INTERVALS) {
			number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_horizontal_htree_nodes - k) + 0.5)) / 2, 1);
			power_addr_datain_htree_node_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment, &power_addr_datain_horizontal_htree_node);
		}
        tot_power_addr_horizontal_htree.readOp.dynamic +=
            (number_addr_bits_mat + number_way_select_signals_mat) * 
			power_addr_datain_horizontal_htree_node.readOp.dynamic * 
            htree_seg_multiplier;
        tot_power_datain_horizontal_htree.readOp.dynamic += number_datain_bits_mat * 
			number_mats_horizontal_direction *	
            power_addr_datain_horizontal_htree_node.readOp.dynamic;
        tot_power_dataout_horizontal_htree.readOp.dynamic += number_dataout_bits_mat * 
            number_mats_horizontal_direction *	
            power_addr_datain_horizontal_htree_node.readOp.dynamic;
        htree_seg_multiplier *= 2;
    }

	tot_power.readOp.dynamic += tot_power_addr_horizontal_htree.readOp.dynamic +
         + tot_power_dataout_horizontal_htree.readOp.dynamic;	

	//Add energy consumed in address/datain vertical htree
	htree_seg_multiplier = 1;
	for(k = 0; k < number_vertical_htree_nodes; ++k){
		power_addr_datain_htree_node(&vertical_addr_din_htree_node[k], 
                    &power_addr_datain_vertical_htree_node, parameters, &bank, 
                    number_vertical_htree_nodes, number_horizontal_htree_nodes);
		if(HTREE_NODES_AT_MAT_INTERVALS){
			number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_vertical_htree_nodes - k) + 0.5)) / 2, 1);
			power_addr_datain_htree_node_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment, &power_addr_datain_vertical_htree_node);
		}
        tot_power_addr_vertical_htree.readOp.dynamic +=
            (number_addr_bits_mat + number_way_select_signals_mat) * 
			power_addr_datain_vertical_htree_node.readOp.dynamic * htree_seg_multiplier * 
			number_mats_horizontal_direction;
        tot_power_datain_vertical_htree.readOp.dynamic += number_datain_bits_mat * 
			power_addr_datain_vertical_htree_node.readOp.dynamic * htree_seg_multiplier * 
			number_mats_horizontal_direction;
//        tot_power_datain_vertical_htree.readOp.dynamic;
//        tot_power_addr_vertical_htree.readOp.dynamic;
		if(BROADCAST_ADDR_DATAIN_OVER_VERTICAL_HTREES){
			htree_seg_multiplier *= 2;
		}
		else{
			htree_seg_multiplier = 1;
		}
        if (parameters->wire_inter_mats == Low_swing) break;
#ifdef USE_GLOBAL
            break;
#endif
    }

	if((number_mats_vertical_direction > 1)&&(BROADCAST_ADDR_DATAIN_OVER_VERTICAL_HTREES)){
		tot_power_addr_vertical_htree.readOp.dynamic *= 2; 
		tot_power_datain_vertical_htree.readOp.dynamic *= 2;
	}
	
	tot_power.readOp.dynamic += tot_power_addr_vertical_htree.readOp.dynamic;

    //Add power consumed in predecoder drivers
    number_addr_bits_nand2_row_decode_path_1 = 
        row_predec_blk_driver_1.number_parallel_instances_driving_1_nand2_load +
        row_predec_blk_driver_1.number_parallel_instances_driving_2_nand2_load +
        row_predec_blk_driver_1.number_parallel_instances_driving_4_nand2_load;
    number_addr_bits_nand3_row_decode_path_1 = 
        row_predec_blk_driver_1.number_parallel_instances_driving_2_nand3_load +
        row_predec_blk_driver_1.number_parallel_instances_driving_8_nand3_load;
    number_addr_bits_nand2_row_decode_path_2 = 
        row_predec_blk_driver_2.number_parallel_instances_driving_1_nand2_load +
        row_predec_blk_driver_2.number_parallel_instances_driving_2_nand2_load +
        row_predec_blk_driver_2.number_parallel_instances_driving_4_nand2_load;
    number_addr_bits_nand3_row_decode_path_2 = 
        row_predec_blk_driver_2.number_parallel_instances_driving_2_nand3_load +
        row_predec_blk_driver_2.number_parallel_instances_driving_8_nand3_load;
    number_addr_bits_nand2_bit_mux_decode_path_1 = 
        bit_mux_predec_blk_driver_1.number_parallel_instances_driving_1_nand2_load +
        bit_mux_predec_blk_driver_1.number_parallel_instances_driving_2_nand2_load +
        bit_mux_predec_blk_driver_1.number_parallel_instances_driving_4_nand2_load;
    number_addr_bits_nand3_bit_mux_decode_path_1 = 
        bit_mux_predec_blk_driver_1.number_parallel_instances_driving_2_nand3_load +
        bit_mux_predec_blk_driver_1.number_parallel_instances_driving_8_nand3_load;
    number_addr_bits_nand2_bit_mux_decode_path_2 = 
        bit_mux_predec_blk_driver_2.number_parallel_instances_driving_1_nand2_load +
        bit_mux_predec_blk_driver_2.number_parallel_instances_driving_2_nand2_load +
        bit_mux_predec_blk_driver_2.number_parallel_instances_driving_4_nand2_load;
    number_addr_bits_nand3_bit_mux_decode_path_2 = 
        bit_mux_predec_blk_driver_2.number_parallel_instances_driving_2_nand3_load +
        bit_mux_predec_blk_driver_2.number_parallel_instances_driving_8_nand3_load;
    number_addr_bits_nand2_senseamp_mux_decode_path_1 = 
        senseamp_mux_predec_blk_driver_1.number_parallel_instances_driving_1_nand2_load +
        senseamp_mux_predec_blk_driver_1.number_parallel_instances_driving_2_nand2_load +
        senseamp_mux_predec_blk_driver_1.number_parallel_instances_driving_4_nand2_load;
    number_addr_bits_nand3_senseamp_mux_decode_path_1 = 
        senseamp_mux_predec_blk_driver_1.number_parallel_instances_driving_2_nand3_load +
        senseamp_mux_predec_blk_driver_1.number_parallel_instances_driving_8_nand3_load;
    number_addr_bits_nand2_senseamp_mux_decode_path_2 = 
        senseamp_mux_predec_blk_driver_2.number_parallel_instances_driving_1_nand2_load +
        senseamp_mux_predec_blk_driver_2.number_parallel_instances_driving_2_nand2_load +
        senseamp_mux_predec_blk_driver_2.number_parallel_instances_driving_4_nand2_load;
    number_addr_bits_nand3_senseamp_mux_decode_path_2 = 
        senseamp_mux_predec_blk_driver_2.number_parallel_instances_driving_2_nand3_load +
        senseamp_mux_predec_blk_driver_2.number_parallel_instances_driving_8_nand3_load;

    tot_power_row_predecode_block_drivers.readOp.dynamic = 
    (number_addr_bits_nand2_row_decode_path_1 *	row_predec_blk_driver_1.power_nand2_path.readOp.dynamic + 
    number_addr_bits_nand3_row_decode_path_1 * row_predec_blk_driver_1.power_nand3_path.readOp.dynamic +
    number_addr_bits_nand2_row_decode_path_2 *	row_predec_blk_driver_2.power_nand2_path.readOp.dynamic + 
    number_addr_bits_nand3_row_decode_path_2 * row_predec_blk_driver_2.power_nand3_path.readOp.dynamic) *
        number_activated_mats_horizontal_direction;

    tot_power_bit_mux_predecode_block_drivers.readOp.dynamic  = 
        (number_addr_bits_nand2_bit_mux_decode_path_1 *	bit_mux_predec_blk_driver_1.power_nand2_path.readOp.dynamic + 
         number_addr_bits_nand3_bit_mux_decode_path_1 * bit_mux_predec_blk_driver_1.power_nand3_path.readOp.dynamic +
         number_addr_bits_nand2_bit_mux_decode_path_2 *	bit_mux_predec_blk_driver_2.power_nand2_path.readOp.dynamic + 
         number_addr_bits_nand3_bit_mux_decode_path_2 * bit_mux_predec_blk_driver_2.power_nand3_path.readOp.dynamic) *
        number_activated_mats_horizontal_direction;

   
    tot_power_senseamp_mux_predecode_block_drivers.readOp.dynamic = 
		(number_addr_bits_nand2_senseamp_mux_decode_path_1 * senseamp_mux_predec_blk_driver_1.power_nand2_path.readOp.dynamic + 
		number_addr_bits_nand3_senseamp_mux_decode_path_1 * senseamp_mux_predec_blk_driver_1.power_nand3_path.readOp.dynamic + 
		number_addr_bits_nand2_senseamp_mux_decode_path_2 * senseamp_mux_predec_blk_driver_2.power_nand2_path.readOp.dynamic + 
		number_addr_bits_nand3_senseamp_mux_decode_path_2 * senseamp_mux_predec_blk_driver_2.power_nand3_path.readOp.dynamic) *
		number_activated_mats_horizontal_direction;

    tot_power.readOp.dynamic += tot_power_row_predecode_block_drivers.readOp.dynamic +
        tot_power_bit_mux_predecode_block_drivers.readOp.dynamic +
        tot_power_senseamp_mux_predecode_block_drivers.readOp.dynamic;

    //Add power consumed in predecode blocks
    tot_power_row_predecode_blocks.readOp.dynamic = 
        (row_predec_blk_1.power_nand2_path.readOp.dynamic + row_predec_blk_1.power_nand3_path.readOp.dynamic +
         row_predec_blk_1.power_second_level.readOp.dynamic + row_predec_blk_2.power_nand2_path.readOp.dynamic + 
         row_predec_blk_2.power_nand3_path.readOp.dynamic +	row_predec_blk_2.power_second_level.readOp.dynamic) *
        number_activated_mats_horizontal_direction;

    tot_power_bit_mux_predecode_blocks.readOp.dynamic =
        (bit_mux_predec_blk_1.power_nand2_path.readOp.dynamic +	bit_mux_predec_blk_1.power_nand3_path.readOp.dynamic +
         bit_mux_predec_blk_1.power_second_level.readOp.dynamic + bit_mux_predec_blk_2.power_nand2_path.readOp.dynamic +
         bit_mux_predec_blk_2.power_nand3_path.readOp.dynamic +	bit_mux_predec_blk_2.power_second_level.readOp.dynamic) *
        number_activated_mats_horizontal_direction;

    tot_power_senseamp_mux_predecode_blocks.readOp.dynamic =
            (senseamp_mux_predec_blk_1.power_nand2_path.readOp.dynamic + senseamp_mux_predec_blk_1.power_nand3_path.readOp.dynamic +
             senseamp_mux_predec_blk_1.power_second_level.readOp.dynamic + senseamp_mux_predec_blk_2.power_nand2_path.readOp.dynamic + 
             senseamp_mux_predec_blk_2.power_nand3_path.readOp.dynamic +	senseamp_mux_predec_blk_2.power_second_level.readOp.dynamic) *
            number_activated_mats_horizontal_direction;
   
    tot_power.readOp.dynamic += tot_power_row_predecode_blocks.readOp.dynamic +
        tot_power_bit_mux_predecode_blocks.readOp.dynamic +	tot_power_senseamp_mux_predecode_blocks.readOp.dynamic;

    //Add power consumed in decoders
    tot_power_row_decoders.readOp.dynamic = row_dec.power.readOp.dynamic * 	4 * 
        number_activated_mats_horizontal_direction;

    //If DRAM, add contribution of power spent in row predecoder drivers, blocks and decoders to refresh power
    if(is_dram){
        refresh_power += (tot_power_row_predecode_block_drivers.readOp.dynamic  +
		tot_power_row_predecode_blocks.readOp.dynamic + row_dec.power.readOp.dynamic) * 
		number_rows_subarray * number_subarrays /  dram_refresh_period;
    }

    tot_power_bit_mux_decoders.readOp.dynamic = bit_mux_dec.power.readOp.dynamic * 
		number_activated_mats_horizontal_direction;

   tot_power_senseamp_mux_decoders.readOp.dynamic = senseamp_mux_dec.power.readOp.dynamic * 
	   number_activated_mats_horizontal_direction;

   tot_power.readOp.dynamic += tot_power_row_decoders.readOp.dynamic +	
	   tot_power_bit_mux_decoders.readOp.dynamic + tot_power_senseamp_mux_decoders.readOp.dynamic;

    //Add power consumed in bitlines
    tot_power_bitlines.readOp.dynamic = bitline_data_power.readOp.dynamic;
    tot_power.readOp.dynamic += tot_power_bitlines.readOp.dynamic;

	//Add power consumed in the precharge and equalization driver
	tot_power_bitlines_precharge_eq_driver.readOp.dynamic = bitline_precharge_eq_driver.power.readOp.dynamic *
		4 * number_activated_mats_horizontal_direction;
	tot_power.readOp.dynamic += tot_power_bitlines_precharge_eq_driver.readOp.dynamic;

    //Add power consumed in sense amplifiers
    tot_power_sense_amps.readOp.dynamic = sense_amp_data_power.readOp.dynamic;
    tot_power.readOp.dynamic += tot_power_sense_amps.readOp.dynamic;

    //Add power consumed in subarray output driver circuitry
    //Actually the below equation hold good only for the final stage of the subarray output driver - only
    //number_dataout_bits_mat final stages are activated. However, the earlier stages for all
    //number_sense_amps_or_output_drivers_subarray output drivers are active. FIX this.
    tot_power_subarray_output_drivers.readOp.dynamic = power_output_drivers_at_subarray.readOp.dynamic * 
        number_dataout_bits_mat * number_activated_mats_horizontal_direction;
    tot_power.readOp.dynamic += tot_power_subarray_output_drivers.readOp.dynamic;

    //Add power consumed in dataout tree
    tot_power_dataout_vertical_htree.readOp.dynamic = number_dataout_bits_mat * 
        number_activated_mats_horizontal_direction * power_dout_htree.readOp.dynamic;
    tot_power.readOp.dynamic += tot_power_dataout_vertical_htree.readOp.dynamic;

    tot_power_comparators.readOp.dynamic = comparator_power.readOp.dynamic;
    tot_power.readOp.dynamic += tot_power_comparators.readOp.dynamic;

    //Leakage power 
    /* NOTE: the same function is redundantly called to calculate leakage.TODO - merge
     it with the original function */
    tot_power_routing_to_bank.readOp.leakage += power_routing_to_bank.readOp.leakage *
		number_bits_routed_to_bank * (parameters->rw_ports + parameters->excl_read_ports +
		parameters->excl_write_ports);
    tot_power.readOp.leakage += tot_power_routing_to_bank.readOp.leakage;
    
	htree_seg_multiplier = 1;
	for(k = 0; k < number_horizontal_htree_nodes; ++k){
        /* for use_global and low_swing, power calculation is taken care at vertical htree */
        if (parameters->wire_inter_mats == Low_swing) break;
#ifdef USE_GLOBAL
            break;
#endif
		power_addr_datain_htree_node(&horizontal_addr_din_htree_node[k],
			&power_addr_datain_horizontal_htree_node, parameters, &bank, number_horizontal_htree_nodes, number_horizontal_htree_nodes);
		if(HTREE_NODES_AT_MAT_INTERVALS){
			number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_horizontal_htree_nodes - k) + 0.5)) / 2, 1);
			power_addr_datain_htree_node_at_mat_interval(&horizontal_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment, 
                &power_addr_datain_horizontal_htree_node 
                );
		}
        tot_power_addr_horizontal_htree.readOp.leakage +=
            (number_addr_bits_mat + number_way_select_signals_mat) * 
			power_addr_datain_horizontal_htree_node.readOp.leakage * htree_seg_multiplier *
			(parameters->rw_ports + parameters->excl_read_ports +
			parameters->excl_write_ports);
        tot_power_datain_horizontal_htree.readOp.leakage += number_datain_bits_mat * 
			number_mats_horizontal_direction * 
			power_addr_datain_horizontal_htree_node.readOp.leakage * (parameters->rw_ports + 
			parameters->excl_write_ports);
        tot_power_dataout_horizontal_htree.readOp.leakage += number_dataout_bits_mat * 
            number_mats_horizontal_direction * power_addr_datain_horizontal_htree_node.readOp.leakage *
			(parameters->rw_ports + parameters->excl_read_ports);
        htree_seg_multiplier *= 2;
    }

    tot_power.readOp.leakage += tot_power_addr_horizontal_htree.readOp.leakage +
		tot_power_datain_horizontal_htree.readOp.leakage +
		tot_power_dataout_horizontal_htree.readOp.leakage;

	htree_seg_multiplier = 1;
	for(k = 0; k < number_vertical_htree_nodes; ++k){
		power_addr_datain_htree_node(&vertical_addr_din_htree_node[k], &power_addr_datain_vertical_htree_node, parameters, &bank, number_vertical_htree_nodes, number_horizontal_htree_nodes);

		if(HTREE_NODES_AT_MAT_INTERVALS){
			number_mats_to_cover_in_this_htree_segment = 
				MAX(((int)(pow(2.0, number_vertical_htree_nodes - k) + 0.5)) / 2, 1);
			power_addr_datain_htree_node_at_mat_interval(&vertical_addr_datain_htree_at_mat_interval, 
				number_mats_to_cover_in_this_htree_segment, &power_addr_datain_vertical_htree_node);
		}
        tot_power_addr_vertical_htree.readOp.leakage +=
            (number_addr_bits_mat + number_way_select_signals_mat) * 
			power_addr_datain_vertical_htree_node.readOp.leakage * htree_seg_multiplier * 
			number_mats_horizontal_direction * (parameters->rw_ports + 
			parameters->excl_read_ports + parameters->excl_write_ports);
        tot_power_datain_vertical_htree.readOp.leakage += number_datain_bits_mat * 
			power_addr_datain_vertical_htree_node.readOp.leakage * htree_seg_multiplier * 
			number_mats_horizontal_direction * (parameters->rw_ports + 
			parameters->excl_write_ports);
		htree_seg_multiplier *= 2;
        if (parameters->wire_inter_mats == Low_swing) break;
#ifdef USE_GLOBAL
            break;
#endif
    }

	if(number_mats_vertical_direction > 1){
		tot_power_addr_vertical_htree.readOp.leakage *= 2; 
		tot_power_datain_vertical_htree.readOp.leakage *= 2;
	}
	
	tot_power.readOp.leakage += tot_power_addr_vertical_htree.readOp.leakage +
		tot_power_datain_vertical_htree.readOp.leakage;

	if(number_mats_vertical_direction > 1){
		tot_power_dataout_vertical_htree.readOp.leakage = number_dataout_bits_mat * 
			number_mats_horizontal_direction *	2 * power_dout_htree.readOp.leakage *
			(parameters->rw_ports + parameters->excl_read_ports);
	}
	else{//number_mats_vertical_direction = 1
		tot_power_dataout_vertical_htree.readOp.leakage = number_dataout_bits_mat * 
			number_mats_horizontal_direction *	1 * power_dout_htree.readOp.leakage * 
			(parameters->rw_ports + parameters->excl_read_ports);
	}

    tot_power.readOp.leakage += tot_power_dataout_vertical_htree.readOp.leakage;

    tot_power_bitlines.readOp.leakage = bitline_data_power.readOp.leakage *
		number_rows_subarray * number_cols_subarray * number_subarrays; //This is actually
	//leakage in the memory cells.
    tot_power.readOp.leakage += tot_power_bitlines.readOp.leakage;

	tot_power_bitlines_precharge_eq_driver.readOp.leakage = 
		bitline_precharge_eq_driver.power.readOp.leakage * number_subarrays;
	tot_power.readOp.leakage += tot_power_bitlines_precharge_eq_driver.readOp.leakage;

    tot_power_sense_amps.readOp.leakage = sense_amp_data_power.readOp.leakage * 
		(parameters->rw_ports + parameters->excl_read_ports);
    tot_power.readOp.leakage += tot_power_sense_amps.readOp.leakage;

    tot_power_subarray_output_drivers.readOp.leakage = power_output_drivers_at_subarray.readOp.leakage *
        number_output_drivers_subarray * number_subarrays * 
		(parameters->rw_ports + parameters->excl_read_ports);
    tot_power.readOp.leakage += tot_power_subarray_output_drivers.readOp.leakage;

    tot_power_comparators.readOp.leakage = comparator_power.readOp.leakage * 
		(parameters->rw_ports + parameters->excl_read_ports);
    tot_power.readOp.leakage += tot_power_comparators.readOp.leakage;

    //If DRAM, add refresh power to total leakage
    if(is_dram){
        tot_power.readOp.leakage += refresh_power;
    }
    if(flag_results_populate){
        ptr_results->Ndwl = Ndwl;
        ptr_results->Ndbl = Ndbl;
        ptr_results->Nspd = Nspd;
        ptr_results->deg_bitline_muxing = deg_bitline_muxing;
        ptr_results->Ndsam = Ndsam;
        ptr_results->number_activated_mats_horizontal_direction = number_activated_mats_horizontal_direction;
        ptr_results->delay_route_to_bank = delay_route_to_bank;
        ptr_results->delay_addr_din_horizontal_htree = delay_addr_din_horizontal_htree;
        ptr_results->delay_addr_din_vertical_htree = delay_addr_din_vertical_htree;
        ptr_results->delay_row_predecode_driver_and_block = max_delay_before_row_decoder -
            (delay_route_to_bank + delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree);
        ptr_results->delay_row_decoder = row_dec.delay;
        ptr_results->delay_bitlines = bitline_data;
        ptr_results->delay_sense_amp = delay_sense_amp;
        ptr_results->delay_subarray_output_driver = delay_subarray_output_driver;
        ptr_results->delay_bit_mux_predecode_driver_and_block = max_delay_before_bit_mux_decoder -
            (delay_route_to_bank + delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree);
        ptr_results->delay_bit_mux_decoder = bit_mux_dec.delay;
        ptr_results->delay_senseamp_mux_predecode_driver_and_block = max_delay_before_senseamp_mux_decoder -
            (delay_route_to_bank + delay_addr_din_horizontal_htree + delay_addr_din_vertical_htree);
        ptr_results->delay_senseamp_mux_decoder = senseamp_mux_dec.delay;
        ptr_results->delay_dout_vertical_htree = delay_dout_vertical_htree;
        ptr_results->delay_dout_horizontal_htree = delay_dout_horizontal_htree;
        ptr_results->delay_comparator = comparator_delay;
        ptr_results->access_time = access_time;
        ptr_results->cycle_time = cycle_time;
        ptr_results->multisubbank_interleave_cycle_time = multisubbank_interleave_cycle_time;
		ptr_results->dram_refresh_period = dram_refresh_period;
		ptr_results->dram_array_availability = dram_array_availability;
        copy_powerDef(&ptr_results->power_routing_to_bank, tot_power_routing_to_bank);
        copy_powerDef(&ptr_results->power_addr_horizontal_htree, 
                tot_power_addr_horizontal_htree);
        copy_powerDef(&ptr_results->power_datain_horizontal_htree, 
                tot_power_datain_horizontal_htree);
        copy_powerDef(&ptr_results->power_dataout_horizontal_htree, 
                tot_power_dataout_horizontal_htree);
        copy_powerDef(&ptr_results->power_addr_vertical_htree, tot_power_addr_vertical_htree);
		copy_powerDef(&ptr_results->power_datain_vertical_htree, tot_power_datain_vertical_htree);
        copy_powerDef(&ptr_results->power_row_predecoder_drivers, tot_power_row_predecode_block_drivers);
        copy_powerDef(&ptr_results->power_row_predecoder_blocks, tot_power_row_predecode_blocks);
        copy_powerDef(&ptr_results->power_row_decoders, tot_power_row_decoders);
        copy_powerDef(&ptr_results->power_bit_mux_predecoder_drivers, tot_power_bit_mux_predecode_block_drivers);
        copy_powerDef(&ptr_results->power_bit_mux_predecoder_blocks, tot_power_bit_mux_predecode_blocks);
        copy_powerDef(&ptr_results->power_bit_mux_decoders, tot_power_bit_mux_decoders);
        copy_powerDef(&ptr_results->power_senseamp_mux_predecoder_drivers, tot_power_senseamp_mux_predecode_block_drivers);
        copy_powerDef(&ptr_results->power_senseamp_mux_predecoder_blocks, tot_power_senseamp_mux_predecode_blocks);
        copy_powerDef(&ptr_results->power_senseamp_mux_decoders, tot_power_senseamp_mux_decoders);
        copy_powerDef(&ptr_results->power_bitlines, tot_power_bitlines);
        copy_powerDef(&ptr_results->power_sense_amps, tot_power_sense_amps);
        copy_powerDef(&ptr_results->power_output_drivers_at_subarray, tot_power_subarray_output_drivers);
        copy_powerDef(&ptr_results->power_dataout_vertical_htree, tot_power_dataout_vertical_htree);
        copy_powerDef(&ptr_results->power_comparators, tot_power_comparators);
        copy_powerDef(&ptr_results->total_power, tot_power);
        ptr_results->area = all_banks.height * all_banks.width * 1e-6;//in mm2
        ptr_results->all_banks_height = all_banks.height;
        ptr_results->all_banks_width = all_banks.width;		
        ptr_results->bank_height = bank.height;
        ptr_results->bank_width = bank.width;
        ptr_results->subarray_memory_cell_area_height = subarray.height;
        ptr_results->subarray_memory_cell_area_width = subarray.width;
        ptr_results->mat_height = mat.height;
        ptr_results->mat_width = mat.width;
        ptr_results->routing_area_height_within_bank = bank.height - number_mats_vertical_direction * 
            mat.height;
		ptr_results->routing_area_width_within_bank = bank.width - number_mats_horizontal_direction * 
            mat.width;
        ptr_results->area_efficiency = area_all_dataramcells * 100 / 
            (all_banks.height * all_banks.width);
        ptr_results->refresh_power = refresh_power;
        ptr_results->delay_before_subarray_output_driver = delay_before_subarray_output_driver;
        ptr_results->delay_from_subarray_output_driver_to_output = delay_from_subarray_output_driver_to_output;
    }
    else{
        ptr_array->Ndwl = Ndwl;
        ptr_array->Ndbl = Ndbl;
        ptr_array->Nspd = Nspd;
        ptr_array->deg_bitline_muxing = deg_bitline_muxing;
        ptr_array->Ndsam = Ndsam;
        ptr_array->access_time = access_time;
        ptr_array->cycle_time = cycle_time;
        ptr_array->multisubbank_interleave_cycle_time = multisubbank_interleave_cycle_time;
        ptr_array->area_ram_cells = area_all_dataramcells;
        ptr_array->area = all_banks.height * all_banks.width;
        copy_powerDef(&ptr_array->power, tot_power);
        ptr_array->delay_senseamp_mux_decoder = max_delay_before_senseamp_mux_decoder + senseamp_mux_dec.delay;
        ptr_array->delay_before_subarray_output_driver = delay_before_subarray_output_driver;
        ptr_array->delay_from_subarray_output_driver_to_output = delay_from_subarray_output_driver_to_output;
    }

    return(TRUE);
partition_not_valid:
    return(FALSE);
}

//void fa_tag_delay (input_params_t *parameters, int C,int B,int Ntwl,int Ntbl,int Ntspd, mem_array *tag_array)
void fa_tag_delay (input_params_t *parameters, int Ntwl,int Ntbl,int Ntspd, results_mem_array *tag_array)
{
    int C = parameters->cache_size/parameters->uca_banks;
    int B = parameters->block_size;

	double Tagdrive, Tag1, Tag2, Tag3, Tag4, Tag5, outrisetime;
    //int nor_inputs;
    //double Ceq, Req, Rwire, rows, tf, nextinputtime, vth, tstep;
    double Ceq, Req, Rwire, rows, tf, nextinputtime;
    //int numstack;
    double Tagdrive1 = 0, Tagdrive2 = 0;
    double Tagdrive0a = 0, Tagdrive0b = 0;
    double TagdriveA = 0, TagdriveB = 0;
    double TagdriveA1 = 0, TagdriveB1 = 0;
    double TagdriveA2 = 0, TagdriveB2 = 0;


	double FACwordmetal, FACbitmetal, FARbitmetal, FARwordmetal, dynPower;

	double WIREWIDTH, WIRESPACING, WIREPITCH;


	WIREWIDTH = 1.6*FEATURESIZE;
	WIRESPACING = 1.6*FEATURESIZE;
	WIREPITCH = WIREWIDTH + WIRESPACING;

	FACwordmetal = wire_cap(
						(32*FEATURESIZE + BitWidth_sram + 
						2*WIREPITCH * (parameters->excl_write_ports + parameters->rw_ports - 1) + 
						WIREPITCH * (parameters->excl_read_ports))*1e-6 /*m*/,
						2, 2);
	FACbitmetal = wire_cap(
						(32*FEATURESIZE + BitHeight_sram +
						2*WIREPITCH * (parameters->excl_write_ports + parameters->excl_read_ports - 1 + parameters->excl_read_ports))*1e-6 /*m*/,
						2, 2);
	FARbitmetal = wire_res(
						(32*FEATURESIZE + BitHeight_sram +
						2*WIREPITCH * (parameters->excl_write_ports + parameters->excl_read_ports - 1 + parameters->excl_read_ports))*1e-6 /* in m */,
						2, 2);		
	FARwordmetal = wire_res (
		(32*FEATURESIZE + BitWidth_sram + 
						2*WIREPITCH * (parameters->excl_write_ports + parameters->rw_ports - 1) + 
						WIREPITCH * (parameters->excl_read_ports))*1e-6/*m*/,
						2, 2);
    dynPower = 0.0;



    rows = C / (B * Ntbl);

    /* Calculate rise time.  Consider two inverters */

    Ceq = draincap (Wdecdrivep, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wdecdriven, NCH, 1, 1, DEFAULTHEIGHTCELL) +
        gatecap (Wdecdrivep + Wdecdriven, 0.0);
    tf = Ceq * transreson (Wdecdriven, NCH, 1);
    nextinputtime = horowitz (0.0, tf, VTHINV360x240, VTHINV360x240, FALL) /
        (VTHINV360x240);

    Ceq = draincap (Wdecdrivep, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wdecdriven, NCH, 1, 1, DEFAULTHEIGHTCELL) +
        gatecap (Wdecdrivep + Wdecdriven, 0.0);
    tf = Ceq * transreson (Wdecdrivep, PCH, 1);
    nextinputtime = horowitz (nextinputtime, tf, VTHINV360x240, VTHINV360x240,
            RISE) / (1.0 - VTHINV360x240);

    // If tag bitline divisions, add extra driver

    if (Ntbl > 1)
    {
        Ceq = draincap (Wfadrivep, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadriven, NCH, 1, 1, DEFAULTHEIGHTCELL) +
            gatecap (Wfadrive2p + Wfadrive2n, 0.0);
        tf = Ceq * transreson (Wfadriven, NCH, 1);
        TagdriveA = horowitz (nextinputtime, tf, VSINV, VSINV, FALL);
        nextinputtime = TagdriveA / (VSINV);
        dynPower += .5 * Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS;

        if (Ntbl <= 4)
        {
            Ceq =
                draincap (Wfadrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadrive2n, NCH, 1, 1, DEFAULTHEIGHTCELL) +
                gatecap (Wfadecdrive1p + Wfadecdrive1n,
                        10.0 / FUDGEFACTOR) * 2 + +FACwordmetal * sqrt ((rows +
                                1) * Ntbl) / 2 +
                            FACbitmetal * sqrt ((rows + 1) * Ntbl) / 2;
            Rwire =
                FARwordmetal * sqrt ((rows + 1) * Ntbl) * .5 / 2 +
                FARbitmetal * sqrt ((rows + 1) * Ntbl) * .5 / 2;
            tf = Ceq * (transreson (Wfadrive2p, PCH, 1) + Rwire);
            TagdriveB = horowitz (nextinputtime, tf, VSINV, VSINV, RISE);
            nextinputtime = TagdriveB / (1.0 - VSINV);
            dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * 2;
        }
        else
        {
            Ceq =
                draincap (Wfadrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadrive2n, NCH,
                        1, 1, DEFAULTHEIGHTCELL) +
                gatecap (Wfadrivep + Wfadriven,
                        10.0 / FUDGEFACTOR) * 2 + +FACwordmetal * sqrt ((rows +
                                1) * Ntbl) / 2 +
                            FACbitmetal * sqrt ((rows + 1) * Ntbl) / 2;
            Rwire =
                FARwordmetal * sqrt ((rows + 1) * Ntbl) * .5 / 2 +
                FARbitmetal * sqrt ((rows + 1) * Ntbl) * .5 / 2;
            tf = Ceq * (transreson (Wfadrive2p, PCH, 1) + Rwire);
            TagdriveB = horowitz (nextinputtime, tf, VSINV, VSINV, RISE);
            nextinputtime = TagdriveB / (1.0 - VSINV);
            dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * 4;

            Ceq = draincap (Wfadrivep, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadriven, NCH, 1, 1, DEFAULTHEIGHTCELL) +
                gatecap (Wfadrive2p + Wfadrive2n, 10.0 / FUDGEFACTOR);
            tf = Ceq * transreson (Wfadriven, NCH, 1);
            TagdriveA1 = horowitz (nextinputtime, tf, VSINV, VSINV, FALL);
            nextinputtime = TagdriveA1 / (VSINV);
            dynPower += .5 * Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS;

            if (Ntbl <= 16)
            {
                Ceq =
                    draincap (Wfadrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadrive2n, NCH,
                            1, 1, DEFAULTHEIGHTCELL) +
                    gatecap (Wfadecdrive1p + Wfadecdrive1n,
                            10.0 / FUDGEFACTOR) * 2 + +FACwordmetal * sqrt ((rows +
                                    1) * Ntbl) / 4 +
                                FACbitmetal * sqrt ((rows + 1) * Ntbl) / 4;
                Rwire =
                    FARwordmetal * sqrt ((rows + 1) * Ntbl) * .5 / 4 +
                    FARbitmetal * sqrt ((rows + 1) * Ntbl) * .5 / 4;
                tf = Ceq * (transreson (Wfadrive2p, PCH, 1) + Rwire);
                TagdriveB1 = horowitz (nextinputtime, tf, VSINV, VSINV, RISE);
                nextinputtime = TagdriveB1 / (1.0 - VSINV);
                dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * 8;
            }
            else
            {
                Ceq =
                    draincap (Wfadrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadrive2n, NCH,
                            1, 1, DEFAULTHEIGHTCELL) +
                    gatecap (Wfadrivep + Wfadriven,
                            10.0 / FUDGEFACTOR) * 2 + +FACwordmetal * sqrt ((rows +
                                    1) * Ntbl) / 4 +
                                FACbitmetal * sqrt ((rows + 1) * Ntbl) / 4;
                Rwire =
                    FARwordmetal * sqrt ((rows + 1) * Ntbl) * .5 / 4 +
                    FARbitmetal * sqrt ((rows + 1) * Ntbl) * .5 / 4;
                tf = Ceq * (transreson (Wfadrive2p, PCH, 1) + Rwire);
                TagdriveB1 = horowitz (nextinputtime, tf, VSINV, VSINV, RISE);
                nextinputtime = TagdriveB1 / (1.0 - VSINV);
                dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * 8;

               Ceq =
                    draincap (Wfadrivep, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadriven, NCH,
                            1, 1, DEFAULTHEIGHTCELL) +
                    gatecap (Wfadrive2p + Wfadrive2n, 10.0 / FUDGEFACTOR);
                tf = Ceq * transreson (Wfadriven, NCH, 1);
                TagdriveA2 = horowitz (nextinputtime, tf, VSINV, VSINV, FALL);
                nextinputtime = TagdriveA2 / (VSINV);
                dynPower += .5 * Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS;

                Ceq =
                    draincap (Wfadrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadrive2n, NCH,
                            1, 1, DEFAULTHEIGHTCELL) +
                    gatecap (Wfadecdrive1p + Wfadecdrive1n,
                            10.0 / FUDGEFACTOR) * 2 + +FACwordmetal * sqrt ((rows +
                                    1) * Ntbl) / 8 +
                                FACbitmetal * sqrt ((rows + 1) * Ntbl) / 8;
                Rwire =
                    FARwordmetal * sqrt ((rows + 1) * Ntbl) * .5 / 8 +
                    FARbitmetal * sqrt ((rows + 1) * Ntbl) * .5 / 8;
                tf = Ceq * (transreson (Wfadrive2p, PCH, 1) + Rwire);
                TagdriveB2 = horowitz (nextinputtime, tf, VSINV, VSINV, RISE);
                nextinputtime = TagdriveB2 / (1.0 - VSINV);
                dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * 16;
            }
        }
    }

    /* Two more inverters for enable delay */

    Ceq = draincap (Wfadecdrive1p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadecdrive1n, NCH, 1, 1, DEFAULTHEIGHTCELL)
        + gatecap (Wfadecdrive2p + Wfadecdrive2n, 0.0);
    tf = Ceq * transreson (Wfadecdrive1n, NCH, 1);
    Tagdrive0a = horowitz (nextinputtime, tf, VSINV, VSINV, FALL);
    nextinputtime = Tagdrive0a / (VSINV);
    dynPower += .5 * Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * Ntbl;

    Ceq = draincap (Wfadecdrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadecdrive2n, NCH, 1, 1, DEFAULTHEIGHTCELL) +
        +gatecap (Wfadecdrivep + Wfadecdriven, 10.0 / FUDGEFACTOR)
        + gatecap (Wfadecdrive2p + Wfadecdrive2n, 10.0 / FUDGEFACTOR);
    tf = Ceq * transreson (Wfadecdrive2p, PCH, 1);
    Tagdrive0b = horowitz (nextinputtime, tf, VSINV, VSINV, RISE);
    nextinputtime = Tagdrive0b / (VSINV);
    dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * Ntbl;
   /* First stage */

    Ceq = draincap (Wfadecdrive2p, PCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfadecdrive2n, NCH, 1, 1, DEFAULTHEIGHTCELL) +
        gatecap (Wfadecdrivep + Wfadecdriven, 10.0 / FUDGEFACTOR);
    //tf = Ceq * transresswitch (Wfadecdrive2n, NCH, 1);
	tf = Ceq * 3*transreson (Wfadecdrive2n, NCH, 1);
    Tagdrive1 = horowitz (nextinputtime, tf, VSINV, VTHFA1, FALL);
    nextinputtime = Tagdrive1 / VTHFA1;
    dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * .5 * Ntbl;

    Ceq = draincap (Wfadecdrivep, PCH, 2, 1, DEFAULTHEIGHTCELL) + draincap (Wfadecdriven, NCH, 2, 1, DEFAULTHEIGHTCELL)
        + draincap (Wfaprechn, NCH, 1, 1, DEFAULTHEIGHTCELL)
        + gatecap (Wdummyn, 0/*FIXME*/) * (rows + 1) + FACbitmetal * (rows + 1);

    Rwire = FARbitmetal * (rows + 1) * .5;
    tf = (Rwire + transreson (Wfadecdrivep, PCH, 1) +
            //transresswitch (Wfadecdrivep, PCH, 1)) * Ceq;
			3 *transreson (Wfadecdrivep, PCH, 1)) * Ceq;
    Tagdrive2 = horowitz (nextinputtime, tf, VTHFA1, VTHFA2, RISE);
    nextinputtime = Tagdrive2 / (1 - VTHFA2);

    Tagdrive =
        Tagdrive1 + Tagdrive2 + TagdriveA + TagdriveB + TagdriveA1 + TagdriveA2 +
        TagdriveB1 + TagdriveB2 + Tagdrive0a + Tagdrive0b;
    dynPower += Ceq * vdd_periph_global * vdd_periph_global * ADDRESS_BITS * Ntbl;

    /* second stage */

    Ceq = .5 * ADDRESS_BITS * 2 * draincap (Wdummyn, NCH, 2, 1, DEFAULTHEIGHTCELL)
        + draincap (Wfaprechp, PCH, 1, 1, DEFAULTHEIGHTCELL)
        + gatecap (Waddrnandn + Waddrnandp, 10.0 / FUDGEFACTOR) + FACwordmetal * ADDRESS_BITS;
    Rwire = FARwordmetal * ADDRESS_BITS * .5;
    tf =
        Ceq * (Rwire + transreson (Wdummyn, NCH, 1) +
                transreson (Wdummyn, NCH, 1));
    Tag1 = horowitz (nextinputtime, tf, VTHFA2, VTHFA3, FALL);
    nextinputtime = Tag1 / VTHFA3;
    dynPower += Ceq * vdd_periph_global * vdd_periph_global * rows * Ntbl;

    /* third stage */

    Ceq = draincap (Waddrnandn, NCH, 2, 1, DEFAULTHEIGHTCELL)
        + 2 * draincap (Waddrnandp, PCH, 1, 1, DEFAULTHEIGHTCELL)
        + gatecap (Wdummyinvn + Wdummyinvp, 10.0 / FUDGEFACTOR);
    //tf = Ceq * (transresswitch (Waddrnandp, PCH, 1));
	tf = Ceq * (3*transreson (Waddrnandp, PCH, 1));
    Tag2 = horowitz (nextinputtime, tf, VTHFA3, VTHFA4, RISE);
    nextinputtime = Tag2 / (1 - VTHFA4);
    dynPower += Ceq * vdd_periph_global * vdd_periph_global * rows * Ntbl;

    /* fourth stage */

    Ceq = (rows) * (gatecap (Wfanorn + Wfanorp, 10.0 / FUDGEFACTOR))
        + draincap (Wdummyinvn, NCH, 1, 1, DEFAULTHEIGHTCELL)
        + draincap (Wdummyinvp, PCH, 1, 1, DEFAULTHEIGHTCELL) + FACbitmetal * rows;
    Rwire = FARbitmetal * rows * .5;
    Req = Rwire + transreson (Wdummyinvn, NCH, 1);
    tf = Req * Ceq;
    Tag3 = horowitz (nextinputtime, tf, VTHFA4, VTHFA5, FALL);
    outrisetime = Tag3 / VTHFA5;
    dynPower += Ceq * vdd_periph_global * vdd_periph_global * Ntbl;

    /* fifth stage */

    Ceq = draincap (Wfanorp, PCH, 2, 1, DEFAULTHEIGHTCELL)
        + 2 * draincap (Wfanorn, NCH, 1, 1, DEFAULTHEIGHTCELL) + gatecap (Wfainvn + Wfainvp, 10.0 / FUDGEFACTOR);
    tf =
        //Ceq * (transresswitch (Wfanorp, PCH, 1) + transreson (Wfanorp, PCH, 1));
		Ceq * (3*transreson (Wfanorp, PCH, 1) + transreson (Wfanorp, PCH, 1));
    Tag4 = horowitz (nextinputtime, tf, VTHFA5, VTHFA6, RISE);
    nextinputtime = Tag4 / (1 - VTHFA6);
    dynPower += Ceq * vdd_periph_global * vdd_periph_global;

    /* final stage */

    Ceq = (gatecap (Wdecinvn + Wdecinvp, 20.0 / FUDGEFACTOR) +
            +draincap (Wfainvn, NCH, 1, 1, DEFAULTHEIGHTCELL) + draincap (Wfainvp, PCH, 1, 1, DEFAULTHEIGHTCELL));
    //Req = transresswitch (Wfainvn, NCH, 1);
	Req = 3*transreson (Wfainvn, NCH, 1);
    tf = Req * Ceq;
    Tag5 = horowitz (nextinputtime, tf, VTHFA6, VSINV, FALL);
    outrisetime = Tag5 / VSINV;
    dynPower += Ceq * vdd_periph_global * vdd_periph_global;

    //      if (Ntbl==32)
    //        fprintf(stderr," 1st - %f %f\n 2nd - %f %f\n 3rd - %f %f\n 4th - %f %f\n 5th - %f %f\nPD : %f\nNAND : %f\n INV : %f\n NOR : %f\n INV : %f\n", TagdriveA*1e9, TagdriveB*1e9, TagdriveA1*1e9, TagdriveB1*1e9, TagdriveA2*1e9, TagdriveB2*1e9, Tagdrive0a*1e9,Tagdrive0b*1e9,Tagdrive1*1e9, Tagdrive2*1e9, *Tag1*1e9, *Tag2*1e9, *Tag3*1e9, *Tag4*1e9, *Tag5*1e9);
    //power->writeOp.dynamic = dynPower;
    //power->readOp.dynamic = dynPower;
	tag_array->access_time = Tagdrive + Tag1 + Tag2 + Tag3 + Tag4 + Tag5;
	tag_array->total_power.readOp.dynamic = dynPower;
	tag_array->Ndbl = Ntbl;
	tag_array->Ndwl = Ntwl;
	tag_array->Nspd = Ntspd;
    return;
}

int
fa_data_organizational_parameters_valid (input_params_t *params,int Ndwl,int Ndbl,double Nspd)
{
    int C = params->cache_size/params->uca_banks;
    int B = params->block_size;

        if (C / (2 * B * Ndbl) <= 0) {return (FALSE);}
        if (Ndbl > MAXDATAN) {return (FALSE);}
		return (TRUE);
}

/* area, power, delay, and cycle time calculation for UCA */
void
find_acc_time(uca_org_t *res)
{
    results_mem_array * data_arr, *tag_arr;
    data_arr = (res->data_array);
    /* 
     * check whether it is just a ram without a tag array 
     */
    if (res->params->pure_sram) {
        res->access_time = data_arr->access_time;
        return;
    }
    tag_arr = (res->tag_array);

    if (res->params->associativity == 1) {
        /* 
         * Direct mapped cache: FIXME Sequential access
         * doesn't make sense since direct mapped
         * caches are always optimized for speed.
         */
        res->access_time = MAX(tag_arr->access_time,
                            data_arr->access_time);
        return;
    }
    if (res->params->fast_access) {
        /* 
         * both tag and data lookup happens
         * in parallel and the entire set is sent over
         * the dara array h-tree without waiting
         * for the way-select signal
         */
        res->access_time = MAX(tag_arr->access_time,
                            data_arr->access_time);
    }
    else if (res->params->sequential_access) {
        /* 
         * tag array is accessed first. On a hit
         * the way-select signal along with the address
         * is sent to read/write the appropriate block
         * in the data array
         */
        res->access_time = tag_arr->access_time +
                            data_arr->access_time;
    }
    else {
        /*
         * normal access: tag array access and data
         * array access happens in parallel, but the
         * data array will wait for the way-select
         * and transfer only the appropriate block
         * in the h-tree
         */
        res->access_time = MAX(tag_arr->access_time +
            data_arr->delay_senseamp_mux_decoder,
            data_arr->delay_before_subarray_output_driver) +
            data_arr->delay_from_subarray_output_driver_to_output;
    }

    return;
}

void
find_power(uca_org_t *res)
{
    double d = 0, 
           l = 0;
    if (!(res->params->pure_sram)) {
        d = res->tag_array->total_power.readOp.dynamic;
        l = res->tag_array->total_power.readOp.leakage;
    }
    res->power.readOp.dynamic = 
            res->data_array->total_power.readOp.dynamic + d;
    res->power.readOp.leakage = 
            res->data_array->total_power.readOp.leakage + l;
}

void
find_area(uca_org_t *res)
{
    if (res->params->pure_sram) {
        res->cache_ht = res->data_array->all_banks_height * 1e-3;
        res->cache_len = res->data_array->all_banks_width * 1e-3;
        res->area =  res->cache_ht * res->cache_len;
        return;
    }
    res->cache_ht = MAX(res->tag_array->all_banks_height, 
                                res->data_array->all_banks_height)*1e-3 /*mm*/;
    res->cache_len = (res->tag_array->bank_width + 
                                res->data_array->all_banks_width) * 1e-3/*mm*/;
    res->area =  res->cache_ht * res->cache_len;
}

void
find_cycle(uca_org_t *res)
{
    if (res->params->pure_sram) {
        res->cycle_time = res->data_array->cycle_time;
    }
    else {
        res->cycle_time = MAX(res->tag_array->cycle_time, 
                res->data_array->cycle_time);
    }
    return;
}

void
update_min_values (uca_org_t *res, double *delay, double *dyn,
                   double *leak, double *cycle, double *area)
{
    if (res->access_time < *delay) {
        *delay = res->access_time;
    }
    if (res->cycle_time < *cycle) {
        *cycle = res->cycle_time;
    }
    if (res->power.readOp.dynamic < *dyn) {
        *dyn = res->power.readOp.dynamic;
    }
    if (res->power.readOp.leakage < *leak) {
        *leak = res->power.readOp.leakage;
    }
    if(res->area < *area) {
        *area = res->area;
    }
}

void
update_min_values2 (nuca_org_t *res, double *delay, double *dyn,
                   double *leak, double *cycle, double *area)
{
    if (res->nuca_pda.delay < *delay) {
        *delay = res->nuca_pda.delay;
    }
    if (res->router.cycle_time < *cycle) {
        *cycle = res->router.cycle_time;
    }
    if (res->nuca_pda.power.dynamic < *dyn) {
        *dyn = res->nuca_pda.power.dynamic;
    }
    if (res->nuca_pda.power.leakage < *leak) {
        *leak = res->nuca_pda.power.leakage;
    }
    if(res->nuca_pda.area_stats.area < *area) {
        *area = res->nuca_pda.area_stats.area;
    }
}

void
update_min (results_mem_array *res, double *delay, double *dyn)
{
    if (res->access_time < *delay) {
        *delay = res->access_time;
    }
    if (res->total_power.readOp.dynamic < *dyn) {
        *dyn = res->total_power.readOp.dynamic;
    }
}

/* 
 * Version - 6.0
 * Performs exhaustive search across different sub-array sizes, 
 * wire types and aspect ratios to find an optimal UCA organization
 * 1. First different valid tag array organizations are calculated
 *    and stored in tag_arr array
 * 2. The exhaustive search is repeated to find valid data array
 *    organizations and stored in data_arr array
 * 3. Cache area, delay, power, and cycle time for different
 *    cache organizations are calculated based on the 
 *    above results and stored in a linked list
 * 4. Cache model with least cost is picked from the ll
 */
uca_res_lentry_t * 
sim_uca(uca_org_t *res)
{
    /* temp variables */
    int t, d, i, j, wt, ram = 0;
    int valid_partition;
    int Ntspd, Ndwl, Ndbl, Ndcm, Ndsam;
    double Nspd; /* NOTE: Nspd can assume fractional values
                          for the data array */
    int wt_min, wt_max;

    input_params_t *parameters = res->params;
    uca_org_t *uca_res;
    double min_area = BIGNUM;
    double min_cycle = BIGNUM;
    double min_delay = BIGNUM;
    double min_dyn = BIGNUM;
    double min_leak = BIGNUM;

    /* Exhaustive search range */
	int maxdnspd = MAXDATASPD;
	int maxndwl = MAXDATAN;
	int maxndbl = MAXDATAN;
	int maxndcm = MAX_COL_MUX;

    /* to store the results of the search space */
    uca_res_lentry_t *home_lentry; /* head of the ll */
    uca_res_lentry_t *lentry, *temp; /* linked list entry */
    long int ll_size = 0; /* size of the ll */


    /* 
     * to store the results of 
     * different tag and data array organizations
     */
    results_mem_array **tag_arr, **data_arr;

    lentry = (uca_res_lentry_t *)malloc(sizeof(uca_res_lentry_t));
    memset(lentry, 0, sizeof(uca_res_lentry_t));
    lentry->next_lentry = NULL;
    home_lentry = lentry;
    lentry->fin_res.params = res->params;
    ll_size++;

	tag_arr = (results_mem_array **) calloc(MAX_NUMBER_ARRAY_PARTITIONS, 
                                sizeof(results_mem_array *));
    data_arr = (results_mem_array **) calloc(MAX_NUMBER_ARRAY_PARTITIONS, 
                                sizeof(results_mem_array *));
    tag_arr[0] = (results_mem_array *) malloc(sizeof(results_mem_array));
    data_arr[0] = (results_mem_array *) malloc(sizeof(results_mem_array));
    t = d = 0;

    is_dram = 0;
    if (parameters->fully_assoc) {
        maxdnspd = 1;
        maxndwl = 1;
        maxndcm = 1;
    }
    if (parameters->force_wiretype) {
        if (parameters->wire_inter_mats == 0) {
            wt_min = Low_swing;
            wt_max = Low_swing;
        }
        else {
            wt_min = Global;
            wt_max = Global_30;
        }
    }
    else {
        wt_min = Global;
        wt_max = Low_swing;
    }
    /*
     * Search the design space for different valid
     * tag array organizations
     */
    if(!parameters->pure_sram){/* tag array */
        for (wt=wt_min; wt<=wt_max; wt++) {
            parameters->wire_inter_mats = wt;
            for(Ntspd = 1; Ntspd <= maxdnspd; Ntspd = Ntspd * 2){
                for(Ndwl = 1; Ndwl <= maxndwl; Ndwl = Ndwl * 2){
                    for(Ndbl = 2; Ndbl <= maxndbl; Ndbl = Ndbl * 2){
                        for(Ndcm = 1; Ndcm <= maxndcm; Ndcm = Ndcm * 2){
                            for(Ndsam = 1; Ndsam <= maxndcm; Ndsam = Ndsam * 2){
                                is_dram = 0; 

                                /* 
                                 * The tag array of the fully associative model
                                 * is a CAM structure. Since FA caches 
                                 * are typically small
                                 * the exhaustive search is limited to different
                                 * Ndbl values. ref:CACTI-2 Tech report
                                 * NOTE: If ndwl or nspd values are modified
                                 * then fa_tag_delay function and the circuit 
                                 * model should also be changed
                                 */
                                if (parameters->fully_assoc) {
                                    valid_partition = 
                                        fa_data_organizational_parameters_valid (
                                                parameters, 1 /* Ndwl */,Ndbl,
                                                1 /* Ntspd */);

                                    if (valid_partition) {
                                        fa_tag_delay(parameters, Ndwl, Ndbl,
                                                Ntspd, tag_arr[t]);
                                    }
                                    else continue;
                                }
                                else{
                                    valid_partition = is_partition_valid(
                                            parameters, 1 /* model tag array */, 
                                            Ndwl, Ndbl,Ntspd, Ndcm, Ndsam); 
                                    if (valid_partition) { 
                                        valid_partition = calculate_time(
                                                parameters, parameters->uca_banks, 
                                                1 /* model tag array */, 
                                                parameters->pure_sram, 
                                                Ntspd, Ndwl, Ndbl, Ndcm, Ndsam, 
                                                NULL, 1, tag_arr[t], NULL);
                                    }
                                }
                                if (!valid_partition) continue;
                                tag_arr[t]->wt = parameters->wire_inter_mats;
                                
                                update_min(tag_arr[t], &min_delay, &min_dyn);
                                ++t;
                                if (t>=100000) {ERROR_EXIT(stderr, 
                                        "Tag array overflow!!\n");}
                                tag_arr[t] = (results_mem_array *) 
                                    malloc(sizeof(results_mem_array));
                                memset(tag_arr[t], 0, sizeof(results_mem_array));
                            } /* search loop */
                        }
                    }
                }
            }
        }
    }
    free(tag_arr[t]);
    PRINTD(fprintf(stderr, "sim_uca: No. of valid tag array org. = %d\n", t));

    if ((t == 0)) {
        if (!parameters->pure_sram) {
            fprintf(stderr, "\n\nERROR: No valid organization for the following input"
                    " parameters!\n");
            dump_input_args(parameters);
            exit(-1);
        }
        else t = 1;
    }
    
    /* remove sub-optimal tag-array organizations */
    filter_arr(tag_arr, min_delay, min_dyn, t, 0/*tag*/);
    min_delay = BIGNUM;
    min_dyn = BIGNUM;

    if (parameters->fully_assoc) {
        maxdnspd = MAXDATAN;
        maxndwl = MAXDATAN;
        maxndcm = MAXDATAN;
    }
    /* 
     * Search the design space for different 
     * data array organizations
     */
    for (wt=wt_min; wt<=wt_max; wt++) {
        parameters->wire_inter_mats = wt;
        for(Nspd = (double)(parameters->output_width)/
                (double)(parameters->block_size*8); 
                Nspd <= MAXDATASPD; Nspd = Nspd * 2){
            for(Ndwl = 2; Ndwl <= maxndwl; Ndwl = Ndwl * 2){
                for(Ndbl = 2; Ndbl <= maxndbl; Ndbl = Ndbl * 2){
                    for(Ndcm = 1; Ndcm <= maxndcm; Ndcm = Ndcm * 2){
                        for(Ndsam = 1; Ndsam <= maxndcm; Ndsam = Ndsam * 2){
                            is_dram = parameters->dram; 

                            valid_partition = is_partition_valid (
                                    parameters, 0,
                                    Ndwl, Ndbl, Nspd, Ndcm, Ndsam);
                            if (valid_partition) {
                                valid_partition = calculate_time(
                                        parameters, parameters->uca_banks, 
                                        0 /* model data array */, parameters->pure_sram,
                                        Nspd, Ndwl, Ndbl, Ndcm, Ndsam, NULL, 1,
                                        data_arr[d], NULL);
                            }
                            if (!valid_partition) continue;
                            data_arr[d]->wt = parameters->wire_inter_mats;
                            
                            update_min(data_arr[d], &min_delay, &min_dyn);

                            ++d;
                            if (d>=100000) {
                                ERROR_EXIT(stderr, "Data array overflow!!\n");}
                            data_arr[d] = (results_mem_array *) 
                                malloc(sizeof(results_mem_array));
                            memset(data_arr[d], 0, sizeof(results_mem_array));
                        } /* search loop */
                    }
                }
            }
        }
    }
    free(data_arr[d]);
    PRINTD(fprintf(stderr, "sim_uca: No. of valid data array org. = %d\n", d));
    if (d == 0) {
        fprintf(stderr, "\n\nERROR: No valid organization for the following input"
            " parameters!\n");
        dump_input_args(parameters);
        exit(-1);
    }

    /* remove sub-optimal data array organizations */
    filter_arr(data_arr, min_delay, min_dyn, d, 1/*data*/);
    min_delay = BIGNUM;
    min_dyn = BIGNUM;

    /* 
     * Calculate the delay, power, area, and cycle time for different
     * cache organizations
     */
    for (i=0; i<t; i++) {
        if(!tag_arr[i]) {
            if (!(parameters->pure_sram)) continue;
            else {
                ram = 1;
            }
        }
        for (j=0; j<d; j++) {
            if(!data_arr[j]) continue;

            if(!ram) lentry->fin_res.tag_array = tag_arr[i];
            lentry->fin_res.data_array = data_arr[j];
            if (data_arr[j]->access_time == 0) {
                fprintf(stderr, "sim_uca: i %d, j%d\n", i, j);
                exit(-1);
            }
                        
            find_acc_time(&(lentry->fin_res));
            find_power(&(lentry->fin_res));
            find_area(&(lentry->fin_res));
            find_cycle(&(lentry->fin_res));

            /* 
             * check whether the current org. gives min. 
             * delay/power/area/cycle time
             */
            update_min_values (&(lentry->fin_res), &min_delay, &min_dyn,
                    &min_leak, &min_cycle, &min_area);

            lentry->next_lentry = (uca_res_lentry_t *) 
                            malloc(sizeof(uca_res_lentry_t));
            temp = lentry;
            lentry = lentry->next_lentry;
            memset(lentry, 0, sizeof(uca_res_lentry_t));
            lentry->fin_res.params = res->params;
            lentry->next_lentry = NULL;
            ll_size++;
        }
    }
    
    if (!temp) {
        ERROR_EXIT(stdout, "ERROR: sim_uca: lentry is NULL!\n");
    }
    temp->next_lentry = NULL;
    free(lentry);

    PRINTD(fprintf(stderr, "sim_uca: No. of valid cache org. = %d\n", ll_size));

    uca_res = find_optimal_uca(home_lentry, min_delay, min_dyn, 
            min_leak, min_cycle, min_area, parameters);

    if(uca_res == NULL) {
        printf("ERROR: There are no valid organization that meets the"
                " current input area/delay/power/cycle time constraints!"
                " Try changing the input deviation"
                " values.\n");
        exit(-1);
    }
    
    update_output(res, uca_res);

    if(parameters->pure_sram) t=0;
    free_mem(tag_arr, t, data_arr, d, home_lentry);

    PRINTD(fprintf(stderr, "sim_uca: Min delay %g, Min dynamic power %g," 
        "Min leakage power %g, Min cycle time %g, Min area %g\n",
        min_delay, min_dyn, min_leak, min_cycle, min_area));
    //fig_out(fin_res);
    return NULL;
}

void
update_output (uca_org_t *res, uca_org_t *uca_res) 
{
    results_mem_array *t, *d;
    t = res->tag_array;
    d = res->data_array;
    if (t) bcopy(uca_res->tag_array, res->tag_array, sizeof(results_mem_array));
    bcopy(uca_res->data_array, res->data_array, sizeof(results_mem_array));
    bcopy(uca_res, res, sizeof(uca_org_t));
    res->tag_array = t;
    res->data_array = d;
}

/* converts latency (in s) to cycles depending upon the FREQUENCY (in GHz) */
int  
calc_cycles(double lat, double oper_freq) 
{
    //TODO: convert latch delay to FO4 */
    double cycle_time = (1.0/(oper_freq*1e9)); /*s*/
    cycle_time -= LATCH_DELAY;
    cycle_time -= FIXED_OVERHEAD;

    PRINTDW(fprintf(stderr, "calc_cycles: Lat  = %e\n", lat));
    PRINTDW(fprintf(stderr, "calc_cycles: Cycle  = %f\n", ceil(lat/cycle_time)));
    return (int)ceil(lat/cycle_time);
}

void
extract_bank_stats (pda_res_t *bank_pda, uca_org_t *bank_res, 
            pda_res_t *router_pda)
{    
    bank_pda->power.dynamic = bank_res->power.readOp.dynamic;
    bank_pda->power.leakage = bank_res->power.readOp.leakage;
    bank_pda->delay = bank_res->access_time;
//    bank_pda->cycle_time = calc_cycles(bank_res->cycle_time); 
    bank_pda->area_stats.area = bank_res->area;
    /* 
     * FIXME: the net area of a NUCA bank = area of the rectangular
     * region bounding the tag & data arrays, and the router 
     * NOTE: In CACTI 5 area is sum of tag array area and data
     * array area.
     */
    bank_pda->area_stats.height = MAX(bank_res->tag_array->all_banks_height +
            router_pda->area_stats.height, 
            bank_res->data_array->all_banks_height) * 1e-3 /* m */;
    bank_pda->area_stats.width = MAX( bank_res->tag_array->all_banks_width + 
            bank_res->data_array->all_banks_width,
            bank_res->data_array->all_banks_width + 
            router_pda->area_stats.width) * 1e-3; // m
}

void
calculate_nuca_area (nuca_org_t *nuca) 
{
    nuca->nuca_pda.area_stats.height = 
            nuca->rows * ((nuca->h_wire.wire_width +
                          nuca->h_wire.wire_spacing)
            * nuca->router.flit_size +
            nuca->bank_pda.area_stats.height);

    nuca->nuca_pda.area_stats.width = 
            nuca->columns * ((nuca->v_wire.wire_width +
                          nuca->v_wire.wire_spacing)
            * nuca->router.flit_size +
            nuca->bank_pda.area_stats.width);

    nuca->nuca_pda.area_stats.area = nuca->nuca_pda.area_stats.width *
            nuca->nuca_pda.area_stats.height;
}

/* 
 * Version - 6.0 
 *
 * Perform exhaustive search across different bank organizatons,
 * router configurations, grid organizations, and wire models and
 * find an optimal NUCA organization
 * For different bank count values
 * 1. Optimal bank organization is calculated  
 * 2. For each bank organization, find different NUCA organizations
 *    using various router configurations, grid organizations,
 *    and wire models. 
 * 3. NUCA model with the least cost is picked for 
 *    this particular bank count 
 * Finally include contention statistics and find the optimal
 *    NUCA configuration
 */
void 
sim_nuca(nuca_org_t *res)
{
    /* temp variables */
    int it, ro, wr;
    int num_cyc;
    double dn_wt;
    int i, j, k;
    int r, c;
    int l2_c;
    int bank_count = 0;
    uca_org_t ures;
    nuca_org_t *opt_n;
    results_mem_array tag, data;
	long int ll_size;
    int wt_min, wt_max;

    /* to search diff grid organizations */
    float curr_hop, totno_hops, totno_hhops, totno_vhops, tot_lat,
        curr_acclat;
    double avg_lat, avg_hop, avg_hhop, avg_vhop, avg_dyn_power,
        avg_leakage_power;

    double opt_acclat = INF, opt_avg_lat = INF, opt_tot_lat = INF;
    int opt_rows = 0;
    int opt_columns = 0;
    double opt_totno_hops = 0;
    double opt_avg_hop = 0;
    double opt_dyn_power = 0, opt_leakage_power = 0;

    double min_delay = BIGNUM, 
           min_dyn = BIGNUM,
           min_area = BIGNUM,
           min_cycle = BIGNUM,
           min_leak = BIGNUM;

    int stride;
    int bank_start = 0;

    /* cummulative_access_stats.dat 
     * file has details at a granularity of 32KB 
     */


    wire_stats_t wire_vertical[WIRE_TYPES],
                 wire_horizontal[WIRE_TYPES];
    int flit_width = 0;

    /* vertical and horizontal hop latency values */
    int ver_hop_lat, hor_hop_lat; /* in cycles */

    uca_res_lentry_t *lentry, *temp; /* linked list entry */

    //results
    nuca_res_lentry_t *home_nuca_lentry;
    nuca_res_lentry_t *nuca_lentry, *nuca_temp;

    /* no. of different bank sizes to consider */
    int iterations;
    
    input_params_t *params = (input_params_t *) malloc(sizeof(input_params_t));
    bcopy(res->params, params, sizeof(input_params_t)); 
    
	lentry = (uca_res_lentry_t *) malloc(sizeof(uca_res_lentry_t)); 
    temp = lentry;
    ll_size = 0; 

    nuca_lentry = (nuca_res_lentry_t *)malloc(sizeof(nuca_res_lentry_t));
    memset(nuca_lentry, 0, sizeof(nuca_res_lentry_t));
    home_nuca_lentry = nuca_lentry;
    nuca_lentry->next_lentry = NULL;

    if (params->cache_level == 0) l2_c = 1;
    else l2_c = 0;
    
    if (params->cores <= 4) core_in = 2;
    else if (params->cores <= 8) core_in = 3;
    else if (params->cores <= 16) core_in = 4;
    else {ERROR_EXIT (stderr, "Number of cores should be <= 16!\n");}


    if (params->associativity > 2) {
        i = 2;
        while (i != params->associativity) {
            MIN_BANKSIZE *= 2;
            i *= 2;
        }
    }

    stride = params->cache_size/MIN_BANKSIZE; 

    iterations = (int)logbasetwo((int)params->cache_size/MIN_BANKSIZE);

    if (params->force_wiretype) {
        if (params->wire_inter_mats == 0) {
            wt_min = Low_swing;
            wt_max = Low_swing;
        }
        else {
            wt_min = Global;
            wt_max = Low_swing-1;
        }
    }
    else {
        wt_min = Global;
        wt_max = Low_swing;
    }
    if (params->nuca_bank_count != 0) {
        if (params->nuca_bank_count != 2 && params->nuca_bank_count != 4 &&
            params->nuca_bank_count != 8 && params->nuca_bank_count != 16 &&
            params->nuca_bank_count != 32 && params->nuca_bank_count != 64) {
            fprintf(stderr,"Incorrect bank count value! Please fix the value in cache.cfg\n");
        }
        bank_start = (int)logbasetwo((double)params->nuca_bank_count);
        iterations = bank_start+1;
        params->cache_size = res->params->cache_size/params->nuca_bank_count;
    }
    for (it=bank_start; it<iterations; it++) { /* different bank count values */
        ures.params = params;
        ures.tag_array = &tag;
        ures.data_array = &data;
        /* 
         * find the optimal bank organization 
         */
        sim_uca(&ures); 
        bank_count = res->params->cache_size/params->cache_size;
        /* HACK FIXME*/
        lentry = temp;
        lentry->next_lentry = NULL;
        bcopy(&ures, &(lentry->fin_res), sizeof(uca_org_t));
        PRINTD(printf("Simulating NUCA with bank size %ld\n", ures.params->cache_size));
       
        while (lentry) {
            for (wr=wt_min; wr<=wt_max; wr++) {
//            for (wr=wt_min; wr<=wt_min; wr++) {
                for (ro=0; ro<ROUTER_TYPES; ro++) {
//                for (ro=0; ro<1; ro++) {
                    flit_width = router_s[ro].flit_size;
                    nuca_lentry->nres.nuca_pda.cycle_time = MAX(router_s[ro].cycle_time,
                                router_s[ro].max_cyc);

                    /* calculate router and wire parameters */

                    wire_vertical[wr].wire_length = 
                        lentry->fin_res.cache_ht*1e-3; /* length of the wire (m)*/ 
                    wire_horizontal[wr].wire_length = 
                        lentry->fin_res.cache_len*1e-3; /* length of the wire (m)*/ 
                    wire_vertical[wr].nsense = wire_horizontal[wr].nsense = 1;
                    wire_vertical[wr].wt = wire_horizontal[wr].wt = wr;

                    /* find delay, area, and power for wires */
                    calc_wire_stats2((enum wire_type) wr,
                            &(wire_vertical[wr]));
                    calc_wire_stats2((enum wire_type) wr,
                            &(wire_horizontal[wr]));

                    hor_hop_lat = calc_cycles(wire_horizontal[wr].wire_pda.delay,
                        1/(nuca_lentry->nres.nuca_pda.cycle_time*.001));
                    ver_hop_lat = calc_cycles(wire_vertical[wr].wire_pda.delay,
                        1/(nuca_lentry->nres.nuca_pda.cycle_time*.001));
                    /*
                     * assume a grid like topology and explore for optimal network
                     * configuration using different row and column count values.
                     */
                    for (c=1; c<=bank_count; c++) {
                        while (bank_count%c != 0) c++;
                        r = bank_count/c;

                        /* 
                         * to find the avg access latency of a NUCA cache, uncontended
                         * access time to each bank from the 
                         * cache controller is calculated. 
                         * avg latency = 
                         * sum of the access latencies to individual banks)/bank 
                         * count value.
                         */
                        totno_hops = totno_hhops = totno_vhops = tot_lat = 0;
                        k = 1;
//                        fprintf(stderr, "------------------------\n");
                        for (i=0; i<r; i++) {
                            for (j=0; j<c; j++) {
                                /* 
                                 * vertical hops including the
                                 * first hop from the cache controller 
                                 */
                                curr_hop = i + 1;
//                                curr_hop = i + 0; 
                                curr_hop += j; /* horizontal hops */
                                totno_hhops += j;
                                totno_vhops += (i+1);
                                curr_acclat = (i * ver_hop_lat + CONTR_2_BANK_LAT +
                                        j * hor_hop_lat);

#ifdef DNUCA
#endif
                                tot_lat += curr_acclat;
                                totno_hops += curr_hop;
                            }
                        }
                        avg_lat = tot_lat/bank_count;
#ifndef DNUCA
                        avg_hop = totno_hops/bank_count;
#else
                        avg_hop = totno_hops;
#endif
                        avg_hhop = totno_hhops/bank_count;
                        avg_vhop = totno_vhops/bank_count;

                        /* net access latency */
                        curr_acclat = 2*avg_lat + 2*(ROUTER_LAT*avg_hop) + 
                            calc_cycles(lentry->fin_res.access_time,
                                1/(nuca_lentry->nres.nuca_pda.cycle_time*.001));
//                        fprintf(stderr, "acc_lat = %g\n", curr_acclat);
//#ifndef DNUCA

                        /* avg access lat of nuca */
                        avg_dyn_power =
                            avg_hop * 
                            (router_s[ro].router_pda.power.dynamic) + avg_hhop *
                            (wire_horizontal[wr].wire_pda.power.dynamic) *
                            (params->block_size*8 + 64) + avg_vhop * 
                            (wire_vertical[wr].wire_pda.power.dynamic) * 
                            (params->block_size*8 + 64) + lentry->fin_res.power.readOp.dynamic;
//#endif
                        avg_leakage_power = 
                         bank_count * router_s[ro].router_pda.power.leakage +
                         (c-1)*r * (wire_horizontal[wr].wire_pda.power.leakage*
                         wire_horizontal[wr].wire_pda.delay) * flit_width +
                         (r-1)*c * (wire_vertical[wr].wire_pda.power.leakage *
                         wire_horizontal[wr].wire_pda.delay) * flit_width;

                        if (curr_acclat < opt_acclat) {
                            opt_acclat = curr_acclat;
                            opt_tot_lat = tot_lat;
                            opt_avg_lat = avg_lat;
                            opt_totno_hops = totno_hops;
                            opt_avg_hop = avg_hop;
                            opt_rows = r;
                            opt_columns = c;
                            opt_dyn_power = avg_dyn_power;
                            opt_leakage_power = avg_leakage_power;
                        }
                        totno_hops = 0;
                        tot_lat = 0;
                        totno_hhops = 0;
                        totno_vhops = 0;
                    }
                    nuca_lentry->nres.wire_pda.power.dynamic =
                        opt_avg_hop * flit_width *
                        (wire_horizontal[wr].wire_pda.power.dynamic +
                        wire_vertical[wr].wire_pda.power.dynamic);
                    nuca_lentry->nres.avg_hops = opt_avg_hop;
                    nuca_lentry->nres.params = ures.params;
                    /* network delay/power */
                    nuca_lentry->nres.h_wire = wire_horizontal[wr];
                    nuca_lentry->nres.v_wire = wire_vertical[wr];
                    nuca_lentry->nres.router = router_s[ro];
                    /* bank delay/power */
                    extract_bank_stats(&(nuca_lentry->nres.bank_pda), 
                    &(lentry->fin_res), &(nuca_lentry->nres.router.router_pda));

                    /* net power, delay, & cycle time */
//                    nuca_lentry->nres.nuca_pda.cycle_time = MAX(router_s[ro].cycle_time,
//                                router_s[ro].max_cyc);
                    num_cyc = calc_cycles(nuca_lentry->nres.bank_pda.delay /*s*/,
                        1/(nuca_lentry->nres.nuca_pda.cycle_time*.001/*GHz*/));
                    if(num_cyc%2 != 0) num_cyc++;
                    if (num_cyc > 16) num_cyc = 16;
#ifndef DNUCA
                    if (it < 7) {
                        nuca_lentry->nres.nuca_pda.delay = opt_acclat + 
                            cont_stats[l2_c][core_in][ro][it][num_cyc/2-1];
                        nuca_lentry->nres.contention = 
                            cont_stats[l2_c][core_in][ro][it][num_cyc/2-1];
                    }
                    else {
                        nuca_lentry->nres.nuca_pda.delay = opt_acclat + 
                            cont_stats[l2_c][core_in][ro][7][num_cyc/2-1];
                        nuca_lentry->nres.contention = 
                            cont_stats[l2_c][core_in][ro][7][num_cyc/2-1];
                    }
#else
                    nuca_lentry->nres.nuca_pda.delay = opt_acclat;
#endif
                    nuca_lentry->nres.nuca_pda.power.dynamic = opt_dyn_power;
                    nuca_lentry->nres.nuca_pda.power.leakage = opt_leakage_power;

                    /* array organization */
                    nuca_lentry->nres.size = ures.params->cache_size;
                    nuca_lentry->nres.bank_count = bank_count;
                    nuca_lentry->nres.rows = opt_rows;
                    nuca_lentry->nres.columns = opt_columns;
                    calculate_nuca_area (&(nuca_lentry->nres));

                    update_min_values2 (&(nuca_lentry->nres), &min_delay, &min_dyn,
                                &min_leak, &min_cycle, &min_area);
                    PRINTD(print_nuca_pda(&(nuca_lentry->nres)));

                    nuca_lentry->next_lentry = (nuca_res_lentry_t *) 
                                            malloc(sizeof(nuca_res_lentry_t));
                    nuca_temp = nuca_lentry;
                    nuca_lentry = nuca_lentry->next_lentry;
                    
                    opt_acclat = BIGNUM;
                }
            }
            temp = lentry; // HACK FIXME
            lentry = lentry->next_lentry;
        }
        params->cache_size /= 2;
    }

    free(lentry); //HACK FIXME
    nuca_temp->next_lentry = NULL;
    free(nuca_lentry);
    ures.params->cache_size /= 2;
    opt_n = find_optimal_nuca(home_nuca_lentry, min_delay, min_dyn, 
            min_leak, min_cycle, min_area);
    assert(opt_n);
    opt_n->params = res->params;
    bcopy(opt_n, res, sizeof(nuca_org_t));
    free_mem2 (home_nuca_lentry);

    free(params);
}





void
print_nuca_pda (nuca_org_t *nres) 
{
    fprintf(stderr, "NUCA stats:\n");
    fprintf(stderr, "\tDelay: %g\n", nres->nuca_pda.delay);
    fprintf(stderr, "\tDynamic power: %g\n", nres->nuca_pda.power.dynamic);
    fprintf(stderr, "\tLeakage power: %g\n", nres->nuca_pda.power.leakage);
    fprintf(stderr, "\tCycle time: %g\n\n", nres->nuca_pda.cycle_time);
    fprintf(stderr, "\tCache height x width (in mm): %g x %g\n\n", 
                    nres->nuca_pda.area_stats.height * 1e3, 
                    nres->nuca_pda.area_stats.width * 1e3);

    fprintf(stderr, "\tBank size %d\n", nres->size);
    fprintf(stderr, "\tNumber of banks %d\n", nres->bank_count);
    fprintf(stderr, "\tRows %d\n", nres->rows);
    fprintf(stderr, "\tColumns %d\n", nres->columns);
    fprintf(stderr, "\tWire type %d\n", nres->h_wire.wt);
    fprintf(stderr, "\tVertical hop latency %d\n", 
                    calc_cycles(nres->v_wire.wire_pda.delay,
                        1/(nres->nuca_pda.cycle_time*.001)));
    fprintf(stderr, "\tHorizontal hop latency %d\n\n",
                    calc_cycles(nres->h_wire.wire_pda.delay,
                        1/(nres->nuca_pda.cycle_time*.001)));
    fprintf(stderr, "\tBank energy - %g(nJ)\n", nres->bank_pda.power.dynamic*1e9);
    fprintf(stderr, "\tWire energy - %g(nJ)\n", nres->wire_pda.power.dynamic*1e9);
    fprintf(stderr, "\tRouter energy - %g(nJ)\n", nres->router.router_pda.power.dynamic*1e9*
        nres->avg_hops);
}

void
free_mem (results_mem_array **tag_arr, int t, results_mem_array **data_arr,
        int d, uca_res_lentry_t *ll) 
{
    int i;
    int te = 0;
    uca_res_lentry_t *temp;

    for (i=0; i<t; i++) {
        if(tag_arr[i]) free(tag_arr[i]);
    }
    free (tag_arr);
    for (i=0; i<d; i++) {
        if(data_arr[i]) free(data_arr[i]);
    }
    free (data_arr);
    while(ll) {
        temp = ll;
        ll = ll->next_lentry;
        free(temp);
        te++;
    }
}

void
free_mem2 (nuca_res_lentry_t *ll) 
{
    nuca_res_lentry_t *temp;

    while(ll) {
        temp = ll;
        ll = ll->next_lentry;
        free(temp);
    }
}


void
print_wire(wire_stats_t *wire)
{
    fprintf(stderr, "\n\nWire stats:\n");
    if (wire->wt == Global) {
        fprintf(stderr, "\tWire type - Full swing global wires\n");
    }
    else if (wire->wt == Global_10) {
        fprintf(stderr, "\tWire type - Full swing global wires with "
                "10%% delay penalty\n");
    }
    else if (wire->wt == Global_20) {
        fprintf(stderr, "\tWire type - Full swing global wires with "
                "20%% delay penalty\n");
    }
    else if (wire->wt == Global_30) {
        fprintf(stderr, "\tWire type - Full swing global wires with "
                "30%% delay penalty\n");
    }
    else if(wire->wt == Low_swing) {
        fprintf(stderr, "\tWire type - Low swing wires\n");
    }

    fprintf(stderr, "\tWire length - %g (mm)\n",
            wire->wire_length*1e3);

    fprintf(stderr, "\tWire width - %g\n", 
            wire->wire_width);
    fprintf(stderr, "\tWire spacing - %g\n", 
            wire->wire_spacing);
    fprintf(stderr, "\tWire delay - %g\n", 
            wire->wire_pda.delay);
    fprintf(stderr, "\tWire energy -dynamic %g (J)\n"
            "\t            -leakage %g (W)\n", 
            wire->wire_pda.power.dynamic,
            wire->wire_pda.power.leakage);
}

void
dump_input_args(input_params_t *b)
{
    if (b) {
        fprintf(stderr, "dump_input_args: \nMost common parameters\n\n"
                "\t Cache size - %ld\n"
                "\t Block size - %d\n"
                "\t Associativity - %d\n"
                "\t RW Ports - %d\n"
                "\t ER Ports - %d\n"
                "\t EW Ports - %d\n"
                "\t UCA banks - %d\n"
                "\t Technology - %f\n"
                "\t Output Width - %d\n"
                "\t Access mode - %d\n"
                "\t Scratch RAM - %d\n"
                "\t Temperature - %d\n"
                "\t DRAM - %d\n\nObjective Function Params\n\n"
                "\t Weight Delay:Dyn Power:Leak Power:Cyc Time:Area  "
                "::  %d:%d:%d:%d:%d\n"
                "\t Deviate Delay:Dyn Power:Leak Power:Cyc Time:Area "
                "::  %d:%d:%d:%d:%d\n\n Other parameters\n\n"
                "\t SERP - %d\n"
                "\t Cores - %d\n"
                "\t Level - %d\n\n",
                b->cache_size, b->block_size, b->associativity,
                b->rw_ports, b->excl_read_ports, b->excl_write_ports,
                b->uca_banks,
                b->tech_size, b->output_width, 
                b->access_mode, b->pure_sram,
                b->temp, b->dram, b->delay_wt, b->dynamic_power_wt,
                b->leakage_power_wt, b->cycle_time_wt, b->area_wt,
                b->delay_dev, b->dynamic_power_dev, b->leakage_power_dev,
                b->cycle_time_dev, b->area_dev,
                b->single_ended_read_ports, b->cores, b->cache_level+1); 
        if(b->force_nuca_bank)
            fprintf(stderr, "\t NUCA banks - %d\n", b->nuca_bank_count);
        
        if(b->force_wiretype)
            fprintf(stderr, "\t Wire Type - %d\n", (int) b->wire_inter_mats);

        if(b->force_tag)
            fprintf(stderr, "\t Tag width - %d\n", b->tag_size);

        fprintf(stderr, "NUCA wt and dev %d:%d:%d:%d:%d\n%d:%d:%d:%d:%d\n",
                   b->delay_wt_nuca, b->dynamic_power_wt_nuca,
                   b->leakage_power_wt_nuca, b->cycle_time_wt_nuca, b->area_wt_nuca,
                   b->delay_dev_nuca, b->dynamic_power_dev_nuca, b->leakage_power_dev_nuca,
                   b->cycle_time_dev_nuca, b->area_dev_nuca);
    }
}

