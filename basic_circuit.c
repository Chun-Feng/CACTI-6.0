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


#include "basic_circuit.h"

int powers (int base, int n)
{
  int i, p;

  p = 1;
  for (i = 1; i <= n; ++i)
    p *= base;
  return p;
}

/*----------------------------------------------------------------------*/

double logtwo (double x)
{
  if (x <= 0)
    printf ("logtwo: Error %e\n", x);
  return ((double) (log (x) / log (2.0)));
}

/*----------------------------------------------------------------------*/

double gatecap (double width,double  wirelength)	/* returns gate capacitance in Farads */
     /* width: gate width in um (length is Lphy_periph_global) */
     /* wirelength: poly wire length going to gate in lambda */
{
	double cap;
	if((is_dram)&&(is_access_transistor)){//DRAM cell access transistor
		cap = (c_g_ideal_itrs_dram_access_transistor + 2 * c_fringe_itrs_dram_access_transistor) * width +
			Lphy_dram_access_transistor * Cpolywire;
	}
	else if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
		cap = (c_g_ideal_itrs_dram_wordline_transistor + 2 * c_fringe_itrs_dram_wordline_transistor) * width +
			Lphy_dram_wordline_transistor * Cpolywire;
	}
	else if((!is_dram)&&(is_access_transistor)){//SRAM cell access transistor
		cap = (c_g_ideal_itrs_sram_cell_transistor + 2 * c_fringe_itrs_sram_cell_transistor) * width +
			Lphy_sram_cell_transistor * Cpolywire;
	}
	else if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
		cap = (c_g_ideal_itrs_sram_cell_transistor + 2 * c_fringe_itrs_sram_cell_transistor) * width +
			Lphy_sram_cell_transistor * Cpolywire;
	}
	else{
		cap = (c_g_ideal_itrs_periph_global  + 2 * c_fringe_itrs_periph_global ) * width +
			Lphy_periph_global * Cpolywire;
	}
	return(cap);
}

double gatecappass (double width, double  wirelength) /* returns gate capacitance in Farads */
     /* width: gate width in um (length is Lphy_periph_global) */
     /* wirelength: poly wire length going to gate in lambda */
{
	//v5.0
	double cap;
	if((is_dram)&&(is_access_transistor)){//DRAM cell access transistor
		cap = (c_g_ideal_itrs_dram_access_transistor + 2 * c_fringe_itrs_dram_access_transistor) * width +
			Lphy_dram_access_transistor * Cpolywire;
	}
	else if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
		cap = (c_g_ideal_itrs_dram_wordline_transistor + 2 * c_fringe_itrs_dram_wordline_transistor) * width +
			Lphy_dram_wordline_transistor * Cpolywire;
	}
	else if((!is_dram)&&(is_access_transistor)){//SRAM cell access transistor
		cap = (c_g_ideal_itrs_sram_cell_transistor + 2 * c_fringe_itrs_sram_cell_transistor) * width +
			Lphy_sram_cell_transistor * Cpolywire;
	}
	else if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
		cap = (c_g_ideal_itrs_sram_cell_transistor + 2 * c_fringe_itrs_sram_cell_transistor) * width +
			Lphy_sram_cell_transistor * Cpolywire;
	}
	else{
		cap = (c_g_ideal_itrs_periph_global  + 2 * c_fringe_itrs_periph_global ) * width +
			Lphy_periph_global * Cpolywire;
	}
	return(cap);
}


double 
draincap(double width, int nchannel, int stack, 
				int next_arg_thresh_folding_width_or_height_cell, 
				double fold_dimension)
{
	double width_folded_transistor, ratio_p_to_n, height_transistor_region,
		total_drain_width, total_drain_height_for_cap_wrt_gate,
		drain_cap_area, drain_cap_sidewall, drain_cap_wrt_gate, 
		drain_cap_total, c_junc_area, c_junc_sidewall, c_fringe, c_overlap, 
		drain_height_for_sidewall, drain_cap_metal_connecting_folded_transistors;
	int number_folded_transistors;

	c_junc_sidewall = 0;//MASTAR does not have numbers for sidewall capacitance. 
	if((is_dram)&&(is_access_transistor)){//DRAM cell access transistor
	  c_junc_area = c_junc_itrs_dram_access_transistor;
	  c_fringe = c_fringe_itrs_dram_access_transistor;
	  c_overlap = c_overlap_itrs_dram_access_transistor;
	}
	else
		if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
			 c_junc_area = c_junc_itrs_dram_wordline_transistor;
			 c_fringe = c_fringe_itrs_dram_wordline_transistor;
			 c_overlap = c_overlap_itrs_dram_wordline_transistor;
		}
		else
			if((!is_dram)&&(is_sram_cell)){//SRAM cell access transistor
				c_junc_area = c_junc_itrs_sram_cell_transistor;
				c_fringe = c_fringe_itrs_sram_cell_transistor;
				c_overlap = c_overlap_itrs_sram_cell_transistor;
			}
			else
				if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
					c_junc_area = c_junc_itrs_sram_cell_transistor;
					c_fringe = c_fringe_itrs_sram_cell_transistor;
					c_overlap = c_overlap_itrs_sram_cell_transistor;
				}
				else{//DRAM or SRAM all other transistors
					c_junc_area = c_junc_itrs_periph_global;
					c_fringe = c_fringe_itrs_periph_global;
					c_overlap = c_overlap_itrs_periph_global;
				}
	
	drain_cap_metal_connecting_folded_transistors = 0;
	//Determine the width of the transistor after folding (if it is getting folded)
	if(next_arg_thresh_folding_width_or_height_cell == 0){//interpret fold_dimension as the
		//the folding width threshold i.e. the value of transistor width above which the
		//transistor gets folded
		width_folded_transistor = fold_dimension;
	} 
	else{//interpret fold_dimension as the height of the cell that this transistor is part of. 
		height_transistor_region = height_cell - 2 * HPOWERRAIL;
		ratio_p_to_n = 2 / (2 + 1);
		if(nchannel){
			width_folded_transistor = (1 - ratio_p_to_n) * (height_transistor_region - MIN_GAP_BET_P_AND_N_DIFFS);
		}
		else{
			width_folded_transistor = ratio_p_to_n * (height_transistor_region - MIN_GAP_BET_P_AND_N_DIFFS);
		}
	}
	number_folded_transistors = (int) (ceil(width / width_folded_transistor));

	if(number_folded_transistors < 2){
		width_folded_transistor = width;
	}

	total_drain_width = (WIDTHPOLYCONTACT + 2 * SPACINGPOLYTOCONTACT) + 
		(stack - 1) * SPACINGPOLYTOPOLY;
	drain_height_for_sidewall = width_folded_transistor;
	total_drain_height_for_cap_wrt_gate = width_folded_transistor + 2 * width_folded_transistor * 
		(stack - 1);
	if(number_folded_transistors > 1){
		total_drain_width += (number_folded_transistors - 2) * (WIDTHPOLYCONTACT + 
			2 * SPACINGPOLYTOCONTACT) + (number_folded_transistors - 1) * ((stack - 1) * 
			SPACINGPOLYTOPOLY);
		if(number_folded_transistors%2 == 0){
			drain_height_for_sidewall = 0 ;
		}
		total_drain_height_for_cap_wrt_gate *=  number_folded_transistors;
		drain_cap_metal_connecting_folded_transistors = wire_local_c_per_micron * 
			total_drain_width;
	}
	
	drain_cap_area = c_junc_area * total_drain_width * width_folded_transistor;
	drain_cap_sidewall = c_junc_sidewall * (drain_height_for_sidewall + 2 * total_drain_width);
	drain_cap_wrt_gate  = (c_fringe + c_overlap) * total_drain_height_for_cap_wrt_gate;

	drain_cap_total = drain_cap_area + drain_cap_sidewall + drain_cap_wrt_gate +
		drain_cap_metal_connecting_folded_transistors;
	return(drain_cap_total);
}



/*----------------------------------------------------------------------*/

/* The following routines estimate the effective resistance of an
   on transistor as described in the tech report.  The first routine
   gives the "switching" resistance, and the second gives the 
   "full-on" resistance */

//Currently transresswitch not being used anywhere but may be used in the future.
//double transresswitch (double width,int nchannel,int stack)	/* returns resistance in ohms */
     /* width: in um */
     /* nchannel: whether n or p-channel (boolean) */
     /* stack: number of transistors in series */
//{
  //double restrans;
  //restrans = (nchannel) ? (Rnchannelstatic) : (Rpchannelstatic);
  /* calculate resistance of stack - assume all but switching trans
     have 0.8X the resistance since they are on throughout switching */
  //return ((1.0 + ((stack - 1.0) * 0.8)) * restrans / width);
//}

/*----------------------------------------------------------------------*/

double transreson (double width,int nchannel,int stack)	/* returns resistance in ohms */
     /* width: in um */
     /* nchannel: whether n or p-channel (boolean) */
     /* stack: number of transistors in series */
{
  double restrans;
 

 
	   if((is_dram)&&(is_access_transistor)){//DRAM cell access transistor
		   restrans = (nchannel) ? Rnchannelon_itrs_dram_access_transistor : Rpchannelon_itrs_dram_access_transistor;	    
	   }
	   else
		   if((!is_dram)&&(is_access_transistor)){//SRAM cell access transistor
			    restrans = (nchannel) ? Rnchannelon_itrs_sram_cell_transistor : Rpchannelon_itrs_sram_cell_transistor;	    
		   }
		   else
			   if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
				   restrans = (nchannel) ? Rnchannelon_itrs_sram_cell_transistor : Rpchannelon_itrs_sram_cell_transistor;	    
			   }
			   else
				   if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
					   restrans = (nchannel) ? Rnchannelon_itrs_dram_wordline_transistor : Rpchannelon_itrs_dram_wordline_transistor;	    
				   }
				   else{//DRAM or SRAM all other transistors
					   restrans = (nchannel) ? Rnchannelon_itrs_periph_global : Rpchannelon_itrs_periph_global ;	    
				   }

  /* calculate resistance of stack.  Unlike transres, we don't
     multiply the stacked transistors by 0.8 */
  return (stack * restrans / width);

}

/*----------------------------------------------------------------------*/

/* This routine operates in reverse: given a resistance, it finds
 * the transistor width that would have this R.  It is used in the
 * data wordline to estimate the wordline driver size. */

double restowidth (double res,int nchannel)	/* returns width in um */
     /* res: resistance in ohms */
     /* nchannel: whether N-channel or P-channel */
{
  double restrans;


  //v5.0

	   if((is_dram)&&(is_access_transistor)){//DRAM cell access transistor
		   restrans = (nchannel) ? Rnchannelon_itrs_dram_access_transistor : Rpchannelon_itrs_dram_access_transistor;	    
	   }
	   else
		   if((!is_dram)&&(is_access_transistor)){//SRAM cell access transistor
			    restrans = (nchannel) ? Rnchannelon_itrs_sram_cell_transistor : Rpchannelon_itrs_sram_cell_transistor;	    
		   }
		   else
			   if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
				   restrans = (nchannel) ? Rnchannelon_itrs_dram_wordline_transistor : Rpchannelon_itrs_dram_wordline_transistor;	    
			   }
			   else{//DRAM or SRAM all other transistors
				   restrans = (nchannel) ? Rnchannelon_itrs_periph_global : Rpchannelon_itrs_periph_global ;	    
			   }
  return (restrans / res);

}

/*----------------------------------------------------------------------*/

double horowitz (double inputramptime,double  tf,double  vs1,double  vs2,int rise)
	/* inputramptime: input rise time */
    /* tf: time constant of gate */
    /* vs1, vs2: threshold voltages */
    /* rise: whether INPUT rise or fall (boolean) */
{
  double a, b, td;

  a = inputramptime / tf;
  if (rise == RISE)
    {
      b = 0.5;
      td = tf * sqrt (log (vs1) * log (vs1) + 2 * a * b * (1.0 - vs1)) +
	tf * (log (vs1) - log (vs2));
    }
  else
    {
      b = 0.4;
      td = tf * sqrt (log (1.0 - vs1) * log (1.0 - vs1) + 2 * a * b * (vs1)) +
	tf * (log (1.0 - vs1) - log (1.0 - vs2));
    }
  return (td);
}

//v5.0: Returns log to the base 4 of the input. When the log to the base 4 of a number is
//an integer (example: logbasefour(64) = 3), this function guarantees (within assumed 
//floating point representation error bounds) that the output would be the exact integral
//representation (example: 3 would be represented as 3.0 and not as something like
//2.99999....or 3.000...1). This guaranteed is based on the use of EPSILON2 below.  Note that
//even though the function works for any double number presently its usage at several places 
//implicitly assumes that the function input is a power of 2 integer
double logbasefour (double x)
{
	double temp;
	if (x <= 0){
	  printf ("logbasefour:%g\n", x);
	  exit(1);
	  //FIX. ADD BETTER ERROR HANDLING LATER
	}
   temp = log (x) / log (4.0);
   if(fabs(temp - (int) (temp + EPSILON3)) < EPSILON2){
	   return((int) (temp + EPSILON3));
   }
   else return(temp);
}

//v5.0: Returns log to the base 2 of the input. Read comments above logbasefour function above. 
double logbasetwo (double x)
{
	double temp;
	if (x <= 0){
	  printf ("logbasetwo:%g\n", x);
	  exit(1);
	  //FIX. ADD BETTER ERROR HANDLING LATER
	}
   temp = log (x) / log (2.0);
   if(fabs(temp - (int) (temp + EPSILON3)) < EPSILON2){
	   return((int) (temp + EPSILON3));
   }
   else return(temp);
}


void compute_delay_optimal_repeater(double length, double *opt_sizing,
									  int *opt_number_repeaters, double *delay, 
									  powerDef *power)
{
	double d, R_v, C_g, C_d, b, a, r, k, R_w, C_w, FO4, lcap, wcap, lopt, wopt,
		delay_per_micron, best_delay, perc_diff_from_dyn_energy_best_delay,
		perc_diff_in_delay_from_best_delay_repeater_solution, dyn_energy_best_delay, delay_given_length,
		dyn_energy, best_perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution,
		perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution;
	int index, i, max_index, best_index;
	typedef struct{
		double lcap;
		double wcap;
		double delay;
		powerDef power;
	} repeater_solution;

	repeater_solution solution[10000];

	*opt_sizing = 0;
	*opt_number_repeaters = 0;
	*delay = 0;
	power->readOp.dynamic = 0;
	power->readOp.leakage = 0;
	power->writeOp.dynamic = 0;
	power->writeOp.leakage = 0;

	if(length > 0){
		d = 0;
		R_v = transreson(1, NCH, 1);
		C_g = gatecap(1, 0);
		C_d = draincap(1, NCH, 1, 1, DEFAULTHEIGHTCELL);
		a = 3;
		b = C_d * a / C_g;
		r = 1;
		k = 0.69;
		R_w = wire_outside_mat_r_per_micron;
		C_w = wire_outside_mat_c_per_micron;
		lcap = 1;
		FO4 = R_v * a * (C_d + 4 * C_g);		

		lopt = sqrt((18.9 * d + 2 * k * r * (a + b)) / k) * sqrt(R_v * C_g / (R_w * C_w));
		wopt = sqrt(r / a) * sqrt(R_v * C_w / (R_w * C_g));
		lcap = 1.0;
		wcap = 1.0;
		delay_per_micron = (sqrt(d + (k * r * (a + b) / 9.45)) * sqrt(k / 2) * (1 / lcap + lcap) +
			k * sqrt(a * r / 9.45) * (1 / wcap + wcap)) * sqrt(FO4 * R_w * C_w);
		best_delay = delay_per_micron * length;
		dyn_energy_best_delay = C_w * length * vdd_periph_global * vdd_periph_global * 
			(1 + sqrt(r / a) * sqrt(k / (18.9 * d + 2 * k * r * (a + b))) * (a + b) * wcap / lcap);

		//Now find a solution that is better from an energy-delay perspective. Find the solution that
		//reduces energy by about OPT_PERC_DIFF_IN_ENERGY_FROM_BEST_DELAY_REPEATER_SOLUTION. We first find solutions that are 
		//worse in delay but within MAX_PERC_DIFF_IN_DELAY_FROM_BEST_DELAY_REPEATER_SOLUTION of the best 
		//delay delay. 
		index = 0;
		for(lcap = 1.0; lcap < 2.1; lcap += 0.05){
			for(wcap = 1.0; wcap > 0; wcap -= 0.05){
				delay_per_micron = (sqrt(d + (k * r * (a + b) / 9.45)) * sqrt(k / 2) * (1 / lcap + lcap) +
					k * sqrt(a * r / 9.45) * (1 / wcap + wcap)) * sqrt(FO4 * R_w * C_w);
				delay_given_length = delay_per_micron * length;
				dyn_energy = C_w * length * vdd_periph_global * vdd_periph_global * 
					(1 + sqrt(r / a) * sqrt(k / (18.9 * d + 2 * k * r * (a + b))) * (a + b) * wcap / lcap);
				perc_diff_in_delay_from_best_delay_repeater_solution = (delay_given_length - best_delay) * 100 / best_delay;
				if(perc_diff_in_delay_from_best_delay_repeater_solution < MAX_PERC_DIFF_IN_DELAY_FROM_BEST_DELAY_REPEATER_SOLUTION){
					solution[index].lcap = lcap;
					solution[index].wcap = wcap;
					solution[index].delay = delay_given_length;
					solution[index].power.readOp.dynamic = dyn_energy;
					solution[index].power.readOp.leakage = 0;	
					++index;
				}
			}
		}
		max_index = index - 1;

		best_perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution = BIGNUM;
		for(i=0; i<=max_index; ++i){
			perc_diff_from_dyn_energy_best_delay = (double)(abs((solution[i].power.readOp.dynamic - dyn_energy_best_delay) * 100 / dyn_energy_best_delay));
			perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution = (double) (abs(
				perc_diff_from_dyn_energy_best_delay - 
				(double)(OPT_PERC_DIFF_IN_ENERGY_FROM_BEST_DELAY_REPEATER_SOLUTION)));
			if(perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution < 
				best_perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution){
					best_perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution =
						perc_diff_from_opt_perc_diff_in_energy_from_best_delay_repeater_solution;	
					best_index = i;
            }
        }
        *opt_sizing = solution[best_index].wcap * wopt / (minimum_width_nmos);
        *opt_number_repeaters = length / (solution[best_index].lcap * lopt);
        *delay = solution[best_index].delay;
        power->readOp.dynamic = solution[best_index].power.readOp.dynamic;
		power->readOp.leakage = (length / *opt_number_repeaters) * cmos_ileakage(solution[best_index].wcap * wopt, 
			2 * solution[best_index].wcap * wopt, temper) * 0.5 * vdd_periph_global;
    }
}

