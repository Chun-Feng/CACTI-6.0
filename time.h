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

#ifndef _time
#define _time

#include "areadef.h"
#include "basic_circuit.h"

extern float cumm_per[1024];
extern router_stats_t router_s[ROUTER_TYPES];
//extern double contention[4][7];
extern double FREQUENCY;
extern int cont_stats[2][5][ROUTER_TYPES][7][8];
extern int core_in;
double wire_res(double, double, double);
double wire_cap(double, double, double);
void calc_wire_stats2 (enum wire_type wire_model, wire_stats_t *wire_st) ;
double signal_rise_time (); 
void free_mem (results_mem_array **tag_arr, int t, results_mem_array **data_arr,
            int d, uca_res_lentry_t *ll);
void free_mem2 (nuca_res_lentry_t *ll);
uca_res_lentry_t* sim_uca(uca_org_t *res);
void print_nuca_pda (nuca_org_t *nres);
void output_all (uca_org_t *fr);
void update_min_values (uca_org_t *res, double *delay, double *dyn,
                   double *leak, double *cycle, double *area);
void update_output (uca_org_t *res, uca_org_t *uca_res);
void find_cycle(uca_org_t *res);
void find_area(uca_org_t *res);
void find_acc_time(uca_org_t *res);
void find_power(uca_org_t *res);
void sim_nuca(nuca_org_t *);
void dump_input_args (input_params_t *);
#endif
