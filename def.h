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
#ifndef _def
#define _def
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifdef __linux__
#include <pthread.h>
#include <assert.h>
#endif

/*  The following are things you might want to change
 *  when compiling
 */


/*
 * Address bits in a word, and number of output bits from the cache 
 */

/*
was: #define ADDRESS_BITS 32
now: I'm using 42 bits as in the Power4, 
since that's bigger then the 36 bits on the Pentium 4 
and 40 bits on the Opteron
*/
#define ADDRESS_BITS 42
/*
was: #define BITOUT 64
now: making it a commandline parameter
*/
//static int BITOUT;
int BITOUT;

/*dt: In addition to the tag bits, the tags also include 1 valid bit, 1 dirty bit, 2 bits for a 4-state 
  cache coherency protocoll (MESI), 1 bit for MRU (change this to log(ways) for full LRU). 
  So in total we have 1 + 1 + 2 + 1 = 5 */
#define EXTRA_TAG_BITS 5


/* limits on the various N parameters */


#define MAXDATAN 256       /* Maximum for Ndwl,Ndbl */
#define MAXSUBARRAYS 1024    /* Maximum subarrays for data and tag arrays */
#define MAXDATASPD 256         /* Maximum for Nspd */


double FUDGEFACTOR;

double FEATURESIZE;

/*===================================================================*/


/*
was: 
#define WIRESPACING (2*FEATURESIZE)
#define WIREWIDTH (3*FEATURESIZE)
is: width and pitch are taken from the Intel IEDM 2004 paper on their 65nm process.
*/

//#define WIRESPACING (1.6*FEATURESIZE)
//#define WIREWIDTH (1.6*FEATURESIZE)
/*dt: I've taken the aspect ratio from the Intel paper on their 65nm process */
#define WIREHEIGTHRATIO	(1.8)

//#define WIREPITCH (WIRESPACING+WIREWIDTH)
/*
#define Default_Cmetal 275e-18
*/
/*dt: The old Cmetal was calculated using SiO2 as dielectric (k = 3.9). Going by the Intel 65nm paper, 
low-k dielectics are not at k=2.9. This is a very simple adjustment, as lots of other factors also go into the average
capacitance, but I don't have any better/more up to date data on wire distribution, etc. than what's been done for cacti 1.0.
So I'm doing a simple adjustment by 2.9/3.9 */
//#define Default_Cmetal (2.9/3.9*275e-18)
/*dt: changing this to reflect newer data */
/* 2.9: is the k value for the low-k dielectric used in 65nm Intel process 
   3.9: is the k value of normal SiO2
   --> multiply by 2.9/3.9 to get new C values
   the Intel 130nm process paper mentioned 230fF/mm for M1 through M5 with k = 3.6
   So we get
   230*10^-15/mm * 10^-3 mm/1um * 3.9/3.6
*/
#define Default_Cmetal (2.9/3.9*230e-18*3.9/3.6)
#define Default_Rmetal 48e-3
/* dt: I'm assuming that even with all the process improvements, 
copper will 'only' have 2/3 the sheet resistance of aluminum. */
#define Default_CopperSheetResistancePerMicroM 32e-3

/*dt: this number is calculated from the 2004 ITRS tables (RC delay and Wire sizes)*/
#define CRatiolocal_to_interm  (1.0/1.4)
/*dt: from ITRS 2004 using wire sizes, aspect ratios and effective resistivities for local and intermediate*/
#define RRatiolocal_to_interm  (1.0/2.04)

#define CRatiointerm_to_global  (1.0/1.9)
/*dt: from ITRS 2004 using wire sizes, aspect ratios and effective resistivities for local and intermediate*/
#define RRatiointerm_to_global  (1.0/3.05)
double Cwordmetal;
double Cbitmetal;
double Rwordmetal;
double Rbitmetal;

double TagCwordmetal;
double TagCbitmetal;
double TagRwordmetal;
double TagRbitmetal;



double Cgate;	
double Cgatepass;		

/* note that the value of Cgatepass will be different depending on 
   whether or not the source and drain are at different potentials or
   the same potential.  The two values were averaged */

/* fF/um */
//#define Cpolywire	(0.25e-15)	
double Cpolywire;			 



//static double vdd_periph_global;
/* Threshold voltages (as a proportion of Vdd)
   If you don't know them, set all values to 0.5 */

#define SizingRatio   0.33
#define VTHNAND       0.561
#define VTHFA1        0.452
#define VTHFA2        0.304
#define VTHFA3        0.420
#define VTHFA4        0.413
#define VTHFA5        0.405
#define VTHFA6        0.452
#define VSINV         0.452   
#define VTHINV100x60  0.438   /* inverter with p=100,n=60 */
#define VTHINV360x240 0.420   /* inverter with p=360, n=240 */
#define VTHNAND60x90  0.561   /* nand with p=60 and three n=90 */
#define VTHNOR12x4x1  0.503   /* nor with p=12, n=4, 1 input */
#define VTHNOR12x4x2  0.452   /* nor with p=12, n=4, 2 inputs */
#define VTHNOR12x4x3  0.417   /* nor with p=12, n=4, 3 inputs */
#define VTHNOR12x4x4  0.390   /* nor with p=12, n=4, 4 inputs */
#define VTHOUTDRINV    0.437
#define VTHOUTDRNOR   0.379
#define VTHOUTDRNAND  0.63
#define VTHOUTDRIVE   0.425
#define VTHCOMPINV    0.437
#define VTHMUXNAND    0.548
#define VTHMUXDRV1    0.406
#define VTHMUXDRV2    0.334
#define VTHMUXDRV3    0.478
#define VTHEVALINV    0.452
#define VTHSENSEEXTDRV  0.438

#define VTHNAND60x120 0.522


double Gm_sense_amp_transistors;



double Wmemcella_dram, Wmemcella_sram;
double Wmemcellpmos_dram, Wmemcellpmos_sram;
double Wmemcellnmos_dram, Wmemcellnmos_sram;

//#define Wpchmax		(25.0) /* precharge transistor sizes usually do not exceed 25 */
double Wpchmax;


double Wiso;
double WsenseEn;
double WsenseN;
double WsenseP;


double width_nmos_bit_mux, width_nmos_sense_amp_mux, width_pmos_bitline_precharge, 
width_pmos_bitline_equalization;


//#define Wcompinvp1	(10.0)
double Wcompinvp1;
//#define Wcompinvn1	(6.0)
double Wcompinvn1;
//#define Wcompinvp2	(20.0)
double Wcompinvp2;
//#define Wcompinvn2	(12.0)
double Wcompinvn2;
//#define Wcompinvp3	(40.0)
double Wcompinvp3;
//#define Wcompinvn3	(24.0)
double Wcompinvn3;
//#define Wevalinvp	(80.0)
double Wevalinvp;
//#define Wevalinvn	(40.0)
double Wevalinvn;



//#define Wcompn		(10.0)
double Wcompn;
//#define Wcompp		(30.0)
double Wcompp;



double WmuxdrvNANDn    ;
double WmuxdrvNANDp    ;






double Vbitpre_dram;
double Vt_dram;
double Vbitpre_sram;
double Vt_sram;
/*
was: #define Vbitsense	(0.10)
now: 50mV seems to be the norm as of 2005
*/
//v5.0
//#define Vbitsense	(0.05*Vdd)
double Vbitsense;
#define Vbitswing	(0.20*Vdd)



/*===================================================================*/

/*
 * The following are things you probably wouldn't want to change.  
 */


#define TRUE 1
#define FALSE 0
#ifndef NULL
#define NULL 0
#endif
#define OK 1
#define ERROR 0
#define BIGNUM 1e30
#define INF 9999999
#define DIVIDE(a,b) ((b)==0)? 0:(a)/(b)
#define MAX(a,b) (((a)>(b))?(a):(b))
//v5.0
#define MIN(a,b) (((a)<(b))?(a):(b))

#define WAVE_PIPE 3
#define MAX_COL_MUX 128

/* Used to communicate with the horowitz model */

#define RISE 1
#define FALL 0
#define NCH  1
#define PCH  0


/* Used to pass values around the program */

/*dt: maximum numbers of entries in the 
      caching structures of the tag calculations
*/
#define MAX_CACHE_ENTRIES 512

//static int sequential_access_flag;
//static int fast_cache_access_flag;
int sequential_access_flag;
int fast_cache_access_flag;
int pure_sram_flag; //Changed from static int to just int as value wasn't getting passed through to 
//area function in area.c
int is_dram;//v5.0

#define EPSILON 0.5 //v4.1: This constant is being used in order to fix floating point -> integer
//conversion problems that were occuring within CACTI. Typical problem that was occuring was
//that with different compilers a floating point number like 3.0 would get represented as either 
//2.9999....or 3.00000001 and then the integer part of the floating point number (3.0) would 
//be computed differently depending on the compiler. What we are doing now is to replace 
//int (x) with (int) (x+EPSILON) where EPSILON is 0.5. This would fix such problems. Note that
//this works only when x is an integer >= 0. 

#define EPSILON2 0.1  

#define EPSILON3 0.6

//v5.0: The following constants define how small/large a single subarray can be. These
//values of these constants have been chosen somewhat arbitrarily based on expected ranges
//of cache sizes, minimum and maximum bounds of the circuit choices. 

#define MINSUBARRAYROWS 16 //At least one 3-to-8 predecoder. Note that it probably
////starts to make sense to replace the 3-i/p NAND gate in the predecoder circuit with
////a 2-i/p NAND gate in order to alleviate the effects of second-order effects that
////are not being modeled.
#define MAXSUBARRAYROWS 262144 //At most max number of inputs of decode NOR is fixed to be 4.
//Note that second-order effects start playing a role even for a 4-i/p static CMOS NOR gate 
//of the kind that we have chosen at today's and future technologies. Need to replace this
//circuit with something else.
#define MINSUBARRAYCOLS 2
#define MAXSUBARRAYCOLS 262144


double MIN_GAP_BET_P_AND_N_DIFFS, HPOWERRAIL, DEFAULTHEIGHTCELL, MIN_GAP_BET_SAME_TYPE_DIFFS;
double WIDTHPOLYCONTACT, SPACINGPOLYTOPOLY, SPACINGPOLYTOCONTACT;


#define INV 0
#define NOR 1
#define NAND 2

double Wcolmuxdec3to8p, Wcolmuxdec3to8n, Wcolmuxnorp, Wcolmuxnorn;

double c_dram_cell;

double restore_delay;
double refresh_power;
double dram_refresh_period;

double input_tech;



#define NUMBER_TECH_FLAVORS 4
double vdd[NUMBER_TECH_FLAVORS];
double Lphy[NUMBER_TECH_FLAVORS];
double Lelec[NUMBER_TECH_FLAVORS];
double t_ox[NUMBER_TECH_FLAVORS];
double v_th[NUMBER_TECH_FLAVORS];
double c_ox[NUMBER_TECH_FLAVORS];
double mobility_eff[NUMBER_TECH_FLAVORS];
double Vdsat[NUMBER_TECH_FLAVORS];
double c_g_ideal[NUMBER_TECH_FLAVORS];
double c_fringe[NUMBER_TECH_FLAVORS];
double c_junc[NUMBER_TECH_FLAVORS];
double I_on_n[NUMBER_TECH_FLAVORS];
double I_on_p[NUMBER_TECH_FLAVORS];
double Rnchannelon[NUMBER_TECH_FLAVORS];
double Rpchannelon[NUMBER_TECH_FLAVORS];
double I_off_n[NUMBER_TECH_FLAVORS][200];
double I_off_p[NUMBER_TECH_FLAVORS][200];

#define NUMBER_INTERCONNECT_PROJECTION_TYPES 2 //aggressive and conservative
//0 = Ron Ho aggressive projections, 1 = Ron Ho conservative projections
#define NUMBER_WIRE_TYPES 6 //local, semi-global and global
//1 = Ron Ho semi-global wire , 2 = Ron Ho global wire, 2 = Local wire: similar to Ron Ho semi-global wire except with pitch of 2*FEATURESIZE


double wire_inside_mat_pitch, wire_inside_mat_r_per_micron, wire_inside_mat_c_per_micron,
wire_outside_mat_pitch, wire_outside_mat_r_per_micron, wire_outside_mat_c_per_micron,
wire_local_pitch, wire_local_r_per_micron, wire_local_c_per_micron;

double wire_pitch[NUMBER_INTERCONNECT_PROJECTION_TYPES][NUMBER_WIRE_TYPES],
wire_r_per_micron[NUMBER_INTERCONNECT_PROJECTION_TYPES][NUMBER_WIRE_TYPES],
wire_c_per_micron[NUMBER_INTERCONNECT_PROJECTION_TYPES][NUMBER_WIRE_TYPES];
int interconnect_projection_type;//whether aggressive or conservative
int wire_inside_mat_type, wire_outside_mat_type; //whether semi-global or global

double wire_local_pitch_tech_node[2], wire_local_r_per_micron_tech_node[2],
wire_local_c_per_micron_tech_node[2], wire_inside_mat_pitch_tech_node[2],
wire_inside_mat_r_per_micron_tech_node[2], wire_inside_mat_c_per_micron_tech_node[2],
wire_outside_mat_pitch_tech_node[2], wire_outside_mat_r_per_micron_tech_node[2],
wire_outside_mat_c_per_micron_tech_node[2];


int periph_global_tech_flavor;
int sram_cell_and_wordline_tech_flavor;
int dram_cell_tech_flavor;
double Cu_resistivity;
	  

double vdd_periph_global;
double Lphy_periph_global;
double Lelec_periph_global;
double t_ox_periph_global;
double v_th_periph_global;
double c_ox_periph_global;
double mobility_eff_periph_global;
double Vdsat_periph_global;
double c_g_ideal_itrs_periph_global ;
double c_fringe_itrs_periph_global ;
double c_junc_itrs_periph_global ;
double c_overlap_itrs_periph_global ;
double I_on_n_periph_global ;
double I_on_p_periph_global ;
double Rnchannelon_itrs_periph_global ;
double Rpchannelon_itrs_periph_global ;
double  I_off_n_periph_global[200];
double  I_off_p_periph_global[200];

double vdd_sram_cell;
double Lphy_sram_cell_transistor;
double t_ox_sram_cell_transistor;
double v_th_sram_cell_transistor;

double Lmemcella;
double Lmemcellpmos;
double Lmemcellnmos;
double Lphy_sram_cell_transistor;
double Lelec_sram_cell_transistor;
double area_cell_dram;
double asp_ratio_cell_dram; 
double BitWidth_dram;
double BitHeight_dram;
double area_cell_sram;
double asp_ratio_cell_sram; 
double BitWidth_sram;
double BitHeight_sram;
double c_g_ideal_itrs_sram_cell_transistor;
double c_fringe_itrs_sram_cell_transisor;
double c_junc_itrs_sram_cell_transistor;
double c_fringe_itrs_sram_cell_transistor;
double c_overlap_itrs_sram_cell_transistor;
double I_on_n_sram_cell_transistor;
double I_on_p_sram_cell_transistor;
double Rnchannelon_itrs_sram_cell_transistor;
double Rpchannelon_itrs_sram_cell_transistor;
double I_off_n_sram_cell_transistor[200];
double I_off_p_sram_cell_transistor[200];
double t_ox_lstp;
double v_th_lstp;
double vth_lstp_sram_cell;

double vdd_dram_cell;
double vpp;
double v_th_dram_access_transistor;
double Lphy_dram_wordline_transistor;
double Lelec_dram_wordline_transistor;
double width_dram_access_transistor;
double c_g_ideal_itrs_dram_wordline_transistor;
double c_fringe_itrs_dram_wordline_transistor;
double c_junc_itrs_dram_wordline_transistor;
double c_fringe_itrs_dram_wordline_transistor;
double c_overlap_itrs_dram_wordline_transistor;
double I_on_n_dram_wordline_transistor;
double I_on_p_dram_wordline_transistor;
double Rnchannelon_itrs_dram_wordline_transistor;
double Rpchannelon_itrs_dram_wordline_transistor;
double I_off_n_dram_wordline_transistor[200];
double I_off_p_dram_wordline_transistor[200];
double Lphy_dram_access_transistor;
double Lelec_dram_access_transistor;
double c_g_ideal_itrs_dram_access_transistor;
double c_fringe_itrs_dram_access_transistor;
double c_junc_itrs_dram_access_transistor;
double c_fringe_itrs_dram_access_transistor;
double c_overlap_itrs_dram_access_transistor;
double I_on_n_dram_access_transistor;
double I_on_p_dram_access_transistor;
double Rnchannelon_itrs_dram_access_transistor;
double Rpchannelon_itrs_dram_access_transistor;
int is_wordline_transistor;
int is_access_transistor;
int is_sram_cell;
double I_on_dram_cell;
double  I_off_dram_cell_worst_case_length_temp;
double gmn_sense_amp_latch;
double gmp_sense_amp_latch;
double Gm_sense_amp_latch;

double vdd_periph_global_tech_node[2], t_ox_periph_global_tech_node[2],
v_th_periph_global_tech_node[2], c_ox_periph_global_tech_node[2],
c_g_ideal_itrs_periph_global_tech_node[2], c_fringe_itrs_periph_global_tech_node[2],
c_junc_itrs_periph_global_tech_node[2], 
Lphy_periph_global_tech_node[2], Lelec_periph_global_tech_node[2], I_on_n_periph_global_tech_node[2],
I_off_n_periph_global_tech_node[2][200], I_off_p_periph_global_tech_node[2][200]; 

double vdd_sram_cell_tech_node[2], Lphy_sram_cell_transistor_tech_node[2], Lelec_sram_cell_transistor_tech_node[2],
t_ox_sram_cell_transistor_tech_node[2], v_th_sram_cell_transistor_tech_node[2], c_g_ideal_itrs_sram_cell_transistor_tech_node[2],
c_fringe_itrs_sram_cell_transistor_tech_node[2], c_junc_itrs_sram_cell_transistor_tech_node[2],
I_on_n_sram_cell_transistor_tech_node[2],
I_off_n_sram_cell_transistor_tech_node[2][200], I_off_p_sram_cell_transistor_tech_node[2][200];

double  vdd_dram_cell_tech_node[2], v_th_dram_access_transistor_tech_node[2], Lphy_dram_access_transistor_tech_node[2], 
Lelec_dram_access_transistor_tech_node[2], c_g_ideal_itrs_dram_access_transistor_tech_node[2],
c_fringe_itrs_dram_access_transistor_tech_node[2], c_junc_itrs_dram_access_transistor_tech_node[2],
I_on_n_dram_access_transistor_tech_node[2], c_dram_cell_tech_node[2],
vpp_tech_node[2], Lphy_dram_wordline_transistor_tech_node[2], Lelec_dram_wordline_transistor_tech_node[2], c_g_ideal_itrs_dram_wordline_transistor_tech_node[2],
c_fringe_itrs_dram_wordline_transistor_tech_node[2], c_junc_itrs_dram_wordline_transistor_tech_node[2],
I_on_n_dram_wordline_transistor_tech_node[2],
I_off_n_dram_wordline_transistor_tech_node[2][200], I_off_p_dram_wordline_transistor_tech_node[2][200],
 mobility_eff_periph_global_tech_node[2], Vdsat_periph_global_tech_node[2];

double Wmemcella_dram_tech_node[2], Wmemcellpmos_dram_tech_node[2], Wmemcellnmos_dram_tech_node[2],
 area_cell_dram_tech_node[2], asp_ratio_cell_dram_tech_node[2];
double Wmemcella_sram_tech_node[2], Wmemcellpmos_sram_tech_node[2], Wmemcellnmos_sram_tech_node[2],
 area_cell_sram_tech_node[2], asp_ratio_cell_sram_tech_node[2];


double height_cell;
double width_cell;


/* 
 * transistor sizes related to tag array in fully associative
 * cache
 */

#define Wfadecdrivep (125*FEATURESIZE)
#define Wfadecdriven (62.5*FEATURESIZE)
#define Wdecdrivep (450*FEATURESIZE)
#define Wdecdriven (300*FEATURESIZE)
#define Wfadrivep  (125*FEATURESIZE)
#define Wfadriven  (62.5*FEATURESIZE)
#define Wfadrive2p (500*FEATURESIZE)
#define Wfadrive2n (250*FEATURESIZE)
#define Wfadecdrive1p (12.5*FEATURESIZE)
#define Wfadecdrive1n (6.25*FEATURESIZE)
#define Wfadecdrive2p (50*FEATURESIZE)
#define Wfadecdrive2n (25*FEATURESIZE)
#define Wfaprechn  (7.5*FEATURESIZE)
#define Wfaprechp  (12.5*FEATURESIZE)
#define Wdummyn    (12.5*FEATURESIZE)
#define Waddrnandn (62.5*FEATURESIZE)
#define Waddrnandp (62.5*FEATURESIZE)
#define Wdummyinvn (75*FEATURESIZE)
#define Wdummyinvp (100*FEATURESIZE)
#define Wfainvn (12.5*FEATURESIZE)
#define Wfainvp (25*FEATURESIZE)
#define Wfanorn (6.25*FEATURESIZE)
#define Wfanorp (12.5*FEATURESIZE)
#define Wdecinvn (20*FEATURESIZE)
#define Wdecinvp (40*FEATURESIZE)

extern int NUCA;
extern double SENSE_AMP_D, SENSE_AMP_P;
#define VOL_SWING 0.100 /* 100 mv - low swing voltage*/
#define test_components
#define GLOBAL_WIRES    1
#define BANDWIDTH 128 // Grid network link bandwidth
//#define ENABLE_THREAD /* CACTI in multi-threaded mode */
//#define DEBUG
//#define DEBUGWIRE
#define ROUTER_LAT 3
#define FIXED_OVERHEAD 55e-12 /* clock skew and jitter in s. Ref: Hrishikesh et al ISCA 01 */
#define LATCH_DELAY 28e-12 /* latch delay in s (later should use FO4 TODO) */
//#define WIREDELAY   0.145812877893269 /* ns/mm FIXME */
/* in FO4  Ref: Hrishikesh et al ISCA 01*/
#define LATCH_FO4 1 
#define WIRE_TYPES 6
#define ROUTER_TYPES 3

#ifdef DEBUG
#define PRINTD(a);\
a;
#else
#define PRINTD(a);\

#endif

#ifdef DEBUGWIRE
#define PRINTDW(a);\
a;
#else
#define PRINTDW(a);\

#endif

#define ERROR_EXIT(a,b);\
fprintf(a,b);\
exit(-1);

#define VBITSENSEMIN 0.08 //minimum bitline sense voltage is fixed to be 80 mV.


int temper;

double height_mat;
double width_mat;
double height_subarray;
double width_subarray;

#define fopt 4.0

double minimum_width_nmos;
double minimum_width_pmos;

double gnand2;
double gnor2;
double gnand3;
double gpmos;
double gnmos;

double MIN_PERCENT_WITHIN_BEST_AREA;
double MIN_PERCENT_WITHIN_BEST_DELAY;

#define INPUT_WIRE_TO_INPUT_GATE_CAP_RATIO 0
double kinv;
double MAX_NMOS_WIDTH;
double MAX_PMOS_WIDTH;

#define MAX_PERC_DIFF_IN_DELAY_FROM_BEST_DELAY_REPEATER_SOLUTION 30
#define OPT_PERC_DIFF_IN_ENERGY_FROM_BEST_DELAY_REPEATER_SOLUTION 50

#endif
