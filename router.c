/*-----------------------------------------------------------------------------
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
 * Users of this software agree to the terms and conditions set forth herein, 
 * and hereby grant back to Hewlett-Packard Company and its affiliated 
 * companies ("HP") a non-exclusive, unrestricted, royalty-free right and 
 * license under any changes, enhancements or extensions  made to the core 
 * functions of the software, including but not limited to those affording 
 * compatibility with other hardware or software environments, but excluding 
 * applications which incorporate this software. Users further agree to use 
 * their best efforts to return to HP any such changes, enhancements or 
 * extensions that they make and inform HP of noteworthy uses of this software.
 * Correspondence should be provided to HP at:
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
 *----------------------------------------------------------------------------*/

#include "router.h"
#define VTHNUCA 0

//long double Buffer_entries = 16; //Number of entries in the buffer
//long double Flit_size = 128; //Flit or word size that is equal to the 
//bandwidth of the interconnect (change W along with this)
long double Pr = 1; //Number of read port
long double Pw = 1; //Number of write port
long double  L; //Process technology
long double Vdd;// = vdd[VTHNUCA]; 
//long double vc_count = 2;

/* Technology independent parameters */
/*Buffer parameters*/
long double Cpoly = 1.95e-3; //poly capacitance per unit area
long double NCdiff_area = 1.37e-4; //diffusion area capacitance per unit area for p-mos 
long double PCdiff_area = 3.43e-4; //cap for n-mos
long double NCdiff_side = 2.75e-10; //diffusion side capacitance per unit length for p-mos
long double PCdiff_side = 2.75e-10; //cap for n-mos
long double NCdiff_ovlp = 4.01e-10; //diffusion overlap capacitance per unit width for p-mos
long double PCdiff_ovlp = 4.76e-10; //cap for p-mos F/m

long double NTm;   //(m)
long double PTm;   //(m)
long double Tpw;   //(m)
long double Tpr;   //(m)
long double dw;    //(m)
long double hcell; //(m)
long double wcell; //(m)

long double c_w3; // Capacitance of a wire with thrice the minimal spacing (F/m)
long double NTwd; // (m)
long double PTwd; // (m)
long double NTbd; // (m)
long double PTbd; // (m)
long double Tc;   // (m)

/*Crossbar parameters. Transmisson gate is employed for connector*/
long double NTtr; /* Transmission gate's nmos tr. length */
long double PTtr; /* pmos tr. length */
long double wt; /* track width */
long double ht; /* track height */
long double I; /* Number of crossbar input ports */
long double O; /* Number of crossbar output ports */
long double W; /* Crossbar port width in bits */
long double NTi;
long double PTi;

long double NTid; // (m)
long double PTid; // (m)
long double NTod; // (m)
long double PTod; // (m)

/*Arbiter parameters*/
long double R; /*number of requests*/
long double NTn1;
long double PTn1;
long double NTn2;
long double PTn2;


/* 
 * v6.0: Router related functions.
 */

void
free_router() {
    int i;
    for(i=0; i<ROUTER_TYPES; i++) {
        free(router_s[i].ubuffer.params);
        free(router_s[i].ubuffer.data_array);
    }
}

void
init_router_params() {
    long double technology = (long double)FEATURESIZE; /* in u */
    int ro;
    int i, j, k, l, m;
    FILE *cont;
    char line[5000];
    char jk[5000];
    input_params_t *rparams[ROUTER_TYPES];
    results_mem_array *buff_data[ROUTER_TYPES];
    for(i=0; i<ROUTER_TYPES; i++) {
        rparams[i] = (input_params_t *) malloc(sizeof(input_params_t));
        buff_data[i] = (results_mem_array *) malloc(sizeof(results_mem_array));
        bzero(buff_data[i], sizeof(results_mem_array));
        bzero(rparams[i], sizeof(input_params_t));
        router_s[i].ubuffer.params = rparams[i];
        router_s[i].ubuffer.data_array = buff_data[i];
    }
    cont = fopen("contention.dat", "r");
   
    for(i=0; i<2; i++) { 
        for(j=2; j<5; j++) {
            for(k=0; k<ROUTER_TYPES; k++) {
                for(l=0;l<7; l++) {
                    int *temp = cont_stats[i/*l2 or l3*/][j/*core*/][k/*64 or 128 or 256 link bw*/][l /* no banks*/];
                    assert(fscanf(cont, "%[^\n]\n", line) != EOF);
                    sscanf(line, "%[^:]: %d %d %d %d %d %d %d %d",jk, &temp[0], &temp[1], &temp[2], &temp[3], 
                            &temp[4], &temp[5], &temp[6], &temp[7]); 
//                    printf("%d %d %d %d %d %d %d %d\n", cont_stats[i][j][k][l][0], 
//                            cont_stats[i][j][k][l][1],
//                            cont_stats[i][j][k][l][2], cont_stats[i][j][k][l][3], 
//                            cont_stats[i][j][k][l][4], cont_stats[i][j][k][l][5],
//                            cont_stats[i][j][k][l][6], cont_stats[i][j][k][l][7]);
                }
            }
        }
    }

    Vdd = .9;
    L = technology*1e-6; //Process technology
    NTm = 12*technology*1e-6/2; //(m)
    PTm = 6*technology*1e-6/2; //(m)
    Tpw = 5*technology*1e-6/2; //(m)
    Tpr = 10*technology*1e-6/2; //(m)
    dw = 15*technology*1e-6/2; //(m)
    hcell = 40*technology*1e-6/2; //(m)
    wcell = 20*technology*1e-6/2; //(m)

    c_w3 = (long double) wire_cap(1, 1, 3);
//    c_w3 = 0.123e-9; // Capacitance of a wire with thrice the minimal spacing (F/m)
    NTwd = 60*technology*1e-6/2; // (m)
    PTwd = 120*technology*1e-6/2; // (m)
    NTbd = 60*technology*1e-6/2; // (m)
    PTbd = 120*technology*1e-6/2; // (m)
    Tc = 60*technology*1e-6/2; // (m)

    /*Crossbar parameters. Transmisson gate is employed for connector*/
    NTtr = 10*technology*1e-6/2; /*Transmission gate's nmos tr. length*/
    PTtr = 20*technology*1e-6/2; /* pmos tr. length*/
    wt = 15*technology*1e-6/2; /*track width*/
    ht = 15*technology*1e-6/2; /*track height*/
    I = 5; /*Number of crossbar input ports*/
    O = 5; /*Number of crossbar output ports*/
//    W = Flit_size; /*Crossbar port width in bits*/
    NTi = 12.5*technology*1e-6/2;
    PTi = 25*technology*1e-6/2;

    NTid = 60*technology*1e-6/2; //m
    PTid = 120*technology*1e-6/2; // m
    NTod = 60*technology*1e-6/2; // m
    PTod = 120*technology*1e-6/2; // m

    /*Arbiter parameters*/
    R=5; /*number of requests*/
    NTn1 = 13.5*technology*1e-6/2;
    PTn1 = 76*technology*1e-6/2;
    NTn2 = 13.5*technology*1e-6/2;
    PTn2 = 76*technology*1e-6/2;

    for (ro=0; ro<ROUTER_TYPES; ro++) {
        if (ro == 0) {
            //Flit or word size that is equal to the 
            router_s[ro].flit_size = 64;
            //Number of entries in the buffer
            router_s[ro].buffer_ent = 8;
            W = router_s[ro].flit_size; /*Crossbar port width in bits*/
            router_s[ro].vc_count = 4;
        }
        else if (ro == 1) {
            //Flit or word size that is equal to the 
            router_s[ro].flit_size = 128;
            //Number of entries in the buffer
            router_s[ro].buffer_ent = 4;
            W = router_s[ro].flit_size; /*Crossbar port width in bits*/
            router_s[ro].vc_count = 4;
        }
        else if (ro == 2) {
            //Flit or word size that is equal to the 
            router_s[ro].flit_size = 256;
            //Number of entries in the buffer
            router_s[ro].buffer_ent = 2;
            W = router_s[ro].flit_size; /*Crossbar port width in bits*/
            router_s[ro].vc_count = 4;
        }
        calc_router_parameters2 (&(router_s[ro]));
    }

}

long double //wire cap with triple spacing
Cw3(long double length) {
    return (c_w3*length);
}

/*Function to calculate the gate capacitance*/
long double 
gate_cap(long double w) {
    return (long double) gatecap (w*1e6 /*u*/, 0);
}

/*Function to calculate the diffusion capacitance*/
long double 
diff_cap(long double w, int type /*0 for n-mos and 1 for p-mos*/,
                            long double s /*number of stacking transistors*/) {
    return (long double) draincap (w*1e6 /*u*/, type, s, 1, DEFAULTHEIGHTCELL);
}

/*Function to calculate the transistor resistance for a given load 
 * and rise time*/
long double 
tr_res(long double cload /*load capacitance*/,
                         long double rise_t /*rise time*/) {
    return (-(.2e-9*rise_t)/(cload * (long double)log (1-.9)));
}

/*Function that returns the transistor size for a given resistance*/
long double 
tr_size(long double res, long double type) {
    return   (long double) restowidth (res, type)*1e-6;
}

/*crossbar related functions */
// Model for simple transmission gate
long double 
transmission_buf_inpcap() {
    return diff_cap(NTtr, 0, 1)+diff_cap(PTtr, 1, 1);
}

long double 
transmission_buf_outcap() {
    return diff_cap(NTtr, 0, 1)+diff_cap(PTtr, 1, 1);
}

long double 
transmission_buf_ctrcap() {
    return gate_cap(NTtr)+gate_cap(PTtr);
}

long double 
crossbar_inpline(router_stats_t *r) {
    return (Cw3(O*r->flit_size*wt) + O*transmission_buf_inpcap() + gate_cap(NTid) + 
            gate_cap(PTid) + diff_cap(NTid, 0, 1) + diff_cap(PTid, 1, 1));
}

long double 
crossbar_outline(router_stats_t *r) {
    return (Cw3(I*r->flit_size*ht) + I*transmission_buf_outcap() + gate_cap(NTod) + 
           gate_cap(PTod) + diff_cap(NTod, 0, 1) + diff_cap(PTod, 1, 1));
}

long double 
crossbar_ctrline(router_stats_t *r) {
    return (Cw3(0.5*O*r->flit_size*wt) + W*transmission_buf_ctrcap() + 
           diff_cap(NTi, 0, 1) + diff_cap(PTi, 1, 1) + 
           gate_cap(NTi) + gate_cap(PTi));
}

long double 
crossbar_power(router_stats_t *r) {
    return (crossbar_inpline(r)*Vdd*Vdd*W/2 + 
            crossbar_outline(r)*Vdd*Vdd*W/2)*2;
}

/*arbiter related functions*/
long double 
arb_req() {
    return ((R-1)*(2*gate_cap(NTn1)+gate_cap(PTn1)) + 2*gate_cap(NTn2) + 
        gate_cap(PTn2) + gate_cap(NTi) + gate_cap(PTi) + 
        diff_cap(NTi, 0, 1) + diff_cap(PTi, 1, 1));
}

long double 
arb_pri() {
    return 2*(2*gate_cap(NTn1)+gate_cap(PTn1)); /* switching capacitance 
    of flip-flop is ignored */
}


long double 
arb_grant(router_stats_t *r) {
    return diff_cap(NTn1, 0, 1)*2 + diff_cap(PTn1, 1, 1) + crossbar_ctrline(r);
}

long double 
arb_int() {
    return (diff_cap(NTn1, 0, 1)*2 + diff_cap(PTn1, 1, 1) + 
            2*gate_cap(NTn2) + gate_cap(PTn2));
}

long double 
arb_power(router_stats_t *r,int n_req /* no. of requests */) {
    return (n_req*arb_req()*Vdd*Vdd/2 + n_req*arb_pri()*Vdd*Vdd/2 + 
           arb_grant(r)*Vdd*Vdd + arb_int()*0.5*Vdd*Vdd);
}

void
buffer_params( router_stats_t *r) 
{
    input_params_t *params=r->ubuffer.params;
    params->cache_size = r->flit_size*r->vc_count*r->buffer_ent;
    params->block_size = r->flit_size/8;
    params->associativity = 1;
    params->rw_ports = 0;
    params->excl_read_ports = 5;
    params->excl_write_ports = 5;
    params->single_ended_read_ports = 0;
    params->uca_banks = 1;
    params->tech_size = FEATURESIZE;
    params->output_width = r->flit_size*2;
    params->force_tag = 0;
    params->tag_size = 0;
    params->access_mode = 2;
    params->pure_sram = 1;
    params->dram = 0;
    params->delay_wt = 100;
    params->dynamic_power_wt = 100;
    params->leakage_power_wt = 0;
    params->cycle_time_wt = 100;
    params->area_wt = 100;
    params->area_dev = 1000;
    params->delay_dev = 1000;
    params->dynamic_power_dev = 1000;
    params->leakage_power_dev = 1000;
    params->cycle_time_dev = 1000;
    params->force_wiretype = 1;
    params->wire_inter_mats = 1;
    params->nuca = 0;
    params->data_associativity = 1;
    params->sequential_access = 0;
    params->print_detail = 1;

}

void
buffer_stats(router_stats_t *r)
{

    /* get buffer parameters */
    buffer_params (r);

    sim_uca(&(r->ubuffer)); /* simulate buffer */
    r->simple_buffer_read = read_operation(r);
    r->simple_buffer_write = write_operation(r);
    return;
}

void
cb_stats (router_stats_t *r)
{
    r->crossbar.power.dynamic = crossbar_power(r);
    r->crossbar.power.leakage = r->flit_size * 5 * 5 * cmos_ileakage(60*minimum_width_nmos, 
        120*minimum_width_pmos, temper);
}

void
arb_stats (router_stats_t *r)
{
    r->arbiter.power.dynamic = (double) arb_power(r, 5) +
                                (double) arb_power(r, 5* r->vc_count);
    r->arbiter.power.leakage = 0;
}

void
get_router_power(router_stats_t *r)
{
    double M = .6; //network load
    /* calculate buffer stats */
    buffer_stats(r);

    /* calculate cross-bar stats */
    cb_stats(r);

    /* calculate arbiter stats */
    arb_stats(r);
    r->router_pda.power.dynamic = 
    (5.0*((double)read_operation(r) + (double)write_operation(r)) + 
//    r->buffer.total_power.readOp.dynamic +
    r->crossbar.power.dynamic +
    r->arbiter.power.dynamic)*M;
}

void
get_router_delay (router_stats_t *rs)
{
    rs->cycle_time = (1/(double)FREQUENCY)*1e3; //ps
    rs->router_pda.delay = 3;
    rs->max_cyc = 17 * FO4; //ps
    if (rs->cycle_time < rs->max_cyc) {
        FREQUENCY = (1/rs->max_cyc)*1e3; //GHz
    }
}

void
get_router_area(router_stats_t *rs)
{
    rs->router_pda.area_stats.height = 
        rs->ubuffer.data_array->bank_height;
    rs->router_pda.area_stats.width =
        rs->ubuffer.data_array->bank_width;
    rs->router_pda.area_stats.area = 
        rs->ubuffer.data_array->area;
}

void
calc_router_parameters2(router_stats_t *rs)
{
    /* calculate router frequency and pipeline cycles */
    get_router_delay(rs);

    /* router power stats */
    get_router_power(rs);

    /* area stats */
    get_router_area(rs);
}

void
print_router(router_stats_t *r)
{
    fprintf(stderr, "\n\nRouter stats:\n");
    fprintf(stderr, "\tMaximum possible network frequency - %g GHz\n",
                    (1/r->max_cyc)*1e3);
    fprintf(stderr, "\tNetwork frequency - %g GHz\n",
                    FREQUENCY);
    fprintf(stderr, "\tNo. of Virtual channels - %d\n", r->vc_count);
    fprintf(stderr, "\tNo. of pipeline stages - %g\n", r->router_pda.delay);
    fprintf(stderr, "\tLink bandwidth - %d (bits)\n", r->flit_size);
    fprintf(stderr, "\tNo. of buffer entries per virtual channel - %d\n", 
            r->buffer_ent);
//    fprintf(stderr, "\tBuffer access (read) - %g (nJ)\n", 
//                    r->buffer.total_power.readOp.dynamic * 1e9);
    fprintf(stderr, "\tSimple buffer access (read) - %g (nJ)\n", 
            r->simple_buffer_read * 1e9);
//    fprintf(stderr, "\tBuffer access (write) - %g (nJ)\n", 
//                    r->buffer.total_power.readOp.dynamic * 1e9);
    fprintf(stderr, "\tSimple buffer access (write) - %g (nJ)\n", 
            r->simple_buffer_write * 1e9);
//    fprintf(stderr, "\tBuffer access (read) - %g (nJ)\n", 
//            r->ubuffer.data_array->total_power.readOp.dynamic * 1e9);
    fprintf(stderr, "\tCross bar access energy - %g (nJ)\n", 
                                r->crossbar.power.dynamic * 1e9);
    fprintf(stderr, "\tArbiter access energy - %g (nJ)\n", 
                                r->arbiter.power.dynamic * 1e9);
//    output_UCA(&(r->ubuffer));
}

long double 
read_wordline(router_stats_t *r) {
    long double temp, temp2, temp3;
    temp2 = wcell;
    temp3 = 2;
    temp3*= dw;
    temp3*= (Pr+Pw);
    temp2+= temp3;
    temp2*= r->flit_size;
    temp = Cw3(temp2);
    temp +=  2*r->flit_size*gate_cap(Tpr);
    temp += gate_cap(NTwd) + gate_cap(PTwd) + diff_cap(NTwd, 1, 1) + diff_cap(PTwd, 0 ,1);
    return temp;
}
    
long double 
write_wordline(router_stats_t *r) {
    long double temp;
    temp = Cw3(r->flit_size*(wcell+2*dw*(Pr+Pw))) + 2*r->flit_size*gate_cap(Tpw) + gate_cap(NTwd) + gate_cap(PTwd) + diff_cap(NTwd, 1, 1) + diff_cap(PTwd, 0, 1);
    return temp;
}


long double 
read_bitline(router_stats_t *r) {
    return Cw3(r->buffer_ent*(hcell + dw*(Pr+Pw))) + r->buffer_ent*diff_cap(Tpr, 0, 1) + diff_cap(Tc, 1, 1);
}
    
long double 
write_bitline(router_stats_t *r) {
    long double temp;
    temp = Cw3(r->buffer_ent*(hcell + dw*(Pr+Pw))) + r->buffer_ent*diff_cap(Tpw, 0, 1) +  gate_cap(NTbd) + gate_cap(PTbd) + diff_cap(NTbd, 1, 1) + diff_cap(PTbd, 0, 1);
    return temp;
}

long double 
memory_cell(router_stats_t *r) {
    long double temp;
    temp = 2* (gate_cap(NTm)+gate_cap(PTm)+diff_cap(NTm, 0, 1)+diff_cap(PTm, 1, 1)) +
    2*(Pr*diff_cap(Tpr, 0, 1) + Pw*diff_cap(Tpw, 0, 1));
    return temp;
}

long double 
pre_charge(router_stats_t *r) {
    return gate_cap(Tc);
}

long double 
read_operation(router_stats_t *r) {
    return read_wordline(r)*Vdd*Vdd + read_bitline(r)*r->flit_size*Vdd*Vdd + 2*r->flit_size*pre_charge(r)*Vdd*Vdd;
}

long double 
write_operation(router_stats_t *r) {
    return write_wordline(r)*Vdd*Vdd + r->flit_size*write_bitline(r)*Vdd*Vdd + r->flit_size*memory_cell(r)*Vdd*Vdd*0.5;
}

