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

#include "leakage.h"

double SENSE_AMP_P;
double SENSE_AMP_D;

double
wire_cap2 (int projection,
        int layer,
        double FSIZE)
{
    double sidewall, adj, tot_cap;
    double wire_height;
    double wire_width, wire_spacing;
    double m = 1.2; 
    double eps0 = 8.8541878176e-18; //F/u 
    double epsr = 0.0;
    FSIZE*= 1e-3; //u
    if (projection == 0) { /* aggressive */
        if (layer == 0) { /* local wires */
            wire_width = 1.25*FSIZE;
        }
        else if (layer == 1) { /* semi global */
            wire_width = 2 * FSIZE;
        }
        else if (layer == 2) { /* global */
            wire_width = 4 * FSIZE;
        }
        wire_spacing = wire_width;
        wire_height = wire_width * 3.0; 

        if (FSIZE > .065) {
            epsr = 2.7;
        }
        else if (FSIZE > 0.045) {
            epsr = 2.3;
        }
        else if (FSIZE > 0.032) {
            epsr = 1.958;
        }
        else {
            epsr = 1.664;
        }

    }
    else if (projection == 1) { /* conservative */
        if (layer == 0) { /* local wires */
            wire_width = 1.25*FSIZE;
        }
        else if (layer == 1) { /* semi global */
            wire_width = 2 * FSIZE;
        }
        else if (layer == 2) { /* global */
            wire_width = 4 * FSIZE;
        }
        wire_spacing = wire_width;
        wire_height = wire_width * 2.0; 

        if (FSIZE > .065) {
            epsr = 3.038;
        }
        else if (FSIZE > 0.045) {
            epsr = 2.734;
        }
        else if (FSIZE > 0.032) {
            epsr = 2.460;
        }
        else {
            epsr = 2.214;
        }
    }

    /* capacitance between wires in the same level */
    sidewall = m * eps0 * epsr * (wire_height/wire_spacing);

    /* capacitance between wires in adjacent levels */
    adj = m * eps0 * epsr;

    tot_cap =  (sidewall + adj + (c_fringe[0])); 

    //    fprintf(stderr, "%d %d %g\n", projection, layer, tot_cap);
    return (tot_cap); // (F/u)
}


void init_tech_params(double technology)
{
    int i, iter;
    double tech, tech_low, tech_high, alpha;
    int int_tech_low, int_tech_high;
    FEATURESIZE = technology;
    technology = technology * 1000.0;

    if(technology < 91 && technology > 89){
        tech_low = 90;
        tech_high = 90;
    }
    else if(technology < 69 && technology > 64){
        tech_low = 68;
        tech_high = 68;
    }
    else if(technology < 46 && technology > 44){
        tech_low = 45;
        tech_high = 45;
    }
    else if(technology < 33 && technology > 31){
        tech_low = 32;
        tech_high = 32;
    }
    else if(technology < 90 && technology > 68){
        tech_low = 90;
        tech_high = 68;
    }
    else if(technology < 65 && technology > 45){
        tech_low = 68;
        tech_high = 45;
    }
    else if(technology < 45 && technology > 32){
        tech_low = 45;
        tech_high = 32;
    }

    for(iter=0; iter<=1; ++iter){
        if(iter==0){
            tech = tech_low;
        }
        else{
            tech = tech_high;
        }

        if(tech < 91 && tech > 89){
            SENSE_AMP_D = .28e-9; // s
            SENSE_AMP_P = 14.7e-15; // J
            //For 2005, MPU/ASIC stagger-contacted M1 half-pitch is 90 nm (so this is 90 nm
            //technology i.e. FEATURESIZE = 0.09). For some of the DRAM technology parameters we use 
            //values of 90 nm IBM embedded DRAM technology.
            //90 nm HP
            FO4 = 16.66*3;
            vdd[0] = 1.1;
            Lphy[0] = 0.032;//Lphy is the physical gate-length
            Lelec[0] = 0.0232;//Lelec is the electrical gate-length
            t_ox[0] = 1.2e-3;//micron
            v_th[0] = .194;//V
            c_ox[0] = 1.79e-14;//F/micron2
            mobility_eff[0] = 420.88 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[0] = Lelec[0] * 1.1e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[0] = 5.73e-16;//F/micron
            c_fringe[0] = 0.08e-15;//F/micron
            c_junc[0] = 1e-15;//F/micron2
            I_on_n[0] = 1024e-6;//A/micron
            I_on_p[0] = I_on_n[0] / 2;//A/micron
            Rnchannelon[0] = vdd[0] / I_on_n[0];//ohm-micron
            Rpchannelon[0] = vdd[0] / I_on_p[0];//ohm-micron
            I_off_n[0][0] = 5.97e-8;
            I_off_n[0][10] = 7.25e-8;
            I_off_n[0][20] = 8.69e-8;
            I_off_n[0][30] = 1.03e-7;
            I_off_n[0][40] = 1.21e-7;
            I_off_n[0][50] = 1.42e-7;
            I_off_n[0][60] = 1.63e-7;
            I_off_n[0][70] = 1.86e-7;
            I_off_n[0][80] = 2.08e-7;
            I_off_n[0][90] = 2.27e-7;
            I_off_n[0][100] = 2.54e-7;
            for(i = 0; i <= 100; i += 10){
                I_off_p[0][i] = I_off_n[0][i]; 
            }

            //90 nm LSTP
            vdd[1] = 1.2;
            Lphy[1] = 0.065;
            Lelec[1] = 0.0426;//Lelec is the electrical gate-length.
            t_ox[1] = 2.1e-3;//micron
            v_th[1] = 0.48203;//V
            c_ox[1] = 1.26e-14;//F/micron2
            mobility_eff[1] = 371.83 * (1e-2 * 1e6 * 1e-2 * 1e6);//micron2 / Vs
            Vdsat[1] = Lelec[1] * 1.02e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[1] = 8.21e-16;//F/micron
            c_fringe[1] = 0.08e-15;
            c_junc[1] = 1e-15;//F/micron2
            I_on_n[1] = 497.1e-6;//A/micron
            I_on_p[1] = I_on_n[1] / 2;//239.6e-6;
            Rnchannelon[1] = vdd[1] / I_on_n[1];//ohm-micron
            Rpchannelon[1] = vdd[1] / I_on_p[1];
            I_off_n[1][0] = 1.01e-11;
            I_off_n[1][10] = 1.64e-11;
            I_off_n[1][20] = 2.57e-11;
            I_off_n[1][30] = 3.92e-11;
            I_off_n[1][40] = 5.83e-11;
            I_off_n[1][50] = 8.45e-11;
            I_off_n[1][60] = 1.18e-10;
            I_off_n[1][70] = 1.57e-10;
            I_off_n[1][80] = 1.98-10;
            I_off_n[1][90] = 2.44e-10;
            I_off_n[1][100] = 3.05e-10;
            for(i = 0; i <= 100; i += 10){
                I_off_p[1][i] = I_off_n[1][i]; 
            }

            //90 nm LOP
            vdd[2] = 0.9;
            Lphy[2] = 0.045;
            Lelec[2] = 0.029;//Lelec is the electrical gate-length.	
            t_ox[2] = 1.4e-3;//micron
            v_th[2] = 0.28776;//V
            c_ox[2] = 1.68e-14;//F/micron2
            mobility_eff[2] = 441.9 * (1e-2 * 1e6 * 1e-2 * 1e6);//micron2 / Vs
            Vdsat[2] = Lelec[2] * 0.97e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[2] = 7.57e-16;//F/micron
            c_fringe[2] = 0.08e-15;
            c_junc[2] = 1e-15;//F/micron2
            I_on_n[2] = 588.9e-6;//A/micron
            I_on_p[2] = I_on_n[2] / 2;//239.6e-6;
            Rnchannelon[2] = vdd[2] / I_on_n[2];//ohm-micron
            Rpchannelon[2] = vdd[2] / I_on_p[2];
            I_off_n[2][0] = 2.94e-9;
            I_off_n[2][10] = 3.97e-9;
            I_off_n[2][20] = 5.25e-9;
            I_off_n[2][30] = 6.84e-9;
            I_off_n[2][40] = 8.78e-9;
            I_off_n[2][50] = 1.11e-8;
            I_off_n[2][60] = 1.37e-8;
            I_off_n[2][70] = 1.65e-8;
            I_off_n[2][80] = 1.89e-8;
            I_off_n[2][90] = 2.11e-8;
            I_off_n[2][100] = 2.35e-8;
            for(i = 0; i <= 100; i += 10){
                I_off_p[2][i] = I_off_n[2][i]; 
            }

            //90 nm DRAM cell access transistor technology parameters. Obtained using MASTAR. Based off IBM eDRAM
            //parameters
            vdd_dram_cell = 1.2;
            Lphy[3] = 0.12;//micron
            Lelec[3] = Lphy[3] - 0.3 * Lphy[3];//micron. For now, assume Lelec is 30% lesser than Lphy for DRAM access and wordline transistors.
            v_th_dram_access_transistor = 0.4;//V
            width_dram_access_transistor = 0.14;//micron
            I_on_dram_cell = 62.3e-6;//A
            I_off_dram_cell_worst_case_length_temp = 18.7e-12;//A
            Wmemcella_dram = width_dram_access_transistor;
            Wmemcellpmos_dram = 0;
            Wmemcellnmos_dram = 0;
            area_cell_dram = 0.168;//micron2. //IBM scalable eDRAM cell paper
            asp_ratio_cell_dram = 1.46; //Fixing this to the SRAM cell aspect ratio for now
            c_dram_cell = 20e-15;

            //90 nm DRAM wordline parameters obtained using MASTAR. Based off IBM DRAM parameters. 
            vpp = 1.2 + v_th_dram_access_transistor;//vpp. V
            t_ox[3] = 2.2e-3;//micron
            v_th[3] = 0.4;//V
            c_ox[3] = 1.22e-14;//F/micron2
            mobility_eff[3] =  328.19 * (1e-2 * 1e6 * 1e-2 * 1e6);//micron2 / Vs
            Vdsat[3] = Lelec[3] * 1.17e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[3] = 1.47e-15;//F/micron
            c_fringe[3] = 0.08e-15;//F/micron
            c_junc[3] = 1e-15;//F/micron2
            I_on_n[3] = 103.5e-6 / 0.14;//A/micron
            I_on_p[3] = I_on_n[3] / 2;
            Rnchannelon[3] = vdd[3] / I_on_n[3];//ohm-micron
            Rpchannelon[3] = vdd[3] / I_on_p[3];
            I_off_n[3][0] = 2e-12 / 0.14; //A/micron
            I_off_n[3][10] = 3.14e-12 / 0.14;
            I_off_n[3][20] = 4.81e-12 / 0.14;
            I_off_n[3][30] = 7.81e-12 / 0.14;
            I_off_n[3][40] = 10.5e-12 / 0.14;
            I_off_n[3][50] = 14.9e-12 / 0.14;
            I_off_n[3][60] = 20.5e-12 / 0.14;
            I_off_n[3][70] = 27.1e-12 / 0.14;
            I_off_n[3][80] = 33.9e-12 / 0.14;
            I_off_n[3][90] = 41.3e-12 / 0.14;
            I_off_n[3][100] = 50.7e-12 / 0.14;
            for(i = 0; i <= 100; i += 10){
                I_off_p[3][i] = I_off_n[3][i]; 
            }

            //SRAM cell properties. Based on values specified for standard SRAM cell in IBM's Blue/Gene embedded
            //DRAM technology paper. Area of the cell comes out to be 146F2 where F is the feature size. 
            Wmemcella_sram = 1.31 * FEATURESIZE;//1.23 * FEATURESIZE;
            Wmemcellpmos_sram = 1.23 * FEATURESIZE;//1.69 * FEATURESIZE;
            Wmemcellnmos_sram = 2.08 * FEATURESIZE;//1.15 * FEATURESIZE;
            area_cell_sram = 146 * FEATURESIZE * FEATURESIZE;
            asp_ratio_cell_sram = 1.46;//1.42; 
            Cu_resistivity = 2.53e-8;
        }

        if(tech < 69 && tech > 67){
            SENSE_AMP_D = .2e-9; // s
            SENSE_AMP_P = 5.7e-15; // J
            //For 2007, MPU/ASIC stagger-contacted M1 half-pitch is 68 nm (so this is 68 nm
            //technology i.e. FEATURESIZE = 0.068). 
            //68 nm HP
            FO4 = 26.925;
            vdd[0] = 1.1;
            Lphy[0] = 0.025;//Lphy is the physical gate-length.
            Lelec[0] = 0.019;//Lelec is the electrical gate-length.
            t_ox[0] = 1.1e-3;//micron
            v_th[0] = .16491;//V
            c_ox[0] = 1.88e-14;//F/micron2
            mobility_eff[0] = 420.88 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[0] = Lelec[0] * 1.1e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[0] = 4.69e-16;//F/micron
            c_fringe[0] = 0.077e-15;//F/micron
            c_junc[0] = 1.0e-15;//F/micron2
            I_on_n[0] = 1197.2e-6;//A/micron
            I_on_p[0] = I_on_n[0] / 2;//A/micron
            Rnchannelon[0] = vdd[0] / I_on_n[0];//ohm-micron
            Rpchannelon[0] = vdd[0] / I_on_p[0];//ohm-micron
            I_off_n[0][0] = 1.96e-7;
            I_off_n[0][10] = 2.29e-7;
            I_off_n[0][20] = 2.66e-7;
            I_off_n[0][30] = 3.05e-7;
            I_off_n[0][40] = 3.49e-7;
            I_off_n[0][50] = 3.95e-7;
            I_off_n[0][60] = 4.45e-7;
            I_off_n[0][70] = 4.97e-7;
            I_off_n[0][80] = 5.48e-7;
            I_off_n[0][90] = 5.94e-7;
            I_off_n[0][100] = 6.3e-7;
            for(i = 0; i <= 100; i += 10){
                I_off_p[0][i] = I_off_n[0][i]; 
            }

            //65 nm LSTP
            vdd[1] = 1.2;
            Lphy[1] = 0.045;
            Lelec[1] = 0.0298;//Lelec is the electrical gate-length.
            t_ox[1] = 1.9e-3;//micron
            v_th[1] = 0.52354;//V
            c_ox[1] = 1.36e-14;//F/micron2
            mobility_eff[1] = 341.21 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[1] = Lelec[1] * 1.12e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[1] = 6.14e-16;//F/micron
            c_fringe[1] = 0.08e-15;
            c_junc[1] = 1.0e-15;//F/micron2
            I_on_n[1] = 519.2e-6;//A/micron
            I_on_p[1] = I_on_n[1] / 2;//239.6e-6;
            Rnchannelon[1] = vdd[1] / I_on_n[1];//ohm-micron
            Rpchannelon[1] = vdd[1] / I_on_p[1];
            I_off_n[1][0] = 9.12e-12;
            I_off_n[1][10] = 1.49e-11;
            I_off_n[1][20] = 2.36e-11;
            I_off_n[1][30] = 3.64e-11;
            I_off_n[1][40] = 5.48e-11;
            I_off_n[1][50] = 8.05e-11;
            I_off_n[1][60] = 1.15e-10;
            I_off_n[1][70] = 1.59e-10;
            I_off_n[1][80] = 2.1e-10;
            I_off_n[1][90] = 2.62e-10;
            I_off_n[1][100] = 3.21e-10;
            for(i = 0; i <= 100; i += 10){
                I_off_p[1][i] = I_off_n[1][i]; 
            }

            //65 nm LOP
            vdd[2] = 0.8;
            Lphy[2] = 0.032;
            Lelec[2] = 0.0216;//Lelec is the electrical gate-length.	
            t_ox[2] = 1.2e-3;//micron
            v_th[2] = 0.28512;//V
            c_ox[2] = 1.87e-14;//F/micron2
            mobility_eff[2] = 495.19 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[2] = Lelec[2] * 0.97e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[2] = 6e-16;//F/micron
            c_fringe[2] = 0.08e-15;
            c_junc[2] = 1e-15;//F/micron2. 
            I_on_n[2] = 573.1e-6;//A/micron
            I_on_p[2] = I_on_n[2] / 2;//239.6e-6;
            Rnchannelon[2] = vdd[2] / I_on_n[2];//ohm-micron
            Rpchannelon[2] = vdd[2] / I_on_p[2];
            I_off_n[2][0] = 4.9e-9;
            I_off_n[2][10] = 6.49e-9;
            I_off_n[2][20] = 8.45e-9;
            I_off_n[2][30] = 1.08e-8;
            I_off_n[2][40] = 1.37e-8;
            I_off_n[2][50] = 1.71e-8;
            I_off_n[2][60] = 2.09e-8;
            I_off_n[2][70] = 2.48e-8;
            I_off_n[2][80] = 2.84e-8;
            I_off_n[2][90] = 3.13e-8;
            I_off_n[2][100] = 3.42e-8;
            for(i = 0; i <= 100; i += 10){
                I_off_p[2][i] = I_off_n[2][i]; 
            }

            //65 nm DRAM cell access transistor technology parameters. Obtained using MASTAR. Based off IBM eDRAM
            //parameters
            vdd_dram_cell = 1.2;
            Lphy[3] = 0.12;//micron
            Lelec[3] = Lphy[3] - 0.3 * Lphy[3];//micron. For now, assume Lelec is 30% lesser than Lphy for DRAM access and wordline transistors.
            v_th_dram_access_transistor = 0.38437;//V
            width_dram_access_transistor = 0.09;//micron
            I_on_dram_cell = 41.4e-6;//A
            I_off_dram_cell_worst_case_length_temp = 17.2e-12;//A
            Wmemcella_dram = width_dram_access_transistor;
            Wmemcellpmos_dram = 0;
            Wmemcellnmos_dram = 0;
            area_cell_dram = 18.2*FEATURESIZE*FEATURESIZE;//micron2. //IBM scalable eDRAM cell paper //Temp value FIXME remove this
//            area_cell_dram = 0.11;//micron2. //IBM scalable eDRAM cell paper
            asp_ratio_cell_dram = 1.46; //Fixing this to the SRAM cell aspect ratio
            c_dram_cell = 20e-15;

            //65 nm DRAM wordline parameters obtained using MASTAR. Based off IBM DRAM parameters. 
            vpp = 1.2 + v_th_dram_access_transistor;//vpp. V
            t_ox[3] = 2.2e-3;//micron
            v_th[3] = 0.38437;//V
            c_ox[3] = 1.22e-14;//F/micron2
            mobility_eff[3] =  333.52 * (1e-2 * 1e6 * 1e-2 * 1e6);//micron2 / Vs
            Vdsat[3] = Lelec[3] * 1.15e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[3] = 1.47e-15;//F/micron
            c_fringe[3] = 0.08e-15;//F/micron
//            c_junc[3] = 2.34e-17 / 0.09 ;//F/micron. Calculated for a width of 0.14 micron. 
            c_junc[3] = 1e-15;//F/micron2. 
            I_on_n[3] = 67e-6 / 0.09;//A/micron
            I_on_p[3] = I_on_n[3] / 2;
            Rnchannelon[3] = vdd[3] / I_on_n[3];//ohm-micron
            Rpchannelon[3] = vdd[3] / I_on_p[3];
            I_off_n[3][0] = 2e-12 / 0.09; //A/micron
            I_off_n[3][10] = 3.1e-12 / 0.09;
            I_off_n[3][20] = 4.67e-12 / 0.09;
            I_off_n[3][30] = 6.87e-12 / 0.09;
            I_off_n[3][40] = 9.89e-12 / 0.09;
            I_off_n[3][50] = 13.9e-12 / 0.09;
            I_off_n[3][60] = 18.9e-12 / 0.09;
            I_off_n[3][70] = 24.7e-12 / 0.09;
            I_off_n[3][80] = 30.5e-12 / 0.09;
            I_off_n[3][90] = 36.9e-12 / 0.09;
            I_off_n[3][100] = 45.1e-12 / 0.09;
            for(i = 0; i <= 100; i += 10){
                I_off_p[3][i] = I_off_n[3][i]; 
            }

            //SRAM cell properties. Based on values specified for standard SRAM cell in IBM's Blue/Gene embedded
            //DRAM technology paper. Area of the cell comes out to be 146F2 where F is the feature size. 
            Wmemcella_sram = 1.31 * FEATURESIZE;//1.23 * FEATURESIZE;
            Wmemcellpmos_sram = 1.23 * FEATURESIZE;//1.69 * FEATURESIZE;
            Wmemcellnmos_sram = 2.08 * FEATURESIZE;//1.15 * FEATURESIZE;
            area_cell_sram = 146 * FEATURESIZE * FEATURESIZE;
            asp_ratio_cell_sram = 1.46;//1.42; 
            Cu_resistivity = 2.73e-8;
        }

        if(tech < 46 && tech > 44){
            SENSE_AMP_D = .04e-9; // s
            SENSE_AMP_P = 2.7e-15; // J
            //For 2010, MPU/ASIC stagger-contacted M1 half-pitch is 45 nm (so this is 
            //45 nm technology i.e. FEATURESIZE = 0.045). For some of the DRAM technology 
            //parameters we use values of 45 nm IBM embedded DRAM technology.
            //45 nm HP
            FO4 = 16.575;
            vdd[0] = 1;
            Lphy[0] = 0.018;//Lphy is the physical gate-length.
            Lelec[0] = 0.01345;//Lelec is the electrical gate-length.
            t_ox[0] = 0.65e-3;//micron
            v_th[0] = .15091;//V
            c_ox[0] = 3.77e-14;//F/micron2
            mobility_eff[0] = 266.68 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[0] = Lelec[0] * 1.1e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[0] = 6.78e-16;//F/micron
            c_fringe[0] = 0.05e-15;//F/micron
            c_junc[0] = 1.0e-15;//F/micron2
            I_on_n[0] = 2046.6e-6;//A/micron
            I_on_p[0] = I_on_n[0] / 2;//A/micron
            Rnchannelon[0] = vdd[0] / I_on_n[0];//ohm-micron
            Rpchannelon[0] = vdd[0] / I_on_p[0];//ohm-micron
            I_off_n[0][0] = 2.8e-7;
            I_off_n[0][10] = 3.28e-7;
            I_off_n[0][20] = 3.81e-7;
            I_off_n[0][30] = 4.39e-7;
            I_off_n[0][40] = 5.02e-7;
            I_off_n[0][50] = 5.69e-7;
            I_off_n[0][60] = 6.42e-7;
            I_off_n[0][70] = 7.2e-7;
            I_off_n[0][80] = 8.03e-7;
            I_off_n[0][90] = 8.91e-7;
            I_off_n[0][100] = 9.84e-7;
            for(i = 0; i <= 100; i += 10){
                I_off_p[0][i] = I_off_n[0][i]; 
            }

            //45 nm LSTP
            vdd[1] = 1.1;
            Lphy[1] =  0.028;
            Lelec[1] = 0.0212;//Lelec is the electrical gate-length.
            t_ox[1] = 1.4e-3;//micron
            v_th[1] = 0.50245;//V
            c_ox[1] = 2.01e-14;//F/micron2
            mobility_eff[1] =  363.96 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[1] = Lelec[1] * 1.35e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[1] = 5.64e-16;//F/micron
            c_fringe[1] = 0.08e-15;
            c_junc[1] = 1.0e-15;//F/micron2
            I_on_n[1] = 666.2e-6;//A/micron
            I_on_p[1] = I_on_n[1] / 2;//239.6e-6;
            Rnchannelon[1] = vdd[1] / I_on_n[1];//ohm-micron
            Rpchannelon[1] = vdd[1] / I_on_p[1];
            I_off_n[1][0] = 1.01e-11;
            I_off_n[1][10] = 1.65e-11;
            I_off_n[1][20] = 2.62e-11;
            I_off_n[1][30] = 4.06e-11;
            I_off_n[1][40] = 6.12e-11;
            I_off_n[1][50] = 9.02e-11;
            I_off_n[1][60] = 1.3e-10;
            I_off_n[1][70] = 1.83e-10;
            I_off_n[1][80] = 2.51e-10;
            I_off_n[1][90] = 3.29e-10;
            I_off_n[1][100] = 4.1e-10;
            for(i = 0; i <= 100; i += 10){
                I_off_p[1][i] = I_off_n[1][i]; 
            }

            //45 nm LOP
            vdd[2] = 0.7;
            Lphy[2] = 0.022;
            Lelec[2] = 0.016;//Lelec is the electrical gate-length.	
            t_ox[2] = 0.9e-3;//micron
            v_th[2] = 0.22599;//V
            c_ox[2] = 2.82e-14;//F/micron2
            mobility_eff[2] = 508.9 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[2] = Lelec[2] * 0.93e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[2] = 6.2e-16;//F/micron
            c_fringe[2] = 0.073e-15;
            c_junc[2] = 1.0e-15;//F/micron2
            I_on_n[2] = 748.9e-6;//A/micron			
            I_on_p[2] = I_on_n[2] / 2;//239.6e-6;
            Rnchannelon[2] = vdd[2] / I_on_n[2];//ohm-micron
            Rpchannelon[2] = vdd[2] / I_on_p[2];
            I_off_n[2][0] = 4.03e-9;
            I_off_n[2][10] = 5.02e-9;
            I_off_n[2][20] = 6.18e-9;
            I_off_n[2][30] = 7.51e-9;
            I_off_n[2][40] = 9.04e-9;
            I_off_n[2][50] = 1.08e-8;
            I_off_n[2][60] = 1.27e-8;
            I_off_n[2][70] = 1.47e-8;
            I_off_n[2][80] = 1.66e-8;
            I_off_n[2][90] = 1.84e-8;
            I_off_n[2][100] = 2.03e-8;
            for(i = 0; i <= 100; i += 10){
                I_off_p[2][i] = I_off_n[2][i]; 
            }

            //45 nm DRAM cell access transistor technology parameters. Obtained using MASTAR. Based off scaling 
            //IBM eDRAM parameters
            vdd_dram_cell = 1.1;
            Lphy[3] = 0.078;//micron
            Lelec[3] = Lphy[3] - 0.3 * Lphy[3];//micron. For now, assume Lelec is 30% lesser than Lphy for DRAM access and wordline transistors.
            v_th_dram_access_transistor = 0.41743;//V
            width_dram_access_transistor = 0.079;//micron
            I_on_dram_cell = 35e-6;//A
            I_off_dram_cell_worst_case_length_temp = 18.2e-12;//A
            Wmemcella_dram = width_dram_access_transistor;
            Wmemcellpmos_dram = 0;
            Wmemcellnmos_dram  = 0;
            area_cell_dram = width_dram_access_transistor * Lphy[3] * 10.0;//micron2. //Based on fitted area model for
            //IBM eDRAM cell
            asp_ratio_cell_dram = 1.46; //Fixing this to the SRAM cell aspect ratio
            c_dram_cell = 20e-15;

            //45 nm DRAM wordline parameters obtained using MASTAR. 
            vpp = 1.1 + v_th_dram_access_transistor;//vpp. V
            t_ox[3] = 2.1e-3;//micron
            v_th[3] = 0.41743;//V
            c_ox[3] = 1.27e-14;//F/micron2
            mobility_eff[3] =   326.19 * (1e-2 * 1e6 * 1e-2 * 1e6);//micron2 / Vs
            Vdsat[3] = Lelec[3] * 1.18e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[3] = 9.91e-16;//F/micron
            c_fringe[3] = 0.08e-15;//F/micron
            c_junc[3] = 1.0e-15;//F/micron2
            I_on_n[3] = 62.6e-6 / 0.078;//A/micron
            I_on_p[3] = I_on_n[3] / 2;
            Rnchannelon[3] = vdd[3] / I_on_n[3];//ohm-micron
            Rpchannelon[3] = vdd[3] / I_on_p[3];
            I_off_n[3][0] = 2e-12 / 0.078; //A/micron
            I_off_n[3][10] = 3.11e-12 / 0.078;
            I_off_n[3][20] = 4.72e-12 / 0.078;
            I_off_n[3][30] = 7e-12 / 0.078;
            I_off_n[3][40] = 10.2e-12 / 0.078;
            I_off_n[3][50] = 14.4e-12 / 0.078;
            I_off_n[3][60] = 19.9e-12 / 0.078;
            I_off_n[3][70] = 26.4e-12 / 0.078;
            I_off_n[3][80] = 33.5e-12 / 0.078;
            I_off_n[3][90] = 40.7e-12 / 0.078;
            I_off_n[3][100] = 49.3e-12 / 0.078;
            for(i = 0; i <= 100; i += 10){
                I_off_p[3][i] = I_off_n[3][i]; 
            }


            //SRAM cell properties. Based on values specified for standard SRAM cell in IBM's Blue/Gene embedded
            //DRAM technology paper. Area of the cell comes out to be 146F2 where F is the feature size. 
            Wmemcella_sram = 1.31 * FEATURESIZE;//1.23 * FEATURESIZE;
            Wmemcellpmos_sram = 1.23 * FEATURESIZE;//1.69 * FEATURESIZE;
            Wmemcellnmos_sram = 2.08 * FEATURESIZE;//1.15 * FEATURESIZE;
            area_cell_sram = 146 * FEATURESIZE * FEATURESIZE;
            asp_ratio_cell_sram = 1.46;//1.42; 
            Cu_resistivity = 3.1e-8;
        }

        if(tech < 33 && tech > 31){
            SENSE_AMP_D = .03e-9; // s
            SENSE_AMP_P = 2.16e-15; // J
            //For 2013, MPU/ASIC stagger-contacted M1 half-pitch is 32 nm (so this is 32 nm
            //technology i.e. FEATURESIZE = 0.032). Using the SOI process numbers for 
            //HP and LSTP.
            //32 nm HP
            FO4 = 10.879;
            vdd[0] = 0.9;
            Lphy[0] = 0.013;//Lphy is the physical gate-length.
            Lelec[0] = 0.01013;//Lelec is the electrical gate-length.
            t_ox[0] = 0.5e-3;//micron
            v_th[0] = 0.18835;//V
            c_ox[0] = 4.11e-14;//F/micron2
            mobility_eff[0] = 361.84 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[0] = Lelec[0] * 1.1e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[0] = 5.34e-16;//F/micron
            c_fringe[0] = 0.04e-15;//F/micron
            c_junc[0] = 1.0e-15;//F/micron2
            I_on_n[0] =  2211.7e-6;//A/micron
            I_on_p[0] = I_on_n[0] / 2;//A/micron
            Rnchannelon[0] = vdd[0] / I_on_n[0];//ohm-micron
            Rpchannelon[0] = vdd[0] / I_on_p[0];//ohm-micron
            I_off_n[0][0] = 1.52e-7;
            I_off_n[0][10] = 1.55e-7;
            I_off_n[0][20] = 1.59e-7;
            I_off_n[0][30] = 1.68e-7;
            I_off_n[0][40] = 1.9e-7;
            I_off_n[0][50] = 2.69e-7;
            I_off_n[0][60] = 5.32e-7;
            I_off_n[0][70] = 1.02e-6;
            I_off_n[0][80] = 1.62e-6;
            I_off_n[0][90] = 2.73e-6;
            I_off_n[0][100] = 6.12e-6;
            for(i = 0; i <= 100; i += 10){
                I_off_p[0][i] = I_off_n[0][i]; 
            }

            //32 nm LSTP
            vdd[1] = 1;
            Lphy[1] = 0.020;
            Lelec[1] = 0.01689;//Lelec is the electrical gate-length.
            t_ox[1] = 1.1e-3;//micron
            v_th[1] = 0.47168;//V
            c_ox[1] = 2.29e-14;//F/micron2
            mobility_eff[1] =  619.24 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[1] = Lelec[1] * 0.51e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[1] = 4.28e-16;//F/micron
            c_fringe[1] = 0.063e-15;
            c_junc[1] = 1.0e-15;//F/micron2
            I_on_n[1] = 694e-6;//A/micron
            I_on_p[1] = I_on_n[1] / 2;//239.6e-6;
            Rnchannelon[1] = vdd[1] / I_on_n[1];//ohm-micron
            Rpchannelon[1] = vdd[1] / I_on_p[1];
            I_off_n[1][0] = 1.52e-11;
            I_off_n[1][10] = 1.62e-11;
            I_off_n[1][20] = 1.85e-11;
            I_off_n[1][30] = 2.38e-11;
            I_off_n[1][40] = 3.16e-11;
            I_off_n[1][50] = 4.15e-11;
            I_off_n[1][60] = 5.98e-11;
            I_off_n[1][70] = 9.62e-11;
            I_off_n[1][80] = 1.62e-10;
            I_off_n[1][90] = 2.66e-10;
            I_off_n[1][100] = 3.99e-10;
            for(i = 0; i <= 100; i += 10){
                I_off_p[1][i] = I_off_n[1][i]; 
            }

            //32 nm LOP
            vdd[2] = 0.6;
            Lphy[2] = 0.016;
            Lelec[2] = 0.01172;//Lelec is the electrical gate-length.	
            t_ox[2] = 0.8e-3;//micron
            v_th[2] = 0.0521;//V
            c_ox[2] = 1.69e-14;//F/micron2
            mobility_eff[2] =   751.71 * (1e-2 * 1e6 * 1e-2 * 1e6); //micron2 / Vs
            Vdsat[2] = Lelec[2] * 0.42e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[2] = 2.7e-16;//F/micron
            c_fringe[2] = 0.06e-15;
            c_junc[2] = 1.0e-15;//F/micron2
            I_on_n[2] = 843.4e-6;//A/micron
            I_on_p[2] = I_on_n[2] / 2;
            Rnchannelon[2] = vdd[2] / I_on_n[2];//ohm-micron
            Rpchannelon[2] = vdd[2] / I_on_p[2];
            I_off_n[2][0] = 8.41e-6;
            I_off_n[2][10] = 9.52e-5;
            I_off_n[2][20] = 9.52e-5;//MASTAR does not generate numbers for > 320 deg, so simply fixing
            //leakage currents for temp > 320 equal to the 320 deg value. 
            I_off_n[2][30] = 9.52e-5;
            I_off_n[2][40] = 9.52e-5;
            I_off_n[2][50] = 9.52e-5;
            I_off_n[2][60] = 9.52e-5;
            I_off_n[2][70] = 9.52e-5;
            I_off_n[2][80] = 9.52e-5;
            I_off_n[2][90] = 9.52e-5;
            I_off_n[2][100] = 9.52e-5;
            for(i = 0; i <= 100; i += 10){
                I_off_p[2][i] = I_off_n[2][i]; 
            }

            //32 nm DRAM cell access transistor technology parameters. Obtained using MASTAR. Based off scaling 
            //IBM eDRAM parameters. Using an SOI process.
            vdd_dram_cell = 1.1;
            Lphy[3] = 0.056;//micron
            Lelec[3] = Lphy[3] - 0.3 * Lphy[3];//micron. For now, assume Lelec is 30% lesser than Lphy for DRAM access and wordline transistors.
            v_th_dram_access_transistor = 0.34001;//V
            width_dram_access_transistor = 0.056;//micron
            I_on_dram_cell = 34.4e-6;//A
            I_off_dram_cell_worst_case_length_temp = 18.8e-12;//A
            Wmemcella_dram = width_dram_access_transistor;
            Wmemcellpmos_dram = 0;
            Wmemcellnmos_dram = 0;
            area_cell_dram = width_dram_access_transistor * Lphy[3] * 10.0;//micron2. //Based on fitted area model for
            //IBM eDRAM cell
            asp_ratio_cell_dram = 1.46; //Fixing this to the SRAM cell aspect ratio
            c_dram_cell = 20e-15;

            //32 nm DRAM wordline parameters obtained using MASTAR. 
            vpp = vdd_dram_cell + v_th_dram_access_transistor;//vpp. V
            t_ox[3] = 2e-3;//micron
            v_th[3] = 0.34001;//V
            c_ox[3] = 1.27e-14;//F/micron2
            mobility_eff[3] =  338.73 * (1e-2 * 1e6 * 1e-2 * 1e6);//micron2 / Vs
            Vdsat[3] = Lelec[3] * 1.13e+6 / (1e-2 * 1e6); //V/micron
            c_g_ideal[3] = 7.45e-16;//F/micron
            c_fringe[3] = 0.08e-15;//F/micron
            c_junc[3] = 0;//F/micron. 
            I_on_n[3] = 52.6e-6 / 0.056;//A/micron
            I_on_p[3] = I_on_n[3] / 2;
            Rnchannelon[3] = vdd[3] / I_on_n[3];//ohm-micron
            Rpchannelon[3] = vdd[3] / I_on_p[3];
            I_off_n[3][0] = 2e-12 / 0.056; //A/micron
            I_off_n[3][10] = 3.12e-12 / 0.056;
            I_off_n[3][20] = 4.73e-12 / 0.056;
            I_off_n[3][30] = 7.01e-12 / 0.056;
            I_off_n[3][40] = 10.2e-12 / 0.056;
            I_off_n[3][50] = 14.4e-12 / 0.056;
            I_off_n[3][60] = 20.1e-12 / 0.056;
            I_off_n[3][70] = 27.3e-12 / 0.056;
            I_off_n[3][80] = 35.9e-12 / 0.056;
            I_off_n[3][90] = 45.1e-12 / 0.056;
            I_off_n[3][100] = 53.7e-12 / 0.056;
            for(i = 0; i <= 100; i += 10){
                I_off_p[3][i] = I_off_n[3][i]; 
            }

            //SRAM cell properties. Based on values specified for standard SRAM cell in IBM's Blue/Gene embedded
            //DRAM technology paper. Area of the cell comes out to be 146F2 where F is the feature size. 
            Wmemcella_sram = 1.31 * FEATURESIZE;//1.23 * FEATURESIZE;
            Wmemcellpmos_sram = 1.23 * FEATURESIZE;//1.69 * FEATURESIZE;
            Wmemcellnmos_sram = 2.08 * FEATURESIZE;//1.15 * FEATURESIZE;
            area_cell_sram = 146 * FEATURESIZE * FEATURESIZE;
            asp_ratio_cell_sram = 1.46;//1.42;  
            Cu_resistivity = 3.52e-8;
        }

        dram_cell_tech_flavor = 3;

        //Properties of peripheral/global circuitry transistors
        vdd_periph_global_tech_node[iter] = vdd[periph_global_tech_flavor];
        t_ox_periph_global_tech_node[iter] = t_ox[periph_global_tech_flavor];//micron
        v_th_periph_global_tech_node[iter] = v_th[periph_global_tech_flavor];//V
        c_ox_periph_global_tech_node[iter] = c_ox[periph_global_tech_flavor];//F/micron2
        c_g_ideal_itrs_periph_global_tech_node[iter] = c_g_ideal[periph_global_tech_flavor];//F/micron
        c_fringe_itrs_periph_global_tech_node[iter] = c_fringe[periph_global_tech_flavor];//F/micron
        c_junc_itrs_periph_global_tech_node[iter] = c_junc[periph_global_tech_flavor];//F/micron
        Lphy_periph_global_tech_node[iter] = Lphy[periph_global_tech_flavor];//Lphy is the physical gate-length.
        Lelec_periph_global_tech_node[iter] = Lelec[periph_global_tech_flavor];//Lphy is the electrical gate-length.
        I_on_n_periph_global_tech_node[iter] = I_on_n[periph_global_tech_flavor];//A/micron
        for(i = 0; i <= 100; i += 10){
            I_off_n_periph_global_tech_node[iter][i] = I_off_n[periph_global_tech_flavor][i];
            I_off_p_periph_global_tech_node[iter][i] = I_off_n_periph_global[i]; 
        }

        //Properties of transistors within SRAM cell. We assume that these transistors follow the LSTP
        //roadmap i.e. their characteristics are same as those of the LSTP transistors of 
        //the same ITRS technology node/year.
        vdd_sram_cell_tech_node[iter] = vdd[sram_cell_and_wordline_tech_flavor];
        Lphy_sram_cell_transistor_tech_node[iter] = Lphy[sram_cell_and_wordline_tech_flavor];
        Lelec_sram_cell_transistor_tech_node[iter] = Lelec[sram_cell_and_wordline_tech_flavor];
        t_ox_sram_cell_transistor_tech_node[iter] = t_ox[sram_cell_and_wordline_tech_flavor];//micron
        v_th_sram_cell_transistor_tech_node[iter] = v_th[sram_cell_and_wordline_tech_flavor];//V
        c_g_ideal_itrs_sram_cell_transistor_tech_node[iter] = c_g_ideal[sram_cell_and_wordline_tech_flavor];//F/micron
        c_fringe_itrs_sram_cell_transistor_tech_node[iter] = c_fringe[sram_cell_and_wordline_tech_flavor];
        c_junc_itrs_sram_cell_transistor_tech_node[iter] = c_junc[sram_cell_and_wordline_tech_flavor];//F/micron
        I_on_n_sram_cell_transistor_tech_node[iter] = I_on_n[sram_cell_and_wordline_tech_flavor];//A/micron
        for(i = 0; i <= 100; i += 10){
            I_off_n_sram_cell_transistor_tech_node[iter][i] = I_off_n[sram_cell_and_wordline_tech_flavor][i];
            I_off_p_sram_cell_transistor_tech_node[iter][i] = I_off_n_sram_cell_transistor[i]; 
        }

        //Properties of transistors used within the DRAM cells and those used within the DRAM wordline 
        //driver circuitry. We assume that these transistors follow the LSTP roadmap i.e. their 
        //characteristics are those of the LSTP transistors of the same ITRS technology node/year. 
        //For the DRAM access transistor, the length and width of the transistor is based on IBM's 65 nm 
        //scalable DRAM cell.

        vdd_dram_cell_tech_node[iter] = vdd_dram_cell;
        v_th_dram_access_transistor_tech_node[iter] = v_th_dram_access_transistor;
        Lphy_dram_access_transistor_tech_node[iter] = Lphy[dram_cell_tech_flavor];//micron
        Lelec_dram_access_transistor_tech_node[iter] = Lelec[dram_cell_tech_flavor];//micron
        c_g_ideal_itrs_dram_access_transistor_tech_node[iter] = c_g_ideal[dram_cell_tech_flavor];
        c_fringe_itrs_dram_access_transistor_tech_node[iter]  = c_fringe[dram_cell_tech_flavor];
        c_junc_itrs_dram_access_transistor_tech_node[iter] = c_junc[dram_cell_tech_flavor];
        I_on_n_dram_access_transistor_tech_node[iter] = I_on_n[dram_cell_tech_flavor];
        c_dram_cell_tech_node[iter] = c_dram_cell;

        vpp_tech_node[iter] = vpp;
        Lphy_dram_wordline_transistor_tech_node[iter] = Lphy[dram_cell_tech_flavor];
        Lelec_dram_wordline_transistor_tech_node[iter] = Lelec[dram_cell_tech_flavor];
        c_g_ideal_itrs_dram_wordline_transistor_tech_node[iter] = c_g_ideal[dram_cell_tech_flavor];
        c_fringe_itrs_dram_wordline_transistor_tech_node[iter]  = c_fringe[dram_cell_tech_flavor];
        c_junc_itrs_dram_wordline_transistor_tech_node[iter] = c_junc[dram_cell_tech_flavor];
        I_on_n_dram_wordline_transistor_tech_node[iter] = I_on_n[dram_cell_tech_flavor];
        for(i = 0; i <= 100; i += 10){
            I_off_n_dram_wordline_transistor_tech_node[iter][i] = I_off_n[dram_cell_tech_flavor][i]; 
            I_off_p_dram_wordline_transistor_tech_node[iter][i] = I_off_p[dram_cell_tech_flavor][i];
        }

        Wmemcella_dram_tech_node[iter] = Wmemcella_dram;
        Wmemcellpmos_dram_tech_node[iter] = Wmemcellpmos_dram;
        Wmemcellnmos_dram_tech_node[iter] = Wmemcellnmos_dram;
        area_cell_dram_tech_node[iter] = area_cell_dram;
        asp_ratio_cell_dram_tech_node[iter] = asp_ratio_cell_dram;
        Wmemcella_sram_tech_node[iter] = Wmemcella_sram;
        Wmemcellpmos_sram_tech_node[iter] = Wmemcellpmos_sram;
        Wmemcellnmos_sram_tech_node[iter] = Wmemcellnmos_sram;
        area_cell_sram_tech_node[iter] = area_cell_sram;
        asp_ratio_cell_sram_tech_node[iter] = asp_ratio_cell_sram;

        //Sense amplifier latch Gm calculation
        mobility_eff_periph_global_tech_node[iter] = mobility_eff[periph_global_tech_flavor]; //micron2 / Vs
        Vdsat_periph_global_tech_node[iter] = Vdsat[periph_global_tech_flavor]; //V/micron
    }

    //Currently we are not modeling the resistance/capacitance of poly anywhere. But this can 
    //be added in the future.

    Cpolywire = 0;
    Wcompinvp1 = 12.5 * FEATURESIZE;//this was 10 micron for the 0.8 micron process
    Wcompinvn1 = 7.5 * FEATURESIZE;//this was 6 micron for the 0.8 micron process
    Wcompinvp2 = 25 * FEATURESIZE;//this was 20 micron for the 0.8 micron process
    Wcompinvn2 = 15 * FEATURESIZE;//this was 12 micron for the 0.8 micron process
    Wcompinvp3 = 50 * FEATURESIZE;//this was 40 micron for the 0.8 micron process
    Wcompinvn3 = 30 * FEATURESIZE;//this was 24 micron for the 0.8 micron process
    Wevalinvp =  100 * FEATURESIZE;//this was 80 micron for the 0.8 micron process
    Wevalinvn = 50 * FEATURESIZE;//this was 40 micron for the 0.8 micron process
    Wcompn =  12.5 * FEATURESIZE;//this was 10 micron for the 0.8 micron process
    Wcompp =  37.5 * FEATURESIZE;//this was 30 micron for the 0.8 micron process
    //WmuxdrvNANDn and WmuxdrvNANDp are no longer being used but it's part of the old
    //delay_comparator function which we are using exactly as it used to be, so just setting these
    //to 0
    WmuxdrvNANDn = 0;
    WmuxdrvNANDp = 0;

    MIN_GAP_BET_P_AND_N_DIFFS = 5 * FEATURESIZE;//This is kind of a reasonable but arbitrary number 
    //currently. CHECK WITH NORM.
    MIN_GAP_BET_SAME_TYPE_DIFFS = 1.5 * FEATURESIZE;
    HPOWERRAIL = 2 * FEATURESIZE;
    DEFAULTHEIGHTCELL = 50 * FEATURESIZE; //This is kind of a reasonable but arbitrary number 
    //currently. CHECK WITH NORM.
    WIDTHPOLYCONTACT = FEATURESIZE;
    SPACINGPOLYTOCONTACT = FEATURESIZE;
    SPACINGPOLYTOPOLY = 1.5 * FEATURESIZE;
    MAX_NMOS_WIDTH = 100 * FEATURESIZE;
    MAX_PMOS_WIDTH = 200 * FEATURESIZE;
    gnand2 = 4.0 / 3.0;
    gnand3 = 6.0 / 3.0;
    gnor2 = 5.0 / 3.0;
    gpmos = 2.0 / 3.0;
    gnmos = 1.0 / 3.0;
    minimum_width_nmos = 3 * FEATURESIZE / 2;
    minimum_width_pmos = 2 * minimum_width_nmos;
    Wiso = 12.5*FEATURESIZE;//was 10 micron for the 0.8 micron process
    WsenseN = 3.75*FEATURESIZE; // sense amplifier N-trans; was 3 micron for the 0.8 micron process
    WsenseP = 7.5*FEATURESIZE; // sense amplifier P-trans; was 6 micron for the 0.8 micron process
    WsenseEn = 5*FEATURESIZE; // Sense enable transistor of the sense amplifier; was 4 micron for the 0.8 micron process
    width_nmos_bit_mux = 6 * minimum_width_nmos;
    width_nmos_sense_amp_mux = 6 * minimum_width_nmos;
    width_pmos_bitline_precharge = 2 * minimum_width_pmos;
    width_pmos_bitline_equalization = minimum_width_pmos;

    int_tech_low = (int) (floor(tech_low + 0.5));
    int_tech_high = (int) (floor(tech_high+ 0.5));
    if(int_tech_low != int_tech_high){
        alpha = (technology - tech_low) / (tech_high - tech_low);
    }
    else{
        alpha =0;
    }

    vdd_periph_global = vdd_periph_global_tech_node[0] + alpha * (vdd_periph_global_tech_node[1] - 
            vdd_periph_global_tech_node[0]);
    Lphy_periph_global = Lphy_periph_global_tech_node[0] + 
        alpha * (Lphy_periph_global_tech_node[1] - Lphy_periph_global_tech_node[0]);
    Lelec_periph_global = Lelec_periph_global_tech_node[0] + 
        alpha * (Lelec_periph_global_tech_node[1] - Lelec_periph_global_tech_node[0]);
    t_ox_periph_global = t_ox_periph_global_tech_node[0] + 
        alpha * (t_ox_periph_global_tech_node[1] - t_ox_periph_global_tech_node[0]);//micron
    v_th_periph_global = v_th_periph_global_tech_node[0] + 
        alpha * (v_th_periph_global_tech_node[1] - v_th_periph_global_tech_node[0]);//V
    c_ox_periph_global = c_ox_periph_global_tech_node[0] + 
        alpha * (c_ox_periph_global_tech_node[1] - c_ox_periph_global_tech_node[0]);//F/micron2
    c_g_ideal_itrs_periph_global = c_g_ideal_itrs_periph_global_tech_node[0] + 
        alpha * (c_g_ideal_itrs_periph_global_tech_node[1] - c_g_ideal_itrs_periph_global_tech_node[0]);//F/micron
    c_fringe_itrs_periph_global = c_fringe_itrs_periph_global_tech_node[0] + 
        alpha * (c_fringe_itrs_periph_global_tech_node[1] - c_fringe_itrs_periph_global_tech_node[0]);//F/micron
    c_junc_itrs_periph_global = c_junc_itrs_periph_global_tech_node[0] + 
        alpha * (c_junc_itrs_periph_global_tech_node[1] - c_junc_itrs_periph_global_tech_node[0]);//F/micron
    c_overlap_itrs_periph_global = (((Lphy_periph_global - Lelec_periph_global) / Lphy_periph_global) / 2) * c_g_ideal_itrs_periph_global;//F/micron

    I_on_n_periph_global = I_on_n_periph_global_tech_node[0] + 
        alpha * (I_on_n_periph_global_tech_node[1] - I_on_n_periph_global_tech_node[0]);//A/micron
    I_on_p_periph_global = I_on_n_periph_global / 2;//A/micron
    Rnchannelon_itrs_periph_global = vdd_periph_global / I_on_n_periph_global;//ohm-micron
    Rpchannelon_itrs_periph_global = vdd_periph_global / I_on_p_periph_global;//ohm-micron
    for(i = 0; i <= 100; i += 10){
        I_off_n_periph_global[i] = I_off_n_periph_global_tech_node[0][i] + 
            alpha * (I_off_n_periph_global_tech_node[1][i] - I_off_n_periph_global_tech_node[0][i] );//A/micron
        I_off_p_periph_global[i] = I_off_n_periph_global[i]; 
    }

    //Properties of transistors within SRAM cell. We assume that these transistors follow the LSTP
    //roadmap i.e. their characteristics are same as those of the LSTP transistors of 
    //the same ITRS technology node/year.
    vdd_sram_cell = vdd_sram_cell_tech_node[0] + alpha * (vdd_sram_cell_tech_node[1] - vdd_sram_cell_tech_node[0]);
    Lphy_sram_cell_transistor = Lphy_sram_cell_transistor_tech_node[0] + 
        alpha * (Lphy_sram_cell_transistor_tech_node[1] - Lphy_sram_cell_transistor_tech_node[0]);//micron. 
    Lelec_sram_cell_transistor = Lelec_sram_cell_transistor_tech_node[0] + 
        alpha * (Lelec_sram_cell_transistor_tech_node[1] - Lelec_sram_cell_transistor_tech_node[0]);//micron.
    t_ox_sram_cell_transistor = t_ox_sram_cell_transistor_tech_node[0] + 
        alpha * (t_ox_sram_cell_transistor_tech_node[1] - t_ox_sram_cell_transistor_tech_node[0]);//micron
    v_th_sram_cell_transistor = v_th_sram_cell_transistor_tech_node[0] + 
        alpha * (v_th_sram_cell_transistor_tech_node[1] - v_th_sram_cell_transistor_tech_node[0]);//V
    c_g_ideal_itrs_sram_cell_transistor = c_g_ideal_itrs_sram_cell_transistor_tech_node[0] + 
        alpha * (c_g_ideal_itrs_sram_cell_transistor_tech_node[1] - c_g_ideal_itrs_sram_cell_transistor_tech_node[0]);//F/micron
    c_fringe_itrs_sram_cell_transistor = c_fringe_itrs_sram_cell_transistor_tech_node[0] + 
        alpha * (c_fringe_itrs_sram_cell_transistor_tech_node[1] - c_fringe_itrs_sram_cell_transistor_tech_node[0]);//F/micron
    c_junc_itrs_sram_cell_transistor = c_junc_itrs_sram_cell_transistor_tech_node[0] + 
        alpha * (c_junc_itrs_sram_cell_transistor_tech_node[1] - c_junc_itrs_sram_cell_transistor_tech_node[0]);//F/micron
    c_overlap_itrs_sram_cell_transistor = (((Lphy_sram_cell_transistor - Lelec_sram_cell_transistor) / Lphy_sram_cell_transistor) / 2) * c_g_ideal_itrs_sram_cell_transistor;//F/micron
    I_on_n_sram_cell_transistor = I_on_n_sram_cell_transistor_tech_node[0] + 
        alpha * (I_on_n_sram_cell_transistor_tech_node[1] - I_on_n_sram_cell_transistor_tech_node[0]);//A/micron
    I_on_p_sram_cell_transistor = I_on_n_sram_cell_transistor / 2;//A/micron
    Rnchannelon_itrs_sram_cell_transistor = vdd_sram_cell / I_on_n_sram_cell_transistor;//ohm-micron
    Rpchannelon_itrs_sram_cell_transistor = vdd_sram_cell / I_on_p_sram_cell_transistor;//ohm-micron
    for(i = 0; i <= 100; i += 10){
        I_off_n_sram_cell_transistor[i] = I_off_n_sram_cell_transistor_tech_node[0][i] + 
            alpha * (I_off_n_sram_cell_transistor_tech_node[1][i] - I_off_n_sram_cell_transistor_tech_node[0][i] );//A/micron
        I_off_p_sram_cell_transistor[i] = I_off_n_sram_cell_transistor[i]; 
    }
    Wmemcella_sram = Wmemcella_sram_tech_node[0] + alpha * (Wmemcella_sram_tech_node[1] - Wmemcella_sram_tech_node[0]);
    Wmemcellpmos_sram = Wmemcellpmos_sram_tech_node[0] + alpha * (Wmemcellpmos_sram_tech_node[1] - Wmemcellpmos_sram_tech_node[0]);
    Wmemcellnmos_sram = Wmemcellnmos_sram_tech_node[0] + alpha * (Wmemcellnmos_sram_tech_node[1] - Wmemcellnmos_sram_tech_node[0]);
    area_cell_sram = area_cell_sram_tech_node[0] + alpha * (area_cell_sram_tech_node[1] - area_cell_sram_tech_node[0]);
    asp_ratio_cell_sram = asp_ratio_cell_sram_tech_node[0] + alpha * (asp_ratio_cell_sram_tech_node[1] - 
            asp_ratio_cell_sram_tech_node[0]);

    //Properties of transistors used within the DRAM cells and those used within the DRAM wordline 
    //driver circuitry. We assume that these transistors follow the LSTP roadmap i.e. their 
    //characteristics are those of the LSTP transistors of the same ITRS technology node/year. 
    //For the DRAM access transistor, the length and width of the transistor is based on IBM's 65 nm 
    //scalable DRAM cell.

    vdd_dram_cell = vdd_dram_cell_tech_node[0] + 
        alpha * (vdd_dram_cell_tech_node[1] - vdd_dram_cell_tech_node[0]);//V
    v_th_dram_access_transistor = v_th_dram_access_transistor_tech_node[0] + 
        alpha * (v_th_dram_access_transistor_tech_node[1] - v_th_dram_access_transistor_tech_node[0]);//V
    Lphy_dram_access_transistor = Lphy_dram_access_transistor_tech_node[0] + 
        alpha * (Lphy_dram_access_transistor_tech_node[1] - Lphy_dram_access_transistor_tech_node[0]);//micron. 
    Lelec_dram_access_transistor = Lelec_dram_access_transistor_tech_node[0] + 
        alpha * (Lelec_dram_access_transistor_tech_node[1] - Lelec_dram_access_transistor_tech_node[0]);//micron.
    c_g_ideal_itrs_dram_access_transistor = c_g_ideal_itrs_dram_access_transistor_tech_node[0] + 
        alpha * (c_g_ideal_itrs_dram_access_transistor_tech_node[1] - c_g_ideal_itrs_dram_access_transistor_tech_node[0]);//F/micron
    c_fringe_itrs_dram_access_transistor = c_fringe_itrs_dram_access_transistor_tech_node[0] + 
        alpha * (c_fringe_itrs_dram_access_transistor_tech_node[1] - c_fringe_itrs_dram_access_transistor_tech_node[0]);//F/micron
    c_junc_itrs_dram_access_transistor = c_junc_itrs_dram_access_transistor_tech_node[0] + 
        alpha * (c_junc_itrs_dram_access_transistor_tech_node[1] - c_junc_itrs_dram_access_transistor_tech_node[0]);//F/micron
    c_overlap_itrs_dram_access_transistor = (((Lphy_dram_access_transistor - Lelec_dram_access_transistor) / Lphy_dram_access_transistor) / 2) * c_g_ideal_itrs_dram_access_transistor;//F/micron
    I_on_n_dram_access_transistor = I_on_n_dram_access_transistor_tech_node[0] + 
        alpha * (I_on_n_dram_access_transistor_tech_node[1] - I_on_n_dram_access_transistor_tech_node[0]);//A/micron
    I_on_p_dram_access_transistor = I_on_n_dram_access_transistor / 2;//A/micron
    Rnchannelon_itrs_dram_access_transistor = vdd_dram_cell / I_on_n_dram_access_transistor;//ohm-micron
    Rpchannelon_itrs_dram_access_transistor = vdd_dram_cell / I_on_p_dram_access_transistor;//ohm-micron
    Wmemcella_dram = Wmemcella_dram_tech_node[0] + alpha * (Wmemcella_dram_tech_node[1] - Wmemcella_dram_tech_node[0]);
    Wmemcellpmos_dram = Wmemcellpmos_dram_tech_node[0] + alpha * (Wmemcellpmos_dram_tech_node[1] - Wmemcellpmos_dram_tech_node[0]);
    Wmemcellnmos_dram = Wmemcellnmos_dram_tech_node[0] + alpha * (Wmemcellnmos_dram_tech_node[1] - Wmemcellnmos_dram_tech_node[0]);
    area_cell_dram = area_cell_dram_tech_node[0] + alpha * (area_cell_dram_tech_node[1] - area_cell_dram_tech_node[0]);
    asp_ratio_cell_dram = asp_ratio_cell_dram_tech_node[0] + alpha * (asp_ratio_cell_dram_tech_node[1] - 
            asp_ratio_cell_dram_tech_node[0]);
    c_dram_cell = c_dram_cell_tech_node[0] +  alpha * (c_dram_cell_tech_node[1] - c_dram_cell_tech_node[0]);//V 


    vpp = vpp_tech_node[0] +  alpha * (vpp_tech_node[1] - vpp_tech_node[0]);//V 
    Lphy_dram_wordline_transistor = Lphy_dram_wordline_transistor_tech_node[0] + 
        alpha * (Lphy_dram_wordline_transistor_tech_node[1] - Lphy_dram_wordline_transistor_tech_node[0]);//micron. 
    Lelec_dram_wordline_transistor = Lelec_dram_wordline_transistor_tech_node[0] + 
        alpha * (Lelec_dram_wordline_transistor_tech_node[1] - Lelec_dram_wordline_transistor_tech_node[0]);//micron.
    c_g_ideal_itrs_dram_wordline_transistor = c_g_ideal_itrs_dram_wordline_transistor_tech_node[0] + 
        alpha * (c_g_ideal_itrs_dram_wordline_transistor_tech_node[1] - c_g_ideal_itrs_dram_wordline_transistor_tech_node[0]);//F/micron
    c_fringe_itrs_dram_wordline_transistor = c_fringe_itrs_dram_wordline_transistor_tech_node[0] + 
        alpha * (c_fringe_itrs_dram_wordline_transistor_tech_node[1] - c_fringe_itrs_dram_wordline_transistor_tech_node[0]);//F/micron
    c_junc_itrs_dram_wordline_transistor = c_junc_itrs_dram_wordline_transistor_tech_node[0] + 
        alpha * (c_junc_itrs_dram_wordline_transistor_tech_node[1] - c_junc_itrs_dram_wordline_transistor_tech_node[0]);//F/micron
    c_overlap_itrs_dram_wordline_transistor = (((Lphy_dram_wordline_transistor - Lelec_dram_wordline_transistor) / Lphy_dram_wordline_transistor) / 2) * c_g_ideal_itrs_dram_wordline_transistor;//F/micron
    I_on_n_dram_wordline_transistor = I_on_n_dram_wordline_transistor_tech_node[0] + 
        alpha * (I_on_n_dram_wordline_transistor_tech_node[1] - I_on_n_dram_wordline_transistor_tech_node[0]);//A/micron
    I_on_p_dram_wordline_transistor = I_on_n_dram_wordline_transistor / 2;//A/micron
    Rnchannelon_itrs_dram_wordline_transistor = vpp / I_on_n_dram_wordline_transistor;//ohm-micron
    Rpchannelon_itrs_dram_wordline_transistor = vpp / I_on_p_dram_wordline_transistor;//ohm-micron
    for(i = 0; i <= 100; i += 10){
        I_off_n_dram_wordline_transistor[i] = I_off_n_dram_wordline_transistor_tech_node[0][i] + 
            alpha * (I_off_n_dram_wordline_transistor_tech_node[1][i] - I_off_n_dram_wordline_transistor_tech_node[0][i] );//A/micron
        I_off_p_dram_wordline_transistor[i] = I_off_n_dram_wordline_transistor[i]; 
    }

    //Sense amplifier latch Gm calculation
    mobility_eff_periph_global = mobility_eff_periph_global_tech_node[0] + 
        alpha * (mobility_eff_periph_global_tech_node[1] - mobility_eff_periph_global_tech_node[0]);//micron2 / Vs
    Vdsat_periph_global = Vdsat_periph_global_tech_node[0] + 
        alpha * (Vdsat_periph_global_tech_node[1] - Vdsat_periph_global_tech_node[0]);//V/micron
    gmn_sense_amp_latch = (mobility_eff_periph_global / 2) * c_ox_periph_global * (WsenseN / Lelec_periph_global) * Vdsat_periph_global; //gm of latch NMOS transistor.
    gmp_sense_amp_latch = gmn_sense_amp_latch;//gm of latch PMOS transistor. We assume that mobility is half that of NMOS and width is 
    //twice that of NMOS, so gmp equal to gmn. 
    Gm_sense_amp_latch = gmn_sense_amp_latch + gmp_sense_amp_latch;

    BitWidth_dram = sqrt(area_cell_dram / (asp_ratio_cell_dram));
    BitHeight_dram = asp_ratio_cell_dram * BitWidth_dram;
    BitWidth_sram = sqrt(area_cell_sram / (asp_ratio_cell_sram));
    BitHeight_sram = asp_ratio_cell_sram * BitWidth_sram;

    Vt_dram = v_th_dram_access_transistor;
    Vbitpre_dram = vdd_dram_cell; //VDD (GND) precharge
    Vt_sram = v_th[sram_cell_and_wordline_tech_flavor];
    Vbitpre_sram = vdd[sram_cell_and_wordline_tech_flavor];


    //Interconnect parameters
    if(technology < 101 && technology > 99){
        tech_low = 100;
        tech_high = 100;
    }
    else if(technology < 71 && technology > 69){
        tech_low = 70;
        tech_high = 70;
    }
    else if(technology < 51 && technology > 49){
        tech_low = 50;
        tech_high = 50;
    }
    else if(technology < 36 && technology > 34){
        tech_low = 35;
        tech_high = 35;
    }
    else if(technology < 26 && technology > 24){
        tech_low = 25;
        tech_high = 25;
    }
    else if(technology < 100 && technology > 70){
        tech_low = 100;
        tech_high = 70;
    }
    else if(technology < 70 && technology > 50){
        tech_low = 70;
        tech_high = 50;
    }
    else if(technology < 50 && technology > 35){
        tech_low = 50;
        tech_high = 35;
    }
    else if(technology < 35 && technology > 25){
        tech_low = 35;
        tech_high = 25;
    }

    for(iter=0; iter<=1; ++iter){
        if(iter==0){
            tech = tech_low;
        }
        else{
            tech = tech_high;
        }

        wire_c_per_micron[0][0] = wire_cap2(0 , 0, technology);//F/micron.
        wire_c_per_micron[0][1] = wire_cap2(0 , 1, technology);//F/micron. 
        wire_c_per_micron[0][2] = wire_cap2(0 , 2, technology);//F/micron. 
        wire_c_per_micron[1][0] = wire_cap2(1 , 0, technology);//F/micron.
        wire_c_per_micron[1][1] = wire_cap2(1 , 1, technology);//F/micron. 
        wire_c_per_micron[1][2] = wire_cap2(1 , 2, technology);//F/micron. 
        if(tech < 101 && tech > 99){
            //Interconnect parameters from Ron Ho's thesis. 
            //Aggressive projections.
            wire_pitch[0][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[0][1] = 4 * FEATURESIZE;//micron
            wire_pitch[0][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[0][0] = 0.363e-15;//F/micron.
            //					  wire_c_per_micron[0][1] = 0.374e-15;//F/micron. 
            //					  wire_c_per_micron[0][2] = 0.403e-15;//F/micron. 
            wire_r_per_micron[0][0] = 0.72; //ohm/micron
            wire_r_per_micron[0][1] = 0.26; //ohm/micron
            wire_r_per_micron[0][2] = 0.054; //ohm/micron
            //Conservative projections
            wire_pitch[1][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[1][1] = 4 * FEATURESIZE;//micron
            wire_pitch[1][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[1][0] = 0.349e-15;//F/micron. 
            //					  wire_c_per_micron[1][1] = 0.359e-15;//F/micron. 
            //					  wire_c_per_micron[1][2] = 0.377e-15;//F/micron. 
            wire_r_per_micron[1][0] = 0.83; //ohm/micron
            wire_r_per_micron[1][1] = 0.307; //ohm/micron
            wire_r_per_micron[1][2] = 0.073; //ohm/micron
        }
        else if(tech < 71 && tech > 69){
            //Interconnect parameters from Ron Ho's thesis. 
            //70 nm Aggressive projections. 
            wire_pitch[0][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[0][1] = 4 * FEATURESIZE;//micron
            wire_pitch[0][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[0][0] = 0.35e-15;//F/micron.
            //					  wire_c_per_micron[0][1] = 0.359e-15;//F/micron. 
            //					  wire_c_per_micron[0][2] = 0.367e-15;//F/micron. 
            wire_r_per_micron[0][0] = 0.87; //ohm/micron
            wire_r_per_micron[0][1] = 0.34; //ohm/micron
            wire_r_per_micron[0][2] = 0.082; //ohm/micron
            //70 nm Conservative projections
            wire_pitch[1][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[1][1] = 4 * FEATURESIZE;//micron
            wire_pitch[1][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[1][0] = 0.324e-15;//F/micron.
            //					  wire_c_per_micron[1][1] = 0.333e-15;//F/micron. 
            //					  wire_c_per_micron[1][2] = 0.353e-15;//F/micron.
            wire_r_per_micron[1][0] = 1.72; //ohm/micron
            wire_r_per_micron[1][1] = 0.627; //ohm/micron
            wire_r_per_micron[1][2] = 0.15; //ohm/micron
        }
        else if(tech < 51 && tech > 49){
            //Interconnect parameters from Ron Ho's thesis. 
            //Aggressive projections. 
            wire_pitch[0][0] = 2 * FEATURESIZE;//micron
            wire_pitch[0][1] = 4 * FEATURESIZE;//micron
            wire_pitch[0][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[0][0] = 0.337e-15;//F/micron.
            //					  wire_c_per_micron[0][1] = 0.345e-15;//F/micron. 
            //					  wire_c_per_micron[0][2] = 0.345e-15;//F/micron. 
            wire_r_per_micron[0][0] = 1.54; //ohm/micron
            wire_r_per_micron[0][1] = 0.6; //ohm/micron
            wire_r_per_micron[0][2] = 0.15; //ohm/micron
            //Conservative projections
            wire_pitch[1][0] = 2 * FEATURESIZE;//micron
            wire_pitch[1][1] = 4 * FEATURESIZE;//micron
            wire_pitch[1][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[1][0] = 0.303e-15;//F/micron.
            //					  wire_c_per_micron[1][1] = 0.311e-15;//F/micron. 
            //					  wire_c_per_micron[1][2] = 0.332e-15;//F/micron. 
            wire_r_per_micron[1][0] = 3.34; //ohm/micron
            wire_r_per_micron[1][1] = 1.22; //ohm/micron
            wire_r_per_micron[1][2] = 0.292; //ohm/micron

        }
        else if(tech < 36 && tech > 34){
            //Interconnect parameters from Ron Ho's thesis. 
            //Aggressive projections. 
            wire_pitch[0][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[0][1] = 4 * FEATURESIZE;//micron
            wire_pitch[0][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[0][0] = 0.306e-15;//F/micron. 
            //					  wire_c_per_micron[0][1] = 0.315e-15/2;//F/micron. 
            //					  wire_c_per_micron[0][2] = 0.315e-15;//F/micron. 
            wire_r_per_micron[0][0] = 3.13; //ohm/micron
            wire_r_per_micron[0][1] = 1.224; //ohm/micron
            wire_r_per_micron[0][2] = 0.306; //ohm/micron
            //Conservative projections
            wire_pitch[1][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[1][1] = 4 * FEATURESIZE;//micron
            wire_pitch[1][2] = 8 * FEATURESIZE;//micron
            //					  wire_c_per_micron[1][0] = 0.286e-15;//F/micron.
            //					  wire_c_per_micron[1][1] = 0.295e-15/2;//F/micron. 
            //					  wire_c_per_micron[1][2] = 0.313e-15;//F/micron. 
            wire_r_per_micron[1][0] = 6.89; //ohm/micron
            wire_r_per_micron[1][1] = 2.509; //ohm/micron
            wire_r_per_micron[1][2] = 0.598; //ohm/micron
        }

        else if(tech < 26 && tech > 24){
            //Interconnect parameters from Ron Ho's thesis.
            //Aggressive projections. 
            wire_pitch[0][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[0][1] = 4 * FEATURESIZE;//micron
            wire_pitch[0][2] = 8 * FEATURESIZE;//micron
            wire_c_per_micron[0][0] = 0.28e-15;//F/micron. 
            wire_c_per_micron[0][1] = 0.288e-15;//F/micron. 
            wire_c_per_micron[0][2] = 0.288e-15;//F/micron. 
            wire_r_per_micron[0][0] = 6.14; //ohm/micron
            wire_r_per_micron[0][1] = 2.4; //ohm/micron
            wire_r_per_micron[0][2] = 0.6; //ohm/micron
            //Conservative projections
            wire_pitch[1][0] = 2.5 * FEATURESIZE;//micron
            wire_pitch[1][1] = 4 * FEATURESIZE;//micron
            wire_pitch[1][2] = 8 * FEATURESIZE;//micron
            wire_c_per_micron[1][0] = 0.278e-15;//F/micron.
            wire_c_per_micron[1][1] = 0.287e-15;//F/micron. 
            wire_c_per_micron[1][2] = 0.288e-15;//F/micron. 
            wire_r_per_micron[1][0] = 14.64; //ohm/micron
            wire_r_per_micron[1][1] = 5.15; //ohm/micron
            wire_r_per_micron[1][2] = 1.2; //ohm/micron
        }

        //Interconnect parameters. For wires inside and outside a mat. The projection type
        //can be aggressive or conservative (Ron Ho's projections). And the wires inside or outside a 
        //mat may be mapped to semi-global or global wire types. 
        wire_local_pitch_tech_node[iter] = wire_pitch[interconnect_projection_type][0];
        wire_local_r_per_micron_tech_node[iter] = wire_r_per_micron[interconnect_projection_type][0];
        wire_local_c_per_micron_tech_node[iter] = wire_c_per_micron[interconnect_projection_type][0];
        wire_inside_mat_pitch_tech_node[iter] = wire_pitch[interconnect_projection_type][wire_inside_mat_type];
        wire_inside_mat_r_per_micron_tech_node[iter] = wire_r_per_micron[interconnect_projection_type][wire_inside_mat_type];
        wire_inside_mat_c_per_micron_tech_node[iter] = wire_c_per_micron[interconnect_projection_type][wire_inside_mat_type];
        wire_outside_mat_pitch_tech_node[iter] = wire_pitch[interconnect_projection_type][wire_outside_mat_type];
        wire_outside_mat_r_per_micron_tech_node[iter] = wire_r_per_micron[interconnect_projection_type][wire_outside_mat_type];
        wire_outside_mat_c_per_micron_tech_node[iter] = wire_c_per_micron[interconnect_projection_type][wire_outside_mat_type];
    }

    //Interconnect parameters. For wires inside and outside a mat. The projection type
    //can be aggressive or conservative (Ron Ho's projections). And the wires inside or outside a 
    //mat may be mapped to semi-global or global wire types. 

    int_tech_low = (int) (floor(tech_low + 0.5));
    int_tech_high = (int) (floor(tech_high+ 0.5));
    if(int_tech_low != int_tech_high){
        alpha = (technology - tech_low) / (tech_high - tech_low);
    }
    else{
        alpha =0;
    }

    wire_local_pitch = wire_local_pitch_tech_node[0] +  alpha * (wire_local_pitch_tech_node[1] - 
            wire_local_pitch_tech_node[0]);
    wire_local_r_per_micron = wire_r_per_micron[interconnect_projection_type][0];
    wire_local_r_per_micron = wire_local_r_per_micron_tech_node[0] +  alpha * 
        (wire_local_r_per_micron_tech_node[1] - wire_local_r_per_micron_tech_node[0]);
    wire_local_c_per_micron = wire_c_per_micron[interconnect_projection_type][0];
    wire_local_c_per_micron = wire_local_c_per_micron_tech_node[0] +  alpha * 
        (wire_local_c_per_micron_tech_node[1] - wire_local_c_per_micron_tech_node[0]);
    wire_inside_mat_pitch = wire_inside_mat_pitch_tech_node[0] +  alpha * (wire_inside_mat_pitch_tech_node[1] - 
            wire_inside_mat_pitch_tech_node[0]);
    wire_inside_mat_r_per_micron = wire_inside_mat_r_per_micron_tech_node[0] +  alpha * 
        (wire_inside_mat_r_per_micron_tech_node[1] - wire_inside_mat_r_per_micron_tech_node[0]);
    wire_inside_mat_c_per_micron = wire_inside_mat_c_per_micron_tech_node[0] +  alpha * 
        (wire_inside_mat_c_per_micron_tech_node[1] - wire_inside_mat_c_per_micron_tech_node[0]);
    wire_outside_mat_pitch = wire_outside_mat_pitch_tech_node[0] +  alpha * (wire_outside_mat_pitch_tech_node[1] - 
            wire_outside_mat_pitch_tech_node[0]);
    wire_outside_mat_r_per_micron = wire_outside_mat_r_per_micron_tech_node[0] +  alpha * 
        (wire_outside_mat_r_per_micron_tech_node[1] - wire_outside_mat_r_per_micron_tech_node[0]);
    wire_outside_mat_c_per_micron = wire_outside_mat_c_per_micron_tech_node[0] +  alpha * 
        (wire_outside_mat_c_per_micron_tech_node[1] - wire_outside_mat_c_per_micron_tech_node[0]);
}
