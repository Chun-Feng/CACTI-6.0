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

#include "io.h"

#define NEXTINT(a) skip(); scanf("%d",&(a));
#define NEXTFLOAT(a) skip(); scanf("%lf",&(a));
/*---------------------------------------------------------------*/

						
/*
 * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C"
 */
void strreverse(char* begin, char* end) {
    char aux;
    while(end>begin)
        aux=*end, *end--=*begin, *begin++=aux;
}

#ifdef __linux__
void itoa(int value, char* str, int base) {
    static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    char* wstr=str;
    int sign;
    div_t res;

    // Validate base
    if (base<2 || base>35){ *wstr='\0'; return; }

    // Take care of sign
    if ((sign=value) < 0) value = -value;

    // Conversion. Number is reversed.
    do {
        res = div(value,base);
        *wstr++ = num[res.rem];
    }
    while(value=res.quot);
    if(sign<0) *wstr++='-';
    *wstr='\0';

    // Reverse string
    strreverse(str,wstr-1);
}
#endif

int
get_cpucount()
{
    #define CPUINFO    "/proc/cpuinfo"

    int r, proc = 0;
#ifdef __linux__
    char line[5000];
    FILE *fp = fopen(CPUINFO, "r");

    if (fp == NULL) {
        fprintf(stderr, "Unable to open %s\n", CPUINFO);
        /* have atleast 2 threads */
        return 2;
    }

    while(fscanf(fp, "%[^\n]\n", line) != EOF) {
        if (!strncmp("processor", line, strlen("processor")))
            if (sscanf(line, "processor : %d", &r))
                proc++;
    }
    fclose(fp);

    if (!proc)
        return 2;
    return proc;
#endif
/* have atleast 2 threads */
    proc = 2;
    return proc;
}

void
read_file() 
{
    FILE *fp;
    int i;
    fp = fopen ("cummulative.dat", "r");
    if (!fp) fprintf(stderr, "Error opening file cummulative.dat!!\n");
    for (i=0; i<1023; i++) 
    {
        fscanf (fp, "%f", &cumm_per[i]);
    }
    cumm_per[1023] = 100.00;
    fclose(fp);
}

/* Parse the "cache.cfg" config file */
void
parse_cfg(input_params_t *params)
{
    FILE *fp = fopen("cache.cfg", "r");
    char line[5000];
    char jk[5000];
    char temp[5000];

    if(!fp) {
        fprintf(stderr, "cache.cfg file missing!\n");
        exit(-1);
    }

    while(fscanf(fp, "%[^\n]\n", line) != EOF) {

        if (!strncmp("-size", line, strlen("-size"))) {
            sscanf(line, "-size %[(:-~)*]%ld", jk, &(params->cache_size));
            continue;
        }

        if (!strncmp("-block", line, strlen("-block"))) {
            sscanf(line, "-block size (bytes) %d", &(params->block_size)); 
            continue;
        }

        if (!strncmp("-associativity", line, strlen("-associativity"))) {
            sscanf(line, "-associativity %d", &(params->associativity));
            continue;
        }

        if (!strncmp("-read-write", line, strlen("-read-write"))) {
            sscanf(line, "-read-write port %d", &(params->rw_ports));
            continue;
        }

        if (!strncmp("-exclusive read", line, strlen("exclusive read"))) {
            sscanf(line, "-exclusive read port %d", &(params->excl_read_ports));
            continue;
        }

        if(!strncmp("-exclusive write", line, strlen("-exclusive write"))) {
            sscanf(line, "-exclusive write port %d", &(params->excl_write_ports));
            continue;
        }

        if (!strncmp("-single ended", line, strlen("-single ended"))) {
            sscanf(line, "-single %[(:-~)*]%d", jk, 
                            &(params->single_ended_read_ports));
            continue;
        }

        if (!strncmp("-UCA bank", line, strlen("-UCA bank"))) {
            sscanf(line, "-UCA bank%[((:-~)| )*]%d", jk, &(params->uca_banks));
            continue;
        }
        
        if (!strncmp("-tech", line, strlen("-tech"))) {
            sscanf(line, "-technology (u) %lf", &(params->tech_size));
            continue;
        }

        if (!strncmp("-output/input", line, strlen("-output/input"))) {
            sscanf(line, "-output/input bus width %d", &(params->output_width));
            continue;
        }

        if (!strncmp("-operating", line, strlen("-operating"))) {
            sscanf(line, "-operating temperature (K) %d", &(params->temp));
			temper = params->temp;
            continue;
        }

        if (!strncmp("-tag size", line, strlen("-tag size"))) {
            sscanf(line, "-tag size%[^\"]\"%[^\"]\"", jk, temp);
            if (!strncmp("default", temp, sizeof("default"))) {
                params->tag_size = 42; /* the acutal value is calculated
                * later based on the cache size, bank count, and associativity
                */
            }
            else {
                sscanf(line, "-tag size (b) %d", &(params->tag_size));
                params->force_tag = 1;
            }
            continue;
        }

        if (!strncmp("-access mode", line, strlen("-access mode"))) {
            sscanf(line, "-access %[^\"]\"%[^\"]\"", jk, temp);
            if (!strncmp("fast", temp, strlen("fast"))) {
                params->access_mode = 2;
            }
            else if (!strncmp("sequential", temp, strlen("sequential"))) {
                params->access_mode = 0;
            }
            else if(!strncmp("normal", temp, strlen("normal"))) {
                params->access_mode = 1;
            }
            else {
                printf("ERROR: Invalid access mode!\n");
            }
            continue;
        }

        if (!strncmp("-cache type", line, strlen("-cache type"))) {
            sscanf(line, "-cache type %[^\"]\"%[^\"]\"", jk, temp);

            if(!strncmp("SRAM_CACHE", temp, strlen("SRAM_CACHE"))) {
                params->pure_sram = 0;
            }
            else {
                params->pure_sram = 1;
            }

            if(!strncmp("DRAM", temp, strlen("DRAM"))) {
                params->dram = 1;
                params->pure_sram = 0;
            }
            else {
                params->dram = 0;
            }
            continue;
        }

        if(!strncmp("-design", line, strlen("-design"))) {
            sscanf(line, "-%[((:-~)| |,)*]%d:%d:%d:%d:%d", jk,
                &(params->delay_wt), &(params->dynamic_power_wt), 
                &(params->leakage_power_wt),
                &(params->cycle_time_wt), &(params->area_wt));
            continue;
        }

        if(!strncmp("-deviate", line, strlen("-deviate"))) {
            sscanf(line, "-%[((:-~)| |,)*]%d:%d:%d:%d:%d", jk,
                &(params->delay_dev), &(params->dynamic_power_dev), 
                &(params->leakage_power_dev),
                &(params->cycle_time_dev), &(params->area_dev));
            continue;
        }

        if(!strncmp("-Optimize", line, strlen("-Optimize"))) {
            sscanf(line, "-Optimize  %[^\"]\"%[^\"]\"", jk, temp);

            if(!strncmp("ED^2", temp, strlen("ED^2"))) {
                params->ed = 2;
            }
            else if(!strncmp("ED", temp, strlen("ED"))) {
                params->ed = 1;
            }
            else {
                params->ed = 0;
            }
        }

        if(!strncmp("-NUCAdesign", line, strlen("-NUCAdesign"))) {
            sscanf(line, "-%[((:-~)| |,)*]%d:%d:%d:%d:%d", jk,
                &(params->delay_wt_nuca), &(params->dynamic_power_wt_nuca), 
                &(params->leakage_power_wt_nuca),
                &(params->cycle_time_wt_nuca), &(params->area_wt_nuca));
            continue;
        }

        if(!strncmp("-NUCAdeviate", line, strlen("-NUCAdeviate"))) {
            sscanf(line, "-%[((:-~)| |,)*]%d:%d:%d:%d:%d", jk,
                &(params->delay_dev_nuca), &(params->dynamic_power_dev_nuca), 
                &(params->leakage_power_dev_nuca),
                &(params->cycle_time_dev_nuca), &(params->area_dev_nuca));
            continue;
        }

        if(!strncmp("-Cache model", line, strlen("-cache model"))) {
            sscanf(line, "-Cache model %[^\"]\"%[^\"]\"", jk, temp);

            if (!strncmp("UCA", temp, strlen("UCA"))) {
                NUCA = 0;
            }
            else {
                NUCA = 1;
            }
            continue;
        }

        if(!strncmp("-NUCA bank", line, strlen("-NUCA bank"))) {
            sscanf(line, "-NUCA bank count %d", &(params->nuca_bank_count));
            if (params->nuca_bank_count != 0) {
                params->force_nuca_bank = 1;
            }
            continue;
        }

        if(!strncmp("-wire", line, strlen("-wire"))) {
            sscanf(line, "-wire%[^\"]\"%[^\"]\"", jk, temp);

            if (!strncmp("default", temp, strlen("default"))) {
                params->force_wiretype = 0;
                continue;
            }
            else if (!(strncmp("fullswing", temp, strlen("fullswing")))) {
                params->wire_inter_mats = 1;
            }
            else {
                params->wire_inter_mats = 0;
            }
            params->force_wiretype = 1;
            continue;
        }

        if(!strncmp("-Core", line, strlen("-Core"))) {
            sscanf(line, "-Core count %d\n", &(params->cores));
            if (params->cores > 16) {
                printf("No. of cores should be less than 16!\n");
            }
            continue;
        }

        if(!strncmp("-Cache level", line, strlen("-Cache level"))) {
            sscanf(line, "-Cache l%[^\"]\"%[^\"]\"", jk, temp);
            if (!strncmp("L2", temp, strlen("L2"))) {
                params->cache_level = 0;
            }
            else {
                params->cache_level = 1;
            }
        }

        if(!strncmp("-Print level", line, strlen("-Print level"))) {
            sscanf(line, "-Print l%[^\"]\"%[^\"]\"", jk, temp);
            if (!strncmp("DETAILED", temp, strlen("DETAILED"))) {
                params->print_detail = 1;
            }
            else {
                params->print_detail = 0;
            }
        }
    }
    fclose(fp);
}


int 
parse_cmd_args(int argc,char *argv[], input_params_t *params)
{
    int i=0;


	if(argc > 1) {
        if(!strcmp(argv[1], "--h") ||
                !strcmp(argv[1], "-h") ||
                !strcmp(argv[1], "-help") ||
                !strcmp(argv[1], "--help")) {
            goto invalid_args;
        }
       if(!strcmp(argv[1], "-weight")) {
           if (argc < 7) goto invalid_args;
           params->delay_wt = atoi(argv[2]);
           params->dynamic_power_wt = atoi(argv[3]);
           params->leakage_power_wt = atoi(argv[4]);
           params->cycle_time_wt = atoi(argv[5]);
           params->area_wt = atoi(argv[6]);
           if(argv[7] && !strcmp(argv[7], "-deviate")) {
               i = 8;
               if (argc < 13) goto invalid_args;
               params->delay_dev = atoi(argv[i]);
               params->dynamic_power_dev = atoi(argv[i+1]);
               params->leakage_power_dev = atoi(argv[i+2]);
               params->cycle_time_dev = atoi(argv[i+3]);
               params->area_dev = atoi(argv[i+4]);
           }
       }
       else if(!strcmp(argv[1], "-deviate")) {
           i = 2;
           if (argc < 7) goto invalid_args;
           params->delay_dev = atoi(argv[i]);
           params->dynamic_power_dev = atoi(argv[i+1]);
           params->leakage_power_dev = atoi(argv[i+2]);
           params->cycle_time_dev = atoi(argv[i+3]);
           params->area_dev = atoi(argv[i+4]);
       }
       else if (argc >= 6) {
           params->cache_size = atoi(argv[1]);
           params->block_size = atoi(argv[2]);
           params->associativity = atoi(argv[3]);
           params->tech_size = atof(argv[4]);
           /* 
            * Till CACTI 4 tech size input unit is microns and 
            * in CACTI 5 it is nm. 
            * Check the input tech size and change it to microns
            * if necessary
            */
           if(params->tech_size > 1) {
               params->tech_size /= 1000; /* If the entry is invalid
                                           * validate_cache_args will throw an error 
                                           */
           }
           params->uca_banks = atoi(argv[5]);
       }
       else {
           goto invalid_args;
       }
    }
	if(argc > 6) {
        if(!strcmp(argv[6], "-weight")) {
            if (argc < 12) goto invalid_args;
            params->delay_wt = atoi(argv[7]);
            params->dynamic_power_wt = atoi(argv[8]);
            params->leakage_power_wt = atoi(argv[9]);
            params->cycle_time_wt = atoi(argv[10]);
            params->area_wt = atoi(argv[11]);
        }
       if(!strcmp(argv[6], "-deviate")) {
           i = 7;
           if (argc < 12) goto invalid_args;
           params->delay_dev = atoi(argv[i]);
           params->dynamic_power_dev = atoi(argv[i+1]);
           params->leakage_power_dev = atoi(argv[i+2]);
           params->cycle_time_dev = atoi(argv[i+3]);
           params->area_dev = atoi(argv[i+4]);
       }
       else if(argv[12]) {
           if(!strcmp(argv[12], "-deviate")) {
               i = 13;
               if (argc < 17) goto invalid_args;
               params->delay_dev = atoi(argv[i]);
               params->dynamic_power_dev = atoi(argv[i+1]);
               params->leakage_power_dev = atoi(argv[i+2]);
               params->cycle_time_dev = atoi(argv[i+3]);
               params->area_dev = atoi(argv[i+4]);
           }
       }
   }
   return 0;
invalid_args:
    printf("\nThe command line arguments are optional. Instead, the " 
            "cache.cfg file can be used to specify cache "
            "parameters\n\n");

    printf("Valid args:\n\t C B A Tech NoBanks "
            "\nand / or\n"
            "\t-weight <delay> <dynamic> <leakage> <cycle> <area> "
            "\nand / or\n"
            "\t-deviate <delay> <dynamic> <leakage> <cycle> <area>\n");
    printf("\nC - Cache size in bytes\nB - Block size in bytes\n"
           "A - Associativity\nTech - Process technology in u or nm\n"
           "NoBanks - No. of banks\n\n");
    exit(-1);
}

void
sim_cache(int argc,char *argv[])
{
	nuca_org_t *nuca_res;
    input_params_t *params = (input_params_t *) malloc(sizeof(input_params_t));
    memset(params, 0, sizeof(input_params_t));

    /*
     * Read default values from the config file 
     */
    parse_cfg(params); 

    /* 
     * Update input parameters with the command line
     * arguments
     */
    parse_cmd_args(argc, argv, params);

    /* 
     * check for any inconsistencies in the input
     * arguments
     */
    validate_cache_args(params);


    PRINTD(dump_input_args(params));

    /*
     * Initialize all process technology parameteres
     */
    init_tech_params(params->tech_size);

//#define study_wires
#ifdef study_wires
    int i;
    wire_stats_t wst[5]; 
    for(i=0; i<5; i++) {
        wst[i].wire_length = 1e-3; 
        wst[i].nsense = 1;
    }
    wst[0].wt = Global; 
    calc_wire_stats2(Global,&(wst[0])); 
    wst[1].wt = Global_10;
    calc_wire_stats2(Global_10,&(wst[1])); 
    wst[2].wt = Global_20;
    calc_wire_stats2(Global_20,&(wst[2])); 
    wst[3].wt = Global_30;
    calc_wire_stats2(Global_30,&(wst[3])); 
    wst[4].wt = Low_swing;
    calc_wire_stats2(Low_swing,&(wst[4])); 
    for(i=0; i<5; i++) {
        print_wire(&wst[i]);
    }
    exit(0);
#endif

    if(NUCA) {
        /*
         * Initialize router parameters
         */
        init_router_params();

        nuca_res = (nuca_org_t *) malloc(sizeof(nuca_org_t));
        memset(nuca_res, 0, sizeof(nuca_org_t));
        nuca_res->params = params;

        sim_nuca(nuca_res);

        output_NUCA(nuca_res);
        params->cache_size = nuca_res->params->cache_size/nuca_res->bank_count;

        free(nuca_res);
    }
    uca_org_t *uca_res = (uca_org_t *) malloc(sizeof(uca_org_t));
    results_mem_array *tag = (results_mem_array *)
        malloc(sizeof(results_mem_array));
    results_mem_array *data = (results_mem_array *)
        malloc(sizeof(results_mem_array));
    memset(uca_res, 0, sizeof(uca_res));
    memset(tag, 0, sizeof(results_mem_array));
    memset(data, 0, sizeof(results_mem_array));

    uca_res->tag_array = tag;
    uca_res->data_array = data;
    uca_res->params = params;


    sim_uca(uca_res);

    output_UCA (uca_res);

    free(tag);
    free(data);
    free(uca_res);
    free(params);
    free_router();
}

void
validate_cache_args (input_params_t *params)
{
   int B = params->block_size;
   long int C = params->cache_size/((int) (params->uca_banks));
   int RWP = params->rw_ports;
   int ERP = params->excl_read_ports;
   int EWP = params->excl_write_ports;
   int A = params->associativity;
   /* temp var */
   double logbanks, logbanksfloor;
   double assoc, assocfloor;

   if ((B < 1)) {
       printf("ERROR: Block size must >=1\n");
	   exit(1);
   }

   if ((B*8 < params->output_width)) {
       printf("WARNING: Bad design choice: block size is less"
              " than bus width\n");
   }
   
   if ((params->tech_size <= 0)) {
       printf("ERROR: Feature size must be > 0\n");
	   exit(1);
   }

   if ((params->tech_size > 0.09)) {
       printf("Feature size must be <= 0.09 (um)\n");
	   exit(1);
   }

   if ((RWP < 0) || (EWP < 0) || (ERP < 0)) {
       printf("Ports must >=0\n");
       exit(1);
   }
   if (RWP > 2) {
       printf("Maximum of 2 read/write ports\n");
       exit(1);
   }
   if ((RWP+ERP+EWP) < 1) {
       printf("Must have at least one port\n");
       exit(1);
   }

   if (params->uca_banks < 1 ) {
       printf("Number of subbanks should be greater "
               "than or equal to 1 and should be a power of 2\n");
       exit(1);
   }

   logbanks = logtwo((double)(params->uca_banks));
   logbanksfloor = floor(logbanks);

   if(logbanks > logbanksfloor){
       printf("Number of subbanks should be greater "
       "than or equal to 1 and should be a power of 2\n");
       exit(1);
   }

 
   if ((C < 64)) {
       printf("ERROR: Cache size must be >=64\n");
	   exit(1);
   }
 
   if (A < 0) {
             printf("ERROR: Associativity must >= 1\n");
   }
   else if (A == 0) {
       params->fully_assoc = 1;
       A = 1;
   }
   assoc = logtwo((double)(A));
   assocfloor = floor(assoc);

   if(assoc > assocfloor){
       printf("ERROR: Associativity should be a power of 2\n");
       exit(1);
   }

//           //if ((A > 32)) {
//             //printf("Associativity must <= 32\n or try FA (fully associative)\n");
//			 //exit(1);
//           //}
//
   if (C/(B*A)<=1 && !params->fully_assoc) {
     printf("Number of sets is too small:\n  Need to either "
            "increase cache size, or decrease associativity "
            "or block size\n  (or use fully associative cache)\n");
	exit(1); 
   }

   if(params->access_mode == 0) {
       params->tag_associativity = params->associativity;
       params->data_associativity = 1;
       params->sequential_access = 1;
   }
   else if(params->access_mode == 1) {
       params->tag_associativity = params->data_associativity = 
           params->associativity;
       params->sequential_access = 0;
   }
   else {
       params->fast_access = 1;
       params->tag_associativity = params->data_associativity = 
           params->associativity;
       params->sequential_access = 0;
   }
   if (params->tag_size != ADDRESS_BITS) {
       params->force_tag = 1;
   }
   if (params->fully_assoc) {
       params->data_associativity = 1;
   }
}

/*
 * Dumps a .fig file that shows the organization 
 * of the cache for
 * the given input parameters
 * Fig format for version 3.2 can be found at
 * http://www.xfig.org/userman/fig-format.html 
 */
void fig_out (final_results *fs)
{
    /* Cache input parameters and other parameters calculated by cacti */
    int wl, bl, nspd, size, blk_size, assoc;
    int ntwl, ntbl, ntspd; /* tag array param */
    int fa; /* fully assoc? */
    char temp[20];
    char filename[100];
    char cmd[300];
    FILE *fp;
    int cy;
    struct fig_mat_dim figdim;

    /* Usable area in the figure.
    If the page size is changed from letter to a3 or a2 then modify these
    values accordingly */

    int LEFT_x = (150+100);
    int RIGHT_x = 10075;

    ntwl = fs->tag_array.Ndwl;
    ntbl = fs->tag_array.Ndbl;
    ntspd = fs->tag_array.Nspd;

    wl = fs->data_array.Ndwl/2; /* each mat has four sub-arrays */
    bl = fs->data_array.Ndbl/2; /* and sub-arrays are grouped as four */
    nspd = (int) fs->data_array.Nspd;
    size = fs->params->cache_size;
    blk_size = fs->params->block_size;
    assoc = fs->params->data_associativity;
    fa = fs->params->fully_assoc;

    /* 
     * To use this functionality for online CACTI 
     * we should provide a unique output file name for
     * each run. The following fragment constructs a
     * name from the input arguments
     */
    itoa(size, filename, 10);
    strcpy(fs->file_n, filename);
    itoa(wl*2, filename, 10);
    strcat(fs->file_n, filename);
    itoa(bl*2, filename, 10);
    strcat(fs->file_n, filename);
    itoa(ntwl, filename, 10);
    strcat(fs->file_n, filename);
    itoa(ntbl, filename, 10);
    strcat(fs->file_n, filename);
    itoa(fs->params->data_associativity, filename, 10);
    strcat(fs->file_n, filename);
    itoa(fs->params->block_size, filename, 10);
    strcat(fs->file_n, filename);
    itoa(fs->params->temp, filename, 10);
    strcat(fs->file_n, filename);
    strcpy(filename, fs->file_n);
    strcat(filename, ".fig");
    strcat(fs->file_n, ".png");

    /* create a .fig file that shows data/tag array organization */
    printf("File name(s) = %s\n", fs->file_n);
    fp = fopen ((char *) &(filename), "w");
    if (!fp) fprintf (stderr, "Unable to create cache.fig file\n");
    /* .fig headers */
    fprintf(fp, "%s\n", "#FIG 3.2\nPortrait\nCenter\nInches\nLetter\n"
        "100.00\nSingle\n-2\n1200 2");

    /* show a sample Mat organization */
    print_mat (fp);

    /* print basic cache parameters */
    itoa (size, temp, 10);
    strcat (temp, " bytes\\001");
    fprintf (fp, "%s", "4 0 0 50 -1 1 16 0.0000 0 225 3165 750 675 Cache Size - ");
    fprintf (fp, "%s\n", temp);
    itoa (blk_size, temp, 10);
    strcat (temp, " bytes\\001");
    fprintf (fp, "%s", "4 0 0 50 -1 1 16 0.0000 0 225 3150 750 960 Block Size - ");
    fprintf (fp, "%s\n", temp);
    fprintf (fp, "%s", "4 0 0 50 -1 1 16 0.0000 0 225 2100 750 1275 Associativity - ");
    if (fa) {
        fprintf (fp, "%s\n", "Fully associative\\001");
    }
    else {
        itoa (assoc, temp, 10);
        strcat (temp, " \\001");
        fprintf (fp, "%s\n", temp);
    }
    fprintf (fp, "4 0 0 50 -1 1 14 0.0000 0 135 690 750 1575 Ndwl = ");
    itoa (fs->data_array.Ndwl, temp, 10);
    fprintf (fp, "%s\\001\n", temp);
    fprintf (fp, "4 0 0 50 -1 1 14 0.0000 0 135 645 750 1800 Ndbl = ");
    itoa (fs->data_array.Ndbl, temp, 10);
    fprintf (fp, "%s\\001\n", temp);
    fprintf (fp, "4 0 0 50 -1 1 14 0.0000 0 180 690 750 2025 Nspd = ");
    itoa (fs->data_array.Nspd, temp, 10);
    fprintf (fp, "%s\\001\n", temp);
    fprintf (fp, "4 0 0 50 -1 1 14 0.0000 0 135 660 750 2250 Ntwl = ");
    itoa (fs->data_array.Ndwl, temp, 10);
    fprintf (fp, "%s\\001\n", temp);
    fprintf (fp, "4 0 0 50 -1 1 14 0.0000 0 135 615 750 2475 Ntbl = ");
    itoa (fs->data_array.Ndbl, temp, 10);
    fprintf (fp, "%s\\001\n", temp);
    fprintf (fp, "4 0 0 50 -1 1 14 0.0000 0 180 750 750 2700 Ntspd = ");
    itoa (fs->data_array.Nspd, temp, 10);
    fprintf (fp, "%s\\001\n", temp);
    /*Show the data array organization */
    print_headers ("Data" /* data or tag */, 3750 /* y co-ordinate */, LEFT_x, RIGHT_x, fp, fs);
    print_array ("Data" /* data or tag */, 4050 /* top y */, 9000 /* bottom y */,
         LEFT_x, RIGHT_x, fp, fs, &figdim);
    cy = figdim.bottom_y;
    cy += 700;
    print_headers ("Tag", cy,  LEFT_x, RIGHT_x, fp, fs);
    cy += 400;
    print_array ("Tag", cy, 13000,  LEFT_x, RIGHT_x, fp, fs, &figdim);

    if (fclose (fp) != 0) fprintf (stderr, "Error creating cache.fig file\n");

#ifdef __linux__
    strcpy(cmd, "convert ");
    strcat(cmd, filename);
    strcat(cmd, " ");
    strcat(cmd, fs->file_n);
    system(cmd);
    strcpy(cmd, "mv ");
    strcat(cmd, fs->file_n);
    strcat(cmd, " webapps/cacti/docroot/");
#ifdef WEB_CACTI
    system(cmd);
#endif
    /* remove the .fig file */
    strcpy(cmd, "rm ");
    strcat(cmd, filename);
    system(cmd);
#endif
}

/* insert a mat to .fig file at the specified coordinate*/
void
insert_mat (int x, int y, int l, int w, FILE *file)
{
    char temp1[10], temp2[10];
    int xt = x;
    int yt = y;
    int i=0;
    fprintf (file, "%s", "\t");
    while (i != 4) {
        itoa (x, temp1, 10);
        itoa (y, temp2, 10);
        fprintf (file, "%s ", temp1);
        fprintf (file, "%s ", temp2);
        if (i < 2) {
            if (i%2) y = y + w;
            else x = x + l;
        }
        else {
            if (i%2) y = y - w;
            else x = x - l;
        }
        i++;
    }
    /* finish the loop */
    itoa (xt, temp1, 10);
    itoa (yt, temp2, 10);
    fprintf (file, "%s ", temp1);
    fprintf (file, "%s\n", temp2);
}

/* insert a line at the given co-ordinates */
void
insert_line (int x1, int y1, int x2, int y2, FILE *file)
{
    char temp1[10], temp2[10], temp3[10], temp4[10];
    fprintf (file, "%s ", "\t");
    itoa (x1, temp1, 10);
    itoa (y1, temp2, 10);
    itoa (x2, temp3, 10);
    itoa (y2, temp4, 10);
    fprintf (file, "%s %s %s %s\n", temp1, temp2, temp3, temp4);
}

/* Show the organization of a sample mat */
void
print_mat (FILE *fp)
{
    fprintf (fp, "%s\n", "2 2 0 1 0 17 50 -1 20 0.000 0 0 -1 0 0 5");
    fprintf (fp, "\t%s\n", "7725 1875 8250 1875 8250 2400 7725 2400 7725 1875");
    fprintf (fp, "%s\n", "2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5");
    fprintf (fp, "\t%s\n", "6525 825 8475 825 8475 2625 6525 2625 6525 825");
    fprintf (fp, "%s\n", "2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2");
    fprintf (fp, "\t%s\n", "7575 1725 9000 1725");
    fprintf (fp, "%s\n", "2 4 0 1 0 7 50 -1 -1 0.000 0 0 7 0 0 5");
    fprintf (fp, "\t%s\n", "7650 1800 7650 1650 7425 1650 7425 1800 7650 1800");
    fprintf (fp, "%s\n", "2 1 1 1 0 7 50 -1 -1 4.000 0 0 -1 1 0 2");
    fprintf (fp, "\t%s\n", "2 1 2.00 90.00 135.00");
    fprintf (fp, "\t%s\n", "8925 2850 7650 1800");
    fprintf (fp, "%s\n", "2 2 0 1 0 17 50 -1 20 0.000 0 0 -1 0 0 5");
    fprintf (fp, "\t%s\n", "7725 1050 8250 1050 8250 1575 7725 1575 7725 1050");
    fprintf (fp, "%s\n", "2 2 0 1 0 17 50 -1 20 0.000 0 0 -1 0 0 5");
    fprintf (fp, "\t%s\n", "6825 1050 7350 1050 7350 1575 6825 1575 6825 1050");
    fprintf (fp, "%s\n", "2 2 0 1 0 17 50 -1 20 0.000 0 0 -1 0 0 5");
    fprintf (fp, "\t%s\n", "6825 1875 7350 1875 7350 2400 6825 2400 6825 1875");
    fprintf (fp, "%s\n", "2 1 1 1 0 7 50 -1 -1 4.000 0 0 -1 1 0 2");
    fprintf (fp, "\t%s\n", "2 1 3.00 90.00 135.00");
    fprintf (fp, "\t%s\n", "9000 675 8250 1050");
    fprintf (fp, "%s\n", "4 0 0 50 -1 1 14 0.0000 0 195 870 9150 675 Sub-array\\001");
    fprintf (fp, "%s\n", "4 0 0 50 -1 1 14 0.0000 0 150 360 7350 2925 Mat\\001");
    fprintf (fp, "%s\n", "4 0 0 50 -1 1 14 0.0000 0 150 1695 8250 3150 Central "
                         "Predecoder\\001");
}

/* Print data/tag array parameteres */
void
print_headers (char *c, int y, int lx, int rx, FILE *fp, final_results *fs)
{
    char *color1, *color2;
    char temp1[10], temp2[10], temp3[10];
    int x = lx;
    results_mem_array *rma;
        if (strcmp (c, "Data") == 0) {
        rma = &fs->data_array;
        color1 = "8";
        color2 = "1";
    }
    else {
        rma = &fs->tag_array;
        color1 = "18";
        color2 = "4";
    }

    fprintf (fp, "%s %s %s", "4 0", color1, "50 -1 1 20 0.0000 4 255 2925 ");
    itoa (y, temp2, 10);
    lx += 100;
    itoa (lx, temp1, 10);
    fprintf (fp, "%s %s %s %s\n", temp1, temp2, c, "Array Organization\\001");
    lx += 3100;
    itoa (lx, temp1, 10);
    fprintf (fp, "%s %s %s %s %s %s", "4 0", color2, "50 -1 1 14 0.0000 4 180 "
        "3240", temp1, temp2, "(Mat height - ");

    itoa ((int)rma->mat_height, temp1, 10);
    fprintf (fp, "%su) %s", temp1, " (Mat width - ");
    itoa ((int)rma->mat_width, temp1, 10);
    fprintf (fp, "%su) (%s %s", temp1, c, "array dim. - ");
    itoa ((int)rma->bank_height, temp1, 10);
    fprintf (fp, "%su %s", temp1, "x ");
    itoa ((int)rma->bank_width, temp1, 10);
    fprintf (fp, "%su)%s\n", temp1, "\\001");
    fprintf (fp, "%s\n", "2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2");
    itoa (x, temp1, 10);
    itoa (y+150, temp2, 10);
    itoa (rx, temp3, 10);
    fprintf (fp, "\t%s %s %s %s\n", temp1, temp2, temp3, temp2);
}

void
print_array (char *c, int ty, int by, int lx, int rx,
             FILE *fp, final_results *fs, struct fig_mat_dim *matdim)
{
    int max_box = 1600;
    int min_box = 700;
    int spacing = 0;
    int h_exceed = 0;
    int v_exceed = 0;
    int center, h_range, v_range, cur_x, cur_y;
    results_mem_array *rma;
    int wl, bl;
    int h, i, j, l, w;
    double mat_h, mat_w;

    if (strcmp (c, "Data") == 0) {
        rma = &fs->data_array;
        bl = rma->Ndbl/2;
        wl = rma->Ndwl/2;
        center = (rx + 200 - lx)/2;
        h_range = (rx - lx);
        v_range = (by - ty);

        /* calculate l and w of the mat */
        w = (h_range - (wl * spacing))/wl;
        if (w > max_box) {
            w = max_box;
        }
        else if (w < min_box) {
            w = min_box;
            h_exceed = 1;
        }
        l = (v_range - (bl * spacing))/bl;
        if (l > max_box) {
            l = max_box;
        }
        else if (l < min_box) {
            l = min_box;
            v_exceed = 1;
        }
        mat_h = rma->mat_height;
        mat_w = rma->mat_width;
        if (mat_h > mat_w) {
            l = MIN (l, w);
            w = (int)((double)l * (mat_w/mat_h));
        }
        else {
            w = MIN (l, w);
            l = (int)((double)w * (mat_h/mat_w));
        }
    }
    else {
        rma = &fs->tag_array;
        bl = rma->Ndbl/2;
        wl = rma->Ndwl/2;
        center = (rx + 200 - lx)/2;
        h_range = (rx - lx);
        v_range = (by - ty);

        /* calculate l and w of the mat */
        l = ((double) rma->mat_height * ((double) matdim->h/fs->data_array.mat_height));
        w = ((double) rma->mat_width * ((double) matdim->w/fs->data_array.mat_width));
    }

    /* draw the array */
    cur_x = center+spacing/2;
    cur_y = ty;
    j = 0;
    /* regular cases */
    /* draw all the mats */
    for (h=0; h<bl; h++) {
        for (i=0; i<wl; i++) {
            if (i == wl/2) {
                j=1;
                cur_x = center - (w + spacing/2);
            }
            /* if the mat size is big then use thicker lines */
            if (l > 600) {
                fprintf (fp, "%s\n", "2 2 0 3 0 7 50 -1 -1 0.000 0 0 -1 0 0 5");
            }
            else {
                fprintf (fp, "%s\n", "2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5");
            }
            insert_mat (cur_x, cur_y, w, l, fp);
            if (j) {
                cur_x -= (spacing + w);
            }
            else {
                cur_x += (spacing + w);
            }
        }
        cur_y += (spacing + l);
        cur_x = center+spacing/2;
        j = 0;
    }
    matdim->h = l;
    matdim->w = w;
    matdim->top_y = ty;
    matdim->bottom_y = cur_y;
}

void
output_NUCA (nuca_org_t *fr)
{
    fprintf(stderr, "\n---------- CACTI version 6.0, Non-uniform Cache Access "
            "----------\n\n");
    fprintf(stderr, "Optimal number of banks - %d\n", fr->bank_count);
    fprintf(stderr, "Grid organization rows x columns - %d x %d\n", 
            fr->rows, fr->columns);
    fprintf(stderr, "Average access latency to a random bank \n\t"
        "(Bank Access time + Avg. Network Delay + Contention Cycles)- %g cycles\n", 
            fr->nuca_pda.delay);
    fprintf(stderr, "Average dynamic energy/access (nJ) - %g \n",
            fr->nuca_pda.power.dynamic*1e9);
    fprintf(stderr, "Network frequency - %g GHz\n",
            (1/fr->nuca_pda.cycle_time)*1e3);
    fprintf(stderr, "Cache dimension (mm x mm) - %g x %g\n",
            fr->nuca_pda.area_stats.height, 
            fr->nuca_pda.area_stats.width);
//    fprintf(stderr, "\n\nBank stats:\n");
//    fprintf(stderr, "\tBank size (bytes) - %d\n", (int)fr->params->cache_size/fr->bank_count);
//    fprintf(stderr, "\tAccess time (ns)- %g\n", fr->bank_pda.delay, FREQUENCY);
//    fprintf(stderr, "\tDynamic Energy - %g J/access\n", fr->bank_pda.power.dynamic);
//    fprintf(stderr, "\tLeakage power - %g W\n", fr->bank_pda.power.leakage);
//    fprintf(stderr, "\tHeight (mm) x Width (mm) -   %g x %g\n", 
//            fr->bank_pda.area_stats.height, fr->bank_pda.area_stats.width);

    print_router(&(fr->router));

    fprintf(stderr, "\n\nWire stats:\n");
    if (fr->h_wire.wt == Global) {
        fprintf(stderr, "\tWire type - Full swing global wires with least "
                "possible delay\n");
    }
    else if (fr->h_wire.wt == Global_10) {
        fprintf(stderr, "\tWire type - Full swing global wires with "
                "10%% delay penalty\n");
    }
    else if (fr->h_wire.wt == Global_20) {
        fprintf(stderr, "\tWire type - Full swing global wires with "
                "20%% delay penalty\n");
    }
    else if (fr->h_wire.wt == Global_30) {
        fprintf(stderr, "\tWire type - Full swing global wires with "
                "30%% delay penalty\n");
    }
    else if(fr->h_wire.wt == Low_swing) {
        fprintf(stderr, "\tWire type - Low swing wires\n");
    }

    //    fprintf(stderr, "\tWire width - %g nm\n", 
    //                    fr->h_wire.wire_pda.area_stats.width*1e9);
    fprintf(stderr, "\tWire width - %g nm\n", 
            fr->h_wire.wire_width*1e9);
    fprintf(stderr, "\tWire spacing - %g nm\n", 
            fr->h_wire.wire_spacing*1e9);
    fprintf(stderr, "\tHorizontal link delay - %g (ns)\n", 
            fr->h_wire.wire_pda.delay*1e9);
    fprintf(stderr, "\tVertical link delay - %g (ns)\n", 
            fr->v_wire.wire_pda.delay*1e9);
    fprintf(stderr, "\tDelay/length - %g (ns/mm)\n", 
            fr->h_wire.wire_pda.delay*1e9/fr->bank_pda.area_stats.width);
//    fprintf(stderr, "\tDelay/length - %g (ns/mm)\n", 
//            fr->v_wire.wire_pda.delay*1e9/fr->bank_pda.area_stats.height);
    if (fr->h_wire.wt != Low_swing) { //FIXME
//        fprintf(stderr, "\tRepeater size - %g times the smallest value\n", 
//                fr->h_wire.repeater_size);
//        fprintf(stderr, "\tRepeater spacing - %g u\n\n", 
//                fr->h_wire.repeater_spacing*1000000);
    }
    fprintf(stderr, "\tHorizontal link energy -dynamic/access %g (nJ)\n"
            "\t                       -leakage %g (nW)\n\n", 
            fr->h_wire.wire_pda.power.dynamic*1e9,
            fr->h_wire.wire_pda.power.leakage*1e9);
    fprintf(stderr, "\tVertical link energy -dynamic/access %g (nJ)\n"
            "\t                     -leakage %g (nW)\n\n", 
            fr->v_wire.wire_pda.power.dynamic*1e9,
            fr->v_wire.wire_pda.power.leakage*1e9);
    fprintf(stderr, "\tEnergy/length per wire - %g (nJ/mm)\n", 
            fr->h_wire.wire_pda.power.dynamic*1e9/(fr->h_wire.wire_length*1e3));
//    fprintf(stderr, "\tEnergy/length - %g (ns/mm)\n", 
//            fr->v_wire.wire_pda.power.dynamic*1e9/fr->bank_pda.area_stats.height);
}

    void
output_UCA (uca_org_t *fr)
{
    if (NUCA) {
        fprintf(stderr, "\n\n Detailed Bank Stats:\n" );
    fprintf(stderr,"    Bank Size (bytes): %d\n",
            (int) (fr->params->cache_size));
    }
    else {
        if (fr->params->dram) {
            fprintf(stderr, "\n---------- CACTI version 6.0, Uniform Cache Access "
                    "%s ----------\n", "DRAM Model");
        }
        else {
            fprintf(stderr, "\n---------- CACTI version 6.0, Uniform Cache Access "
                    "%s ----------\n", "SRAM Model");
        }
        fprintf(stderr,"\nCache Parameters:\n");
        fprintf(stderr,"    Total cache size (bytes): %d\n",
            (int) (fr->params->cache_size));
    }

    fprintf(stderr,"    Number of banks: %d\n",
            (int)fr->params->uca_banks);
    if (fr->params->fully_assoc)
        fprintf(stderr,"    Associativity: fully associative\n");
    else {
        if (fr->params->tag_associativity==1)
            fprintf(stderr,"    Associativity: direct mapped\n");
        else
            fprintf(stderr,"    Associativity: %d\n",
                    fr->params->tag_associativity);
    }
    fprintf(stderr,"    Block size (bytes): %d\n",fr->params->block_size);
    fprintf(stderr,"    Read/write Ports: %d\n",
            fr->params->rw_ports);
    fprintf(stderr,"    Read ports: %d\n",
            fr->params->excl_read_ports);
    fprintf(stderr,"    Write ports: %d\n",
            fr->params->excl_write_ports);
    fprintf(stderr,"    Technology size: %2.2fum\n", 
            fr->params->tech_size);
    //    fprintf(stderr,"    Vdd: %2.1fV\n", fr->params->vdd_periph_global);

    fprintf(stderr,"\n    Access time (ns): %g\n",fr->access_time*1e9);
    fprintf(stderr,"    Cycle time (ns):  %g\n",fr->cycle_time*1e9);
    fprintf(stderr,"    Total dynamic read energy per access (nJ):%g\n",
            fr->power.readOp.dynamic*1e9);
    fprintf(stderr,"    Total leakage power of a bank"
            " (mW):%g\n", fr->power.readOp.leakage*1e3);

    if (fr->params->dram) {
        fprintf(stderr, "    Refresh power (mW): %g\n", 
                fr->data_array->refresh_power*1e3); //FIXME
    }

    fprintf(stderr, "    Cache height x width (mm):"
            " %f x %f\n", fr->cache_ht, fr->cache_len);


    fprintf(stderr,"\n    Best Ndwl (L1): %d\n",fr->data_array->Ndwl);
    fprintf(stderr,"    Best Ndbl (L1): %d\n",fr->data_array->Ndbl);
    fprintf(stderr,"    Best Nspd (L1): %f\n",fr->data_array->Nspd);
    fprintf(stderr,"    Best Ndcm (L1): %d\n",fr->data_array->Ndsam);
    fprintf(stderr,"    Best Ndsam (L1): %d\n\n",fr->data_array->deg_bitline_muxing);

    if ((!fr->params->pure_sram)) {
        fprintf(stderr,"    Best Ntwl (L1): %d\n",fr->tag_array->Ndwl);
        fprintf(stderr,"    Best Ntbl (L1): %d\n",fr->tag_array->Ndbl);
        fprintf(stderr,"    Best Ntspd (L1): %f\n",fr->tag_array->Nspd);
        fprintf(stderr,"    Best Ntcm (L1): %d\n",fr->tag_array->Ndsam);
        fprintf(stderr,"    Best Ntsam (L1): %d\n\n",fr->tag_array->deg_bitline_muxing);
    }

    switch (fr->data_array->wt) {
        case (0):
            fprintf(stderr, "    Data array, H-tree wire type: %s\n", "Delay "
                    "optimized global wires");
            break;
        case (1):
            fprintf(stderr, "    Data array, H-tree wire type: %s\n", "Global "
                    "wires with 10\% delay penalty");
            break;
        case (2):
            fprintf(stderr, "    Data array, H-tree wire type: %s\n", "Global "
                    "wires with 20\% delay penalty");
            break;
        case (3):
            fprintf(stderr, "    Data array, H-tree wire type: %s\n", "Global "
                    "wires with 30\% delay penalty");
            break;
        case (4):
            fprintf(stderr, "    Data array, wire type: %s\n", "Low "
                    "swing wires");
            break;
        default:
            printf("ERROR - Unknown wire type %d!\n", (int) fr->data_array->wt);
            exit(-1);
    }

    if ((!fr->params->pure_sram)) {
        switch (fr->tag_array->wt) {
            case (0):
                fprintf(stderr, "    Tag array, H-tree wire type: %s\n", "Delay "
                        "optimized global wires");
                break;
            case (1):
                fprintf(stderr, "    Tag array, H-tree wire type: %s\n", "Global "
                        "wires with 10\% delay penalty");
                break;
            case (2):
                fprintf(stderr, "    Tag array, H-tree wire type: %s\n", "Global "
                        "wires with 20\% delay penalty");
                break;
            case (3):
                fprintf(stderr, "    Tag array, H-tree wire type: %s\n", "Global "
                        "wires with 30\% delay penalty");
                break;
            case (4):
                fprintf(stderr, "    Tag array, wire type: %s\n", "Low "
                        "swing wires");
                break;
            default:
                printf("ERROR - Unknown wire type %d!\n", (int) fr->tag_array->wt);
                exit(-1);
        }
    }

    if (fr->params->print_detail) {

        if(fr->params->fully_assoc) return;

        /* Delay stats */
        /* data array stats */ 
        fprintf(stderr,"\nTime Components:\n\n");

        fprintf(stderr,"  Data side (with Output driver) (ns): %g\n",
                fr->data_array->access_time/1e-9);

        fprintf(stderr, "\tH-tree input delay (ns): %g\n", 
                fr->data_array->delay_route_to_bank * 1e9 +
                fr->data_array->delay_addr_din_horizontal_htree * 1e9 +
                fr->data_array->delay_addr_din_vertical_htree * 1e9);

        fprintf(stderr, "\tDecoder + wordline delay (ns): %g\n",
                fr->data_array->delay_row_predecode_driver_and_block * 1e9 +
                fr->data_array->delay_row_decoder * 1e9);

        fprintf(stderr, "\tBitline delay (ns): %g\n",
                fr->data_array->delay_bitlines/1e-9);

        fprintf(stderr, "\tSense Amplifier delay (ns): %g\n",
                fr->data_array->delay_sense_amp * 1e9);


        fprintf(stderr, "\tH-tree output delay (ns): %g\n",
                fr->data_array->delay_subarray_output_driver * 1e9 +
                fr->data_array->delay_dout_vertical_htree * 1e9 +
                fr->data_array->delay_dout_horizontal_htree * 1e9);

        //    if (fr->params->print_detail) {
        //        fprintf(stderr, "    Subarray output driver (ns): %g\n",
        //                fr->data_array->delay_subarray_output_driver * 1e9);
        //        fprintf(stderr, "    Bitline mux predecoder and decoder (ns): %g\n",
        //                fr->data_array->delay_senseamp_mux_predecode_driver_and_block 
        //                * 1e9);
        //        fprintf(stderr, "    Sense amp mux decoder (ns): %g\n", 
        //                fr->data_array->delay_senseamp_mux_decoder * 1e9);
        //    }
        //
        if ((!fr->params->pure_sram)) {
            /* tag array stats */
            fprintf(stderr,"\n  Tag side (with Output driver) (ns): %g\n",
                    fr->tag_array->access_time/1e-9);

            fprintf(stderr, "\tH-tree input delay (ns): %g\n", 
                    fr->tag_array->delay_route_to_bank * 1e9 +
                    fr->tag_array->delay_addr_din_horizontal_htree * 1e9 +
                    fr->tag_array->delay_addr_din_vertical_htree * 1e9);

            fprintf(stderr, "\tDecoder + wordline delay (ns): %g\n",
                    fr->tag_array->delay_row_predecode_driver_and_block * 1e9 +
                    fr->tag_array->delay_row_decoder * 1e9);

            fprintf(stderr, "\tBitline delay (ns): %g\n",
                    fr->tag_array->delay_bitlines/1e-9);

            fprintf(stderr, "\tSense Amplifier delay (ns): %g\n",
                    fr->tag_array->delay_sense_amp * 1e9);

            fprintf(stderr, "\tComparator delay (ns): %g\n",
                    fr->tag_array->delay_comparator * 1e9);

            fprintf(stderr, "\tH-tree output delay (ns): %g\n",
                    fr->tag_array->delay_subarray_output_driver * 1e9 +
                    fr->tag_array->delay_dout_vertical_htree * 1e9 +
                    fr->tag_array->delay_dout_horizontal_htree * 1e9);
        }

        /* Energy/Power stats */
        fprintf(stderr,"\nPower Components:\n");
        fprintf(stderr,"\n  Data array: Total dynamic read energy/access  (nJ): %g\n",
                fr->data_array->total_power.readOp.dynamic * 1e9);
        fprintf(stderr,"\tTotal leakage read/write power all banks at "
                "maximum frequency (mW): %g\n",
                fr->data_array->total_power.readOp.leakage * 1e3);
        fprintf(stderr,"\tTotal energy in H-tree (that includes both "
                "address and data transfer) (nJ): %g\n", 
                (fr->data_array->power_addr_vertical_htree.readOp.dynamic +
                 fr->data_array->power_datain_vertical_htree.readOp.dynamic +
                 fr->data_array->power_routing_to_bank.readOp.dynamic) * 1e9);
        fprintf(stderr, "\tDecoder (nJ): %g\n",
                fr->data_array->power_row_predecoder_drivers.readOp.dynamic * 1e9 +
                fr->data_array->power_row_predecoder_blocks.readOp.dynamic * 1e9);
        fprintf(stderr, "\tWordline (nJ): %g\n",
                fr->data_array->power_row_decoders.readOp.dynamic * 1e9);
        fprintf(stderr, "\tBitline mux & associated drivers (nJ): %g\n",
                fr->data_array->power_bit_mux_predecoder_drivers.readOp.dynamic * 1e9 +
                fr->data_array->power_bit_mux_predecoder_blocks.readOp.dynamic * 1e9 +
                fr->data_array->power_bit_mux_decoders.readOp.dynamic * 1e9);
        fprintf(stderr, "\tSense amp mux & associated drivers (nJ): %g\n",
                fr->data_array->power_senseamp_mux_predecoder_drivers.readOp.dynamic * 1e9 +
                fr->data_array->power_senseamp_mux_predecoder_blocks.readOp.dynamic * 1e9 +
                fr->data_array->power_senseamp_mux_decoders.readOp.dynamic * 1e9);
        fprintf(stderr, "\tBitlines (nJ): %g\n",
                fr->data_array->power_bitlines.readOp.dynamic * 1e9);
        fprintf(stderr, "\tSense amplifier energy (nJ): %g\n",
                fr->data_array->power_sense_amps.readOp.dynamic * 1e9);
        fprintf(stderr, "\tSub-array output driver (nJ): %g\n",
                fr->data_array->power_output_drivers_at_subarray.readOp.dynamic * 1e9);
        if ((!fr->params->pure_sram)) {
            fprintf(stderr, "\n  Tag array:  Total dynamic read energy/access (nJ): %g\n",
                    fr->tag_array->total_power.readOp.dynamic * 1e9);
            fprintf(stderr,"\tTotal leakage read/write power all banks at "
                    "maximum frequency (mW): %g\n",
                    fr->tag_array->total_power.readOp.leakage * 1e3);
            fprintf(stderr,"\tTotal energy in H-tree "
                    " (nJ): %g\n", 
                    fr->tag_array->power_addr_vertical_htree.readOp.dynamic * 1e9 +
                    fr->tag_array->power_datain_vertical_htree.readOp.dynamic * 1e9 +
                    fr->tag_array->power_routing_to_bank.readOp.dynamic * 1e9);
            fprintf(stderr, "\tDecoder (nJ): %g\n",
                    fr->tag_array->power_row_predecoder_drivers.readOp.dynamic * 1e9 +
                    fr->tag_array->power_row_predecoder_blocks.readOp.dynamic * 1e9);
            fprintf(stderr, "\tWordline (nJ): %g\n",
                    fr->tag_array->power_row_decoders.readOp.dynamic * 1e9);
            fprintf(stderr, "\tBitline mux & associated drivers (nJ): %g\n",
                    fr->tag_array->power_bit_mux_predecoder_drivers.readOp.dynamic * 1e9 +
                    fr->tag_array->power_bit_mux_predecoder_blocks.readOp.dynamic * 1e9 +
                    fr->tag_array->power_bit_mux_decoders.readOp.dynamic * 1e9);
            fprintf(stderr, "\tSense amp mux & associated drivers (nJ): %g\n",
                    fr->tag_array->power_senseamp_mux_predecoder_drivers.readOp.dynamic * 1e9 +
                    fr->tag_array->power_senseamp_mux_predecoder_blocks.readOp.dynamic * 1e9 +
                    fr->tag_array->power_senseamp_mux_decoders.readOp.dynamic * 1e9);
            fprintf(stderr, "\tBitlines (nJ): %g\n",
                    fr->tag_array->power_bitlines.readOp.dynamic * 1e9);
            fprintf(stderr, "\tSense amplifier energy (nJ): %g\n",
                    fr->tag_array->power_sense_amps.readOp.dynamic * 1e9);
            fprintf(stderr, "\tSub-array output driver (nJ): %g\n",
                    fr->tag_array->power_output_drivers_at_subarray.readOp.dynamic * 1e9);
        }
        fprintf(stderr, "\nArea Components:\n");
        /* Data array area stats */
        fprintf(stderr, "\n  Data array: Area (mm2): %g\n", fr->data_array->area);
        fprintf(stderr, "\tHeight (mm): %g\n",
                fr->data_array->all_banks_height*1e-3);
        fprintf(stderr, "\tWidth (mm): %g\n",
                fr->data_array->all_banks_width*1e-3);
        if (fr->params->print_detail) {
            fprintf(stderr, "    Area efficiency (Memory cell area/Total area) - %g\n",
                    fr->data_array->area_efficiency);
        }
        /* Tag array area stats */
        if ((!fr->params->pure_sram)) {
            fprintf(stderr, "\n  Tag array: Area (mm2): %g\n", fr->tag_array->area);
            fprintf(stderr, "\tHeight (mm): %g\n",
                    fr->tag_array->all_banks_height*1e-3);
            fprintf(stderr, "\tWidth (mm): %g\n",
                    fr->tag_array->all_banks_width*1e-3);
            if (fr->params->print_detail) {
                fprintf(stderr, "    Area efficiency (Memory cell area/Total area) - %g\n",
                        fr->tag_array->area_efficiency);
            }
        }
    }
}
