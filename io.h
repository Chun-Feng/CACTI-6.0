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

#ifndef _ioh
#define _ioh

#include "time.h"
#include "router.h"

typedef struct fig_mat_dim {
    int h;
    int w;
    int top_y;
    int bottom_y;
}fig_mat_dim_t;

/* Utility functions */
void strreverse(char *, char*);
#ifdef __linux__
void itoa(int, char*, int);
#endif
void read_file();
int get_cpucount();

/* output functions */
void dump_input_args(input_params_t *b);
void output_UCA (uca_org_t *fr);
void output_NUCA (nuca_org_t *fr);

/* output figure */
void insert_mat (int x, int y, int l, int w, FILE *file);
void insert_line (int x1, int y1, int x2, int y2, FILE *file);
void print_mat (FILE *fp);
void print_headers (char *c, int y, int lx, 
    int rx, FILE *fp, final_results *fs);
void print_array (char *c, int ty, int by, int lx, 
        int rx, FILE *fp, final_results *fs, 
        struct fig_mat_dim *matdim);
void fig_out (final_results *fr);

void validate_cache_args (input_params_t *params);
int parse_cmd_args(int argc,char *argv[], input_params_t *params);
void parse_cfg(input_params_t *params);
void sim_cache(int argc,char *argv[]);

#endif

//void * sim_nucaold(bank_out_t *);
//bank_out_t *
//find_optimal_config(bank_out_t **, int);

//void dump_bank_out(bank_out_t *);
