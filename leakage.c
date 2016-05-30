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

double FO4;

double cmos_ileakage(double nWidth, double pWidth, int temp) 
{
	double leakage = 0.0;
	static double norm_nmos_leakage = 0;
	static double norm_pmos_leakage = 0;

	if((!is_dram)&&(is_sram_cell)){//SRAM cell access transistor
		norm_nmos_leakage = I_off_n_sram_cell_transistor[temp-300];
		norm_pmos_leakage = I_off_p_sram_cell_transistor[temp-300];
	}
	else
		if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
			norm_nmos_leakage = I_off_n_sram_cell_transistor[temp-300];
			norm_pmos_leakage = I_off_p_sram_cell_transistor[temp-300];
		}
		else
			if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
				norm_nmos_leakage = I_off_n_dram_wordline_transistor[temp-300];
				norm_pmos_leakage = I_off_p_dram_wordline_transistor[temp-300];
			}
			else{//DRAM or SRAM all other transistors
				norm_nmos_leakage = I_off_n_periph_global[temp-300];
				norm_pmos_leakage = I_off_p_periph_global[temp-300];
			}
			leakage = nWidth * norm_nmos_leakage + pWidth * norm_pmos_leakage;
	return leakage;
}


double simplified_nmos_leakage(double nwidth, int temp)
{
	double nIleak;
	if((!is_dram)&&(is_access_transistor)){//SRAM cell access transistor
		nIleak = nwidth * I_off_n_sram_cell_transistor[temp-300];
	}
	else
		if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
			nIleak = nwidth * I_off_n_sram_cell_transistor[temp-300];
		}
		else
			if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
				nIleak = nwidth * I_off_n_dram_wordline_transistor[temp-300];	    
			}
			else{//DRAM or SRAM all other transistors
				nIleak = nwidth * I_off_n_periph_global[temp-300];	
			}
	return nIleak;
}


double simplified_pmos_leakage(double pwidth, int temp)
{
	double pIleak;
	
	if((!is_dram)&&(is_wordline_transistor)){//SRAM wordline transistor
			pIleak = pwidth * I_off_p_sram_cell_transistor[temp-300];
		}
	else
		if((is_dram)&&(is_wordline_transistor)){//DRAM wordline transistor
			pIleak = pwidth * I_off_p_dram_wordline_transistor[temp-300];	    
		}
		else{//DRAM or SRAM all other transistors
			pIleak = pwidth * I_off_p_periph_global[temp-300];	
		}
		return pIleak;
}



  

   
  



