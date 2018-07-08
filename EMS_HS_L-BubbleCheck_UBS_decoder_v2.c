/* EMS_HS_L-BubbleCheck_UBS_decoder_v2.c

  Extended Min-Sum Decoder
  Horizontal Scheduling
  Variable Node modified (to take into account all LLR values, instead of offset for some GF symbols)
  L-Bubble-Check 2rows, 2 columns, with optimal implementation to speed-up simulation (no IFs)
*/

/* 11/02/2009 L. Conde-Canencia  (UBS)
/* L-bubble check: simplified Bubble check algorithm
/* only the values of the first 2 columns and the first 2 rows of the tab_aux matrix are considered
/* simulation shows that there is no performance loss, and complexity is significantly reduced


/* 10/02/2009 Version by A. Al-Ghouwayel and L. Conde-Canencia  (UBS)
/* Decoder Optimized at Variable Node level
/* in ModelChannel, the LLR values are calculated for code.GF values
/* in DataPassLog (VN processing) the CtoV messages are combined with the LLR values (and not with the offset values for some of them)
/* This optimization is already considered in the Fixed Point model


/* 17/10/2008  Version by A. Al-Ghouwayel and L. Conde-Canencia  (UBS)

/* elementary check processing implemented with the bubble check algorithm (see D6.1.2)
/* check processing in forward-backward (total sum not implemented)
/*


/** ***********************
 * ====>>>>  IMPORTANT !!!!! IMPORTANT !!!!! IMPORTANT !!!!!  <<<<====
 * Some Window platforms will not compile this code, because the uniform random number generator function "drand48()"
 * is NOT a standard C-ANSI function.
 * In this case, you have to replace  the drand48() function with an equivalent function on your platform.
 * Example for Visual C and Dev C, you can use the following general function to replace drand48():
 *
 * // random generator: uniform distribution on [0,1]
 * float My_drand48 (void)
 * {
 *   return((float)(rand())/(float)(RAND_MAX+1.0));
 * }
 *
 * Please note that, according to NUMERCIAL RECIPES:
 * the built-in generators srand() and rand() have no standard implementation and are often badly flawed.
 *
 * A bad random generator can lead to some performance degradation.
 */


/*
 * Module Name       	: EMS_HS_bubbleCheck_ForBack_byUBS.c
 * Short Description 	: Non-binary LDPC reduced complexity decoder with horizontal scheduling and bubble check
 * Parameters
 *
 * Inputs
 *
 *		NbMonteCarlo     : # simulated frames
 *		NbIterMax        : # of maximum decoding iteration
 *		FileMatrix       : File name of the parity-check matrix
 *		EbN              : Eb/No (dB)
 *		NbMax            : size of truncated messages
 *		Offset           : offset correction factor (0.4 -- 1)
 *		NbOper           : Maximum number of operations for sorting
 * Output
 *              Frame Error Rate for the given simulation parameters
 * Input File : 'FileMatrix' is an ASCII file with the parity-check matrix description in aList format.
 * Output File : non
 *
 * Test description: Run the executable with the following parameters
	./EMS_HS_bubbleCheck_ForBack.exe 10000 100 Mat24_N192_Alist 1.75 24 1 48 result

	 NbMonteCarlo     : 10000
	 NbIterMax        : 100
	 FileMatrix       : Mat24_N192_Alist
	 Eb/No (dB)       : 1.75
	 NbMax            : 24
	 Offset           : 1
	 NbOper           : 48
LDPC code parameters: N=192  K=96  M=96  CodeRate=0.5  GF=64  logGF=6
Loading of the binary image of GF(64): Success
<0> FER= 13/ 10000 = 0.0013013

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <err.h>

#define SQR(A) ((A)*(A))   //calcule le carré d'un nombre
#define BPSK(x) (1-2*(x))  // modulation BPSK
#define PI 3.1415926536
#define STR_MAXSIZE 350    //la taille maximale permise pour les chaines de caractères. Utilisée dans la fonction calloc().

/******************************************************/
/****************** LDPC code description**************/
/******************************************************/
// The parity-check matrix has a sparse description.
// ce type contient les paramètres du code LDPC considéré. Ces paramètres sont collectés depuis un fichier texte contenant la matrice
// de parité du code.

typedef struct {
	int N;			/* number of columns in H */
	int M;			/* number of rows in H */
	int K;			/* number of information symbols : K = N-M */
	int GF;			/* Field order (eg. GF=256) */
	int logGF;		/* logGF = log2(GF)  ( given GF = 2^q => logGF = q) : logGF is the number of bits forming GF symbols */
	int **mat;		/* the Parity-check matrix : tableau bidimensionnel contenant les VN (i.e les colonnes) qui participe dans chaque contrainte
					   de parité  (i.e les lignes) */
	int **matValue;		/* Parity-check matrix non-binary coefficients : contient les coefficient GF(q) pour chaque ligne */
	int *rowDegree;		/* rowDegree[i] = the i^th check node degree */
	int *columnDegree; 	/* columnDegree[j] = the j^th variable node degree */
	int *interleaver;	/* Interleaver: label the ones by reading row-wise the parity-check matrix H,
					    and read the labels by reading column-wise H */
	int nbBranch;		/* Size of the interleaver = number of edges in the Tanner graph */
	float rate; 		/* Code rate */
    int **NtoB;         /* N c'est Node (variable node) et B c'est Branche (edge)
						  NtoB est le tableau bidimentionnelle ColonneBranche[][] de l'exemple précedant */
	int **matUT;		/* Upper Triangular form of the parity-check matrix after Gaussian elimination. matUT is used for encoding*/
	int *Perm;		/* Permutation used for Gaussian Elimination  */
} code_t;

/**********************************/
/*** Computation tables in GF(q)***/
/**********************************/
typedef struct {
    int **BINGF;		/* Mapping symbol GFq -> ensemble de symboles binaires */
    int **ADDGF;		/* Addition table in GFq */
    int **MULGF;		/* Multiplication table in GFq */
    int **DIVGF;		/* Division table in GFq */
} table_t;

/** *************************************************************************/
/******************************** EMS decoder********************************/
/*** Structure for the messages on the edges on decoding graph Likelihood ***/
/****************************************************************************/
typedef float softdata_t ; /* la fiabilité d'un message est de type float */
typedef struct {
	int N;
	int nbMax; /* c'est le paramètre nm dans l'algo EMS : Au lieu de transmettre les q symbole dans GF(q), on fait une trancature et on transmet
				uniquement les nm plus fiables symboles */
	int nbBranch; /* nombre de branches dans le graphe de Tanner */
	softdata_t 	**CtoV;		/* An array to store the nbMax largest Check_to_Variable reliabilities for all edges on the graph
							CtoV Array size = (nbBranch x nbMax) */
	int 		**CtoVindex; 	/*An array to store the nbMax symbols Indexes corresponding to the nbMax largest reliabilities stored in CtoV
								CtoVindex Array size = the same as for C_To_V */
	softdata_t 	**VtoC;		/* An array to store the nbMax largest Variable_to_Check reliabilities for all edges on the graph
								Same size as for CtoV and CtoVindex*/
	int 		**VtoCindex; 	/*An array to store the nbMax symbols Indexes corresponding to the nbMax largest reliabilities stored in VtoC
								Same size as for V_To_C */
	softdata_t 	**likelihood;	   /* An array to store intrinsic Log Likelihood Ratios received from channel. Each variable node has GF intrinsic
                                      LLRs. So, the size of likelihood is (N x GF) where N is the number of VNs and GF is the order of the Galois
                                      field. The values are sorted in decreasing way */
	int 		**likelihood_index; /*  Galois field symbols associated to the LLRs values in 'likelihood'
										Same size as for 'likelihood' */
	softdata_t 	**APP;		/* Array to store A Posteriori Probability used to make hard decision on symbol value
						       Same size as for 'likelihood'*/
        softdata_t      **VC_arr; //l'entrée d'un check node processor. Je l'aurais appeler CNPInput_LLR
        int             **IndiceVC_arr; //Les indices correspondant à VC_arr Je l'aurais appeler CNPInput_Index
        softdata_t      **CV_arr ; //la sortie d'un check node processor. Je l'aurais appeler CNPOutput_LLR
        int             **IndiceCV_arr; //Les indices correspondant à CV_arr Je l'aurais appeler CNPOutput_Index
        softdata_t      *CV; //L'entrée d'un VN processor. C'est un vecteur car dv=2. Je l'aurais appeler VNPInput_LLR
        int             *IndiceCV; //les indices correpondant à CV. Je l'aurais appeler VNPInput_Index
        softdata_t      *VC; //la sortie d'un VN processor. Un vecteur car dv=2. Je l'aurais appeler VNPOutput_LLR
        int             *IndiceVC; //les indices correpondant à VC. Je l'aurais appeler VNPOutput_Index
} decoder_t;

/*
 * Binary image of the field GF(64)
 * Primitive polynomial used P(x)=X^6+X+1 (GF(64))
 */

static const int BinGF_256[256][8]={0,	0,	0,	0,	0,	0,	0,	0,
 1,	0,	0,	0,	0,	0,	0,	0,
 0,	1,	0,	0,	0,	0,	0,	0,
 0,	0,	1,	0,	0,	0,	0,	0,
 0,	0,	0,	1,	0,	0,	0,	0,
 0,	0,	0,	0,	1,	0,	0,	0,
 0,	0,	0,	0,	0,	1,	0,	0,
 0,	0,	0,	0,	0,	0,	1,	0,
 0,	0,	0,	0,	0,	0,	0,	1,
 1,	0,	1,	1,	1,	0,	0,	0,
 0,	1,	0,	1,	1,	1,	0,	0,
 0,	0,	1,	0,	1,	1,	1,	0,
 0,	0,	0,	1,	0,	1,	1,	1,
 1,	0,	1,	1,	0,	0,	1,	1,
 1,	1,	1,	0,	0,	0,	0,	1,
 1,	1,	0,	0,	1,	0,	0,	0,
 0,	1,	1,	0,	0,	1,	0,	0,
 0,	0,	1,	1,	0,	0,	1,	0,
 0,	0,	0,	1,	1,	0,	0,	1,
 1,	0,	1,	1,	0,	1,	0,	0,
 0,	1,	0,	1,	1,	0,	1,	0,
 0,	0,	1,	0,	1,	1,	0,	1,
 1,	0,	1,	0,	1,	1,	1,	0,
 0,	1,	0,	1,	0,	1,	1,	1,
 1,	0,	0,	1,	0,	0,	1,	1,
 1,	1,	1,	1,	0,	0,	0,	1,
 1,	1,	0,	0,	0,	0,	0,	0,
 0,	1,	1,	0,	0,	0,	0,	0,
 0,	0,	1,	1,	0,	0,	0,	0,
 0,	0,	0,	1,	1,	0,	0,	0,
 0,	0,	0,	0,	1,	1,	0,	0,
 0,	0,	0,	0,	0,	1,	1,	0,
 0,	0,	0,	0,	0,	0,	1,	1,
 1,	0,	1,	1,	1,	0,	0,	1,
 1,	1,	1,	0,	0,	1,	0,	0,
 0,	1,	1,	1,	0,	0,	1,	0,
 0,	0,	1,	1,	1,	0,	0,	1,
 1,	0,	1,	0,	0,	1,	0,	0,
 0,	1,	0,	1,	0,	0,	1,	0,
 0,	0,	1,	0,	1,	0,	0,	1,
 1,	0,	1,	0,	1,	1,	0,	0,
 0,	1,	0,	1,	0,	1,	1,	0,
 0,	0,	1,	0,	1,	0,	1,	1,
 1,	0,	1,	0,	1,	1,	0,	1,
 1,	1,	1,	0,	1,	1,	1,	0,
 0,	1,	1,	1,	0,	1,	1,	1,
 1,	0,	0,	0,	0,	0,	1,	1,
 1,	1,	1,	1,	1,	0,	0,	1,
 1,	1,	0,	0,	0,	1,	0,	0,
 0,	1,	1,	0,	0,	0,	1,	0,
 0,	0,	1,	1,	0,	0,	0,	1,
 1,	0,	1,	0,	0,	0,	0,	0,
 0,	1,	0,	1,	0,	0,	0,	0,
 0,	0,	1,	0,	1,	0,	0,	0,
 0,	0,	0,	1,	0,	1,	0,	0,
 0,	0,	0,	0,	1,	0,	1,	0,
 0,	0,	0,	0,	0,	1,	0,	1,
 1,	0,	1,	1,	1,	0,	1,	0,
 0,	1,	0,	1,	1,	1,	0,	1,
 1,	0,	0,	1,	0,	1,	1,	0,
 0,	1,	0,	0,	1,	0,	1,	1,
 1,	0,	0,	1,	1,	1,	0,	1,
 1,	1,	1,	1,	0,	1,	1,	0,
 0,	1,	1,	1,	1,	0,	1,	1,
 1,	0,	0,	0,	0,	1,	0,	1,
 1,	1,	1,	1,	1,	0,	1,	0,
 0,	1,	1,	1,	1,	1,	0,	1,
 1,	0,	0,	0,	0,	1,	1,	0,
 0,	1,	0,	0,	0,	0,	1,	1,
 1,	0,	0,	1,	1,	0,	0,	1,
 1,	1,	1,	1,	0,	1,	0,	0,
 0,	1,	1,	1,	1,	0,	1,	0,
 0,	0,	1,	1,	1,	1,	0,	1,
 1,	0,	1,	0,	0,	1,	1,	0,
 0,	1,	0,	1,	0,	0,	1,	1,
 1,	0,	0,	1,	0,	0,	0,	1,
 1,	1,	1,	1,	0,	0,	0,	0,
 0,	1,	1,	1,	1,	0,	0,	0,
 0,	0,	1,	1,	1,	1,	0,	0,
 0,	0,	0,	1,	1,	1,	1,	0,
 0,	0,	0,	0,	1,	1,	1,	1,
 1,	0,	1,	1,	1,	1,	1,	1,
 1,	1,	1,	0,	0,	1,	1,	1,
 1,	1,	0,	0,	1,	0,	1,	1,
 1,	1,	0,	1,	1,	1,	0,	1,
 1,	1,	0,	1,	0,	1,	1,	0,
 0,	1,	1,	0,	1,	0,	1,	1,
 1,	0,	0,	0,	1,	1,	0,	1,
 1,	1,	1,	1,	1,	1,	1,	0,
 0,	1,	1,	1,	1,	1,	1,	1,
 1,	0,	0,	0,	0,	1,	1,	1,
 1,	1,	1,	1,	1,	0,	1,	1,
 1,	1,	0,	0,	0,	1,	0,	1,
 1,	1,	0,	1,	1,	0,	1,	0,
 0,	1,	1,	0,	1,	1,	0,	1,
 1,	0,	0,	0,	1,	1,	1,	0,
 0,	1,	0,	0,	0,	1,	1,	1,
 1,	0,	0,	1,	1,	0,	1,	1,
 1,	1,	1,	1,	0,	1,	0,	1,
 1,	1,	0,	0,	0,	0,	1,	0,
 0,	1,	1,	0,	0,	0,	0,	1,
 1,	0,	0,	0,	1,	0,	0,	0,
 0,	1,	0,	0,	0,	1,	0,	0,
 0,	0,	1,	0,	0,	0,	1,	0,
 0,	0,	0,	1,	0,	0,	0,	1,
 1,	0,	1,	1,	0,	0,	0,	0,
 0,	1,	0,	1,	1,	0,	0,	0,
 0,	0,	1,	0,	1,	1,	0,	0,
 0,	0,	0,	1,	0,	1,	1,	0,
 0,	0,	0,	0,	1,	0,	1,	1,
 1,	0,	1,	1,	1,	1,	0,	1,
 1,	1,	1,	0,	0,	1,	1,	0,
 0,	1,	1,	1,	0,	0,	1,	1,
 1,	0,	0,	0,	0,	0,	0,	1,
 1,	1,	1,	1,	1,	0,	0,	0,
 0,	1,	1,	1,	1,	1,	0,	0,
 0,	0,	1,	1,	1,	1,	1,	0,
 0,	0,	0,	1,	1,	1,	1,	1,
 1,	0,	1,	1,	0,	1,	1,	1,
 1,	1,	1,	0,	0,	0,	1,	1,
 1,	1,	0,	0,	1,	0,	0,	1,
 1,	1,	0,	1,	1,	1,	0,	0,
 0,	1,	1,	0,	1,	1,	1,	0,
 0,	0,	1,	1,	0,	1,	1,	1,
 1,	0,	1,	0,	0,	0,	1,	1,
 1,	1,	1,	0,	1,	0,	0,	1,
 1,	1,	0,	0,	1,	1,	0,	0,
 0,	1,	1,	0,	0,	1,	1,	0,
 0,	0,	1,	1,	0,	0,	1,	1,
 1,	0,	1,	0,	0,	0,	0,	1,
 1,	1,	1,	0,	1,	0,	0,	0,
 0,	1,	1,	1,	0,	1,	0,	0,
 0,	0,	1,	1,	1,	0,	1,	0,
 0,	0,	0,	1,	1,	1,	0,	1,
 1,	0,	1,	1,	0,	1,	1,	0,
 0,	1,	0,	1,	1,	0,	1,	1,
 1,	0,	0,	1,	0,	1,	0,	1,
 1,	1,	1,	1,	0,	0,	1,	0,
 0,	1,	1,	1,	1,	0,	0,	1,
 1,	0,	0,	0,	0,	1,	0,	0,
 0,	1,	0,	0,	0,	0,	1,	0,
 0,	0,	1,	0,	0,	0,	0,	1,
 1,	0,	1,	0,	1,	0,	0,	0,
 0,	1,	0,	1,	0,	1,	0,	0,
 0,	0,	1,	0,	1,	0,	1,	0,
 0,	0,	0,	1,	0,	1,	0,	1,
 1,	0,	1,	1,	0,	0,	1,	0,
 0,	1,	0,	1,	1,	0,	0,	1,
 1,	0,	0,	1,	0,	1,	0,	0,
 0,	1,	0,	0,	1,	0,	1,	0,
 0,	0,	1,	0,	0,	1,	0,	1,
 1,	0,	1,	0,	1,	0,	1,	0,
 0,	1,	0,	1,	0,	1,	0,	1,
 1,	0,	0,	1,	0,	0,	1,	0,
 0,	1,	0,	0,	1,	0,	0,	1,
 1,	0,	0,	1,	1,	1,	0,	0,
 0,	1,	0,	0,	1,	1,	1,	0,
 0,	0,	1,	0,	0,	1,	1,	1,
 1,	0,	1,	0,	1,	0,	1,	1,
 1,	1,	1,	0,	1,	1,	0,	1,
 1,	1,	0,	0,	1,	1,	1,	0,
 0,	1,	1,	0,	0,	1,	1,	1,
 1,	0,	0,	0,	1,	0,	1,	1,
 1,	1,	1,	1,	1,	1,	0,	1,
 1,	1,	0,	0,	0,	1,	1,	0,
 0,	1,	1,	0,	0,	0,	1,	1,
 1,	0,	0,	0,	1,	0,	0,	1,
 1,	1,	1,	1,	1,	1,	0,	0,
 0,	1,	1,	1,	1,	1,	1,	0,
 0,	0,	1,	1,	1,	1,	1,	1,
 1,	0,	1,	0,	0,	1,	1,	1,
 1,	1,	1,	0,	1,	0,	1,	1,
 1,	1,	0,	0,	1,	1,	0,	1,
 1,	1,	0,	1,	1,	1,	1,	0,
 0,	1,	1,	0,	1,	1,	1,	1,
 1,	0,	0,	0,	1,	1,	1,	1,
 1,	1,	1,	1,	1,	1,	1,	1,
 1,	1,	0,	0,	0,	1,	1,	1,
 1,	1,	0,	1,	1,	0,	1,	1,
 1,	1,	0,	1,	0,	1,	0,	1,
 1,	1,	0,	1,	0,	0,	1,	0,
 0,	1,	1,	0,	1,	0,	0,	1,
 1,	0,	0,	0,	1,	1,	0,	0,
 0,	1,	0,	0,	0,	1,	1,	0,
 0,	0,	1,	0,	0,	0,	1,	1,
 1,	0,	1,	0,	1,	0,	0,	1,
 1,	1,	1,	0,	1,	1,	0,	0,
 0,	1,	1,	1,	0,	1,	1,	0,
 0,	0,	1,	1,	1,	0,	1,	1,
 1,	0,	1,	0,	0,	1,	0,	1,
 1,	1,	1,	0,	1,	0,	1,	0,
 0,	1,	1,	1,	0,	1,	0,	1,
 1,	0,	0,	0,	0,	0,	1,	0,
 0,	1,	0,	0,	0,	0,	0,	1,
 1,	0,	0,	1,	1,	0,	0,	0,
 0,	1,	0,	0,	1,	1,	0,	0,
 0,	0,	1,	0,	0,	1,	1,	0,
 0,	0,	0,	1,	0,	0,	1,	1,
 1,	0,	1,	1,	0,	0,	0,	1,
 1,	1,	1,	0,	0,	0,	0,	0,
 0,	1,	1,	1,	0,	0,	0,	0,
 0,	0,	1,	1,	1,	0,	0,	0,
 0,	0,	0,	1,	1,	1,	0,	0,
 0,	0,	0,	0,	1,	1,	1,	0,
 0,	0,	0,	0,	0,	1,	1,	1,
 1,	0,	1,	1,	1,	0,	1,	1,
 1,	1,	1,	0,	0,	1,	0,	1,
 1,	1,	0,	0,	1,	0,	1,	0,
 0,	1,	1,	0,	0,	1,	0,	1,
 1,	0,	0,	0,	1,	0,	1,	0,
 0,	1,	0,	0,	0,	1,	0,	1,
 1,	0,	0,	1,	1,	0,	1,	0,
 0,	1,	0,	0,	1,	1,	0,	1,
 1,	0,	0,	1,	1,	1,	1,	0,
 0,	1,	0,	0,	1,	1,	1,	1,
 1,	0,	0,	1,	1,	1,	1,	1,
 1,	1,	1,	1,	0,	1,	1,	1,
 1,	1,	0,	0,	0,	0,	1,	1,
 1,	1,	0,	1,	1,	0,	0,	1,
 1,	1,	0,	1,	0,	1,	0,	0,
 0,	1,	1,	0,	1,	0,	1,	0,
 0,	0,	1,	1,	0,	1,	0,	1,
 1,	0,	1,	0,	0,	0,	1,	0,
 0,	1,	0,	1,	0,	0,	0,	1,
 1,	0,	0,	1,	0,	0,	0,	0,
 0,	1,	0,	0,	1,	0,	0,	0,
 0,	0,	1,	0,	0,	1,	0,	0,
 0,	0,	0,	1,	0,	0,	1,	0,
 0,	0,	0,	0,	1,	0,	0,	1,
 1,	0,	1,	1,	1,	1,	0,	0,
 0,	1,	0,	1,	1,	1,	1,	0,
 0,	0,	1,	0,	1,	1,	1,	1,
 1,	0,	1,	0,	1,	1,	1,	1,
 1,	1,	1,	0,	1,	1,	1,	1,
 1,	1,	0,	0,	1,	1,	1,	1,
 1,	1,	0,	1,	1,	1,	1,	1,
 1,	1,	0,	1,	0,	1,	1,	1,
 1,	1,	0,	1,	0,	0,	1,	1,
 1,	1,	0,	1,	0,	0,	0,	1,
 1,	1,	0,	1,	0,	0,	0,	0,
 0,	1,	1,	0,	1,	0,	0,	0,
 0,	0,	1,	1,	0,	1,	0,	0,
 0,	0,	0,	1,	1,	0,	1,	0,
 0,	0,	0,	0,	1,	1,	0,	1,
 1,	0,	1,	1,	1,	1,	1,	0,
 0,	1,	0,	1,	1,	1,	1,	1,
 1,	0,	0,	1,	0,	1,	1,	1,
 1,	1,	1,	1,	0,	0,	1,	1,
 1,	1,	0,	0,	0,	0,	0,	1,
 1,	1,	0,	1,	1,	0,	0,	0,
 0,	1,	1,	0,	1,	1,	0,	0,
 0,	0,	1,	1,	0,	1,	1,	0,
 0,	0,	0,	1,	1,	0,	1,	1,
 1,	0,	1,	1,	0,	1,	0,	1,
 1,	1,	1,	0,	0,	0,	1,	0,
 0,	1,	1,	1,	0,	0,	0,	1,};

static const int BinGF_64[64][6]={
	{0,0,0,0,0,0},
	{1,0,0,0,0,0},
	{0,1,0,0,0,0},
	{0,0,1,0,0,0},
	{0,0,0,1,0,0},
	{0,0,0,0,1,0},
	{0,0,0,0,0,1},
	{1,1,0,0,0,0},
	{0,1,1,0,0,0},
	{0,0,1,1,0,0},
	{0,0,0,1,1,0},
	{0,0,0,0,1,1},
	{1,1,0,0,0,1},
	{1,0,1,0,0,0},
	{0,1,0,1,0,0},
	{0,0,1,0,1,0},
	{0,0,0,1,0,1},
	{1,1,0,0,1,0},
	{0,1,1,0,0,1},
	{1,1,1,1,0,0},
	{0,1,1,1,1,0},
	{0,0,1,1,1,1},
	{1,1,0,1,1,1},
	{1,0,1,0,1,1},
	{1,0,0,1,0,1},
	{1,0,0,0,1,0},
	{0,1,0,0,0,1},
	{1,1,1,0,0,0},
	{0,1,1,1,0,0},
	{0,0,1,1,1,0},
	{0,0,0,1,1,1},
	{1,1,0,0,1,1},
	{1,0,1,0,0,1},
	{1,0,0,1,0,0},
	{0,1,0,0,1,0},
	{0,0,1,0,0,1},
	{1,1,0,1,0,0},
	{0,1,1,0,1,0},
	{0,0,1,1,0,1},
	{1,1,0,1,1,0},
	{0,1,1,0,1,1},
	{1,1,1,1,0,1},
	{1,0,1,1,1,0},
	{0,1,0,1,1,1},
	{1,1,1,0,1,1},
	{1,0,1,1,0,1},
	{1,0,0,1,1,0},
	{0,1,0,0,1,1},
	{1,1,1,0,0,1},
	{1,0,1,1,0,0},
	{0,1,0,1,1,0},
	{0,0,1,0,1,1},
	{1,1,0,1,0,1},
	{1,0,1,0,1,0},
	{0,1,0,1,0,1},
	{1,1,1,0,1,0},
	{0,1,1,1,0,1},
	{1,1,1,1,1,0},
	{0,1,1,1,1,1},
	{1,1,1,1,1,1},
	{1,0,1,1,1,1},
	{1,0,0,1,1,1},
	{1,0,0,0,1,1},
	{1,0,0,0,0,1},
};




/**
 * Function name : Bin2GF
 * Description   : compute a GF(q) symbol corresponding to a frame of log_2(GF) bits
 * Parameters    :
 * Inputs        :
 * 	- int *U    : array representing logGF bits
 * 	- int logGF : size of the array U. logGF = log2 (GF)
 * 	- int GF    : order of the field
 * 	- int ** BINGF: binary mapping table
 * Outputs       :
 *      - index of the non-binary symbol
 */

 int Bin2GF(int *U,int GF,int logGF,int **BINGF)
{
	int k;

	for (k=0;k<GF;k++) {
		if (memcmp(U,BINGF[k],sizeof(int)*logGF)==0) break;
	}
	return(k);
}

/**
 * Function name : Table_Add_GF
 * Description   : compute the addition table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
static void Table_Add_GF(table_t *table, int GF, int logGF)
{
	int i,j,k,m;
	int temp[8];

	for(j=0; j<GF; j++) {
		for(k=0; k<GF; k++) {
			for(i=0; i<logGF; i++) {
				temp[i] = (table->BINGF[j][i])^(table->BINGF[k][i]);
			}
			table->ADDGF[j][k] = Bin2GF(temp,GF,logGF,table->BINGF);
		}
	}
	//for(k=0;k<256;k++){for(m=0;m<256;m++){printf("%d ",table->ADDGF[m][k]);}printf("\n");}getchar();

}

/**
 * Function name : Table_Mul_GF
 * Description   : compute the multiplication table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
static void Table_Mul_GF(int **MULGF, int GF)
{
	int i,j,temp;
	for(i=0; i<GF; i++) {
		for(j=0;j<GF;j++) {
			if (i==0 || j==0)
				MULGF[i][j] = 0;
			else if (i==1)
				MULGF[i][j] = j;
			else if (j==1)
				MULGF[i][j] = i;
			else {
				temp=i+j-2;
				if(temp<GF-1)
					MULGF[i][j] = temp+1;
				else
					MULGF[i][j]=(temp%(GF-1))+1;
			}
		}
	}
}

/**
 * Function name : Table_Div_GF
 * Description   : compute the division table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
static void Table_Div_GF(int **DIVGF, int GF)
{
	int i,j,nb;
	nb=GF-1;
	for(i=0;i<GF;i++) {
		for(j=0;j<GF;j++) {
			if(j==0){DIVGF[i][j]=-1;}
			else if(i==0){ DIVGF[i][j]=0;}
			else if(j==1){ DIVGF[i][j]=i;}
			else         { DIVGF[i][j]=nb--;};
			if(nb<1){nb=GF-1;}
		}
	}
}

/**
 * Function name : LoadCode
 * Description   : Open a parity-check matrix file. Allocate and initialize the code structures
 * Parameters    :
 * Inputs        :
 * 	- FileMatrix  : Name of the parity-check matrix file
 * 	- code_t code : Structure that describes the code
 * Outputs       :
 */
void LoadCode (char *FileMatrix, code_t *code)
{
	int M,N,GF,logGF;
	int numColumn;
	int n,m,k;
	char *FileName = malloc(STR_MAXSIZE);
	int *ind, numBranch;
        //int **NtoB,
	double temp;
	FILE *f;
	/*
	 * Load the files corresponding to code (graph, size, GF)
	 */
	strcpy(FileName,FileMatrix);
	f=fopen(FileName,"r");

	fscanf(f,"%d",&N);
	fscanf(f,"%d",&M);
	fscanf(f,"%d",&GF);
	temp = log((double)(GF));
	temp = temp/log((double)2.0);
	temp = rint(temp);
	logGF = (int)temp;
	code->N = N;
	code->M = M;
	code->K = N-M;
	code->rate=(float)(N-M)/N;
	code->GF = GF;
	code->logGF = logGF;


    code->columnDegree = calloc(N,sizeof(int));
	if ( code->columnDegree == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (n=0;n<N;n++) fscanf(f,"%d",&code->columnDegree[n]);

	code->rowDegree = calloc(M,sizeof(int));
	if ( code->rowDegree == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (m=0;m<M;m++) fscanf(f,"%d",&code->rowDegree[m]);

	code->mat = calloc(M,sizeof(int *));
	if ( code->mat == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (m=0; m<M; m++) {
		code->mat[m] = calloc(code->rowDegree[m],sizeof(int));
		if ( code->mat[m] == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	}
	for (m=0; m<M; m++)
		for (k=0; k<code->rowDegree[m]; k++)
			fscanf(f,"%d",&code->mat[m][k]);

	code->matValue = calloc(M,sizeof(int *));
	for (m=0;m<M;m++) {
		code->matValue[m] = calloc(code->rowDegree[m],sizeof(int));
		if ( code->matValue [m] == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	}
	for (m=0;m<M;m++)
		for (k=0;k<code->rowDegree[m];k++)
			fscanf(f,"%d",&code->matValue[m][k]);

	fclose(f);
	/*
	 * Build  the interleaver:
	 * We label of the ones in the parity-check matrix in a 'row per row' direction
	 * The interleaver contains the order of the labels when reading in a 'column per column' direction
	 */
	code->nbBranch=0;
	for (m=0;m<M;m++) code->nbBranch += code->rowDegree[m];

	code->interleaver = calloc(code->nbBranch,sizeof(int));

        code->NtoB=calloc(N,sizeof(int*));
        if ( code->NtoB == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (n=0; n<N; n++) {
		code->NtoB[n] = calloc(code->columnDegree[n],sizeof(int));
		if ( code->NtoB[n] == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	}

	ind = (int *)calloc(N,sizeof(int));
	numBranch=0;
	for (m=0;m<M;m++) {
		for (k=0;k<code->rowDegree[m];k++) {
			numColumn = code->mat[m][k];
			code->NtoB[numColumn][ind[numColumn]++] = numBranch++;
		}
	}
	numBranch=0;
	for (n=0;n<N;n++) {
		for (k=0;k<code->columnDegree[n];k++)
			code->interleaver[numBranch++] =code->NtoB[n][k];
	}

	printf("LDPC code parameters: N=%d  K=%d  M=%d  CodeRate=%g  GF=%d  logGF=%d\n",N,N-M,M,code->rate,GF,logGF);
	fflush(stdout);
	free(FileName);
	free(ind);
}

/**
 * Function name : FreeCode
 * Description   : Free pointers in a code_t structure
 * Inputs        :
 * 	- code_t code : Structure that describes the code
 * Outputs       :
 */
void FreeCode(code_t *code)
{
	int m,n;
	for (m=0; m<code->M; m++) {
		free(code->mat[m]);
		free(code->matValue[m]);
	}
        for (n=0; n<code->N;n++) free(code->NtoB[n]);
        free(code->NtoB);
	free(code->matUT[0]);
	free(code->matUT);
	free(code->mat );
	free(code->matValue );
	free(code->rowDegree );
	free(code->columnDegree );
	free(code->Perm );
	free(code->interleaver );
}

/*
 * Memory allocations for the decoder
 */
void AllocateDecoder (code_t *code, decoder_t *decoder)
{
	const int N = code->N;
	int nbRow, nbCol, k, l, nbRow_arr;

	decoder->nbBranch 	= code->nbBranch;
	decoder->N 		= code->N;

	nbRow = code->nbBranch; nbCol = decoder->nbMax;
        nbRow_arr=code->rowDegree[0];
	/* VtoC [nbBranch][nbMax] */
	decoder->VtoC = calloc((size_t)nbRow,sizeof(softdata_t *));
	if (decoder->VtoC  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->VtoC [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
	if (decoder->VtoC [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->VtoC[k] = decoder->VtoC[0] + k*nbCol;

	/* CtoV [nbBranch][nbMax] */
	decoder->CtoV =calloc((size_t)nbRow,sizeof(softdata_t *));
	if (decoder->CtoV  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->CtoV [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
	if (decoder->CtoV [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->CtoV[k] = decoder->CtoV[0] + k*nbCol;

	/* VtoCindex [nbBranch][nbMax] */
	decoder->VtoCindex =calloc((size_t)nbRow,sizeof(int *));
	if (decoder->VtoCindex  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->VtoCindex [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
	if (decoder->VtoCindex [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->VtoCindex[k] = decoder->VtoCindex[0] + k*nbCol;

	/* CtoVindex [nbBranch][nbMax] */
	decoder->CtoVindex =calloc((size_t)nbRow,sizeof(int *));
	if (decoder->CtoVindex  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->CtoVindex [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
	if (decoder->CtoVindex [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->CtoVindex[k] = decoder->CtoVindex[0] + k*nbCol;


        /* CV_arr [nbBranch][nbMax] */
	decoder->CV_arr =calloc((size_t)nbRow_arr,sizeof(softdata_t *));
	if (decoder->CV_arr  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->CV_arr [0] = calloc((size_t)nbRow_arr*nbCol,sizeof(softdata_t));
	if (decoder->CV_arr [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow_arr; k++) decoder->CV_arr[k] = decoder->CV_arr[0] + k*nbCol;

       /* VC_arr [nbBranch][nbMax] */
	decoder->VC_arr =calloc((size_t)nbRow_arr,sizeof(softdata_t *));
	if (decoder->VC_arr  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->VC_arr [0] = calloc((size_t)nbRow_arr*nbCol,sizeof(softdata_t));
	if (decoder->VC_arr [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow_arr; k++) decoder->VC_arr[k] = decoder->VC_arr[0] + k*nbCol;

        /* IndiceCV_arr [rowDegree[0]][nbMax] */
	decoder->IndiceCV_arr =calloc((size_t)nbRow_arr,sizeof(int *));
	if (decoder->IndiceCV_arr  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->IndiceCV_arr [0] = calloc((size_t)nbRow_arr*nbCol,sizeof(int));
	if (decoder->IndiceCV_arr [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow_arr; k++) decoder->IndiceCV_arr[k] = decoder->IndiceCV_arr[0] + k*nbCol;
        /* IndiceVC_arr [rowDegree[0]][nbMax] */
	decoder->IndiceVC_arr =calloc((size_t)nbRow_arr,sizeof(int *));
	if (decoder->IndiceVC_arr  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->IndiceVC_arr [0] = calloc((size_t)nbRow_arr*nbCol,sizeof(int));
	if (decoder->IndiceVC_arr [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow_arr; k++) decoder->IndiceVC_arr[k] = decoder->IndiceVC_arr[0] + k*nbCol;


        /* CV[nbMax] */

        decoder->CV = calloc(nbCol,sizeof(softdata_t));
	if ( decoder->CV == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);

       /* VC[nbMax] */
	decoder->VC = calloc(nbCol,sizeof(softdata_t));
	if ( decoder->VC == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);

        /* IndiceCV[nbMax] */
	decoder->IndiceCV = calloc(nbCol,sizeof(int));
	if ( decoder->IndiceCV == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);


        /* IndiceVC[nbMax] */
	decoder->IndiceVC = calloc(nbCol,sizeof(int));
	if ( decoder->IndiceVC == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);


	nbRow = N; nbCol = code->GF;
	/* APP [N][GF] */
	decoder->APP =calloc((size_t)nbRow,sizeof(softdata_t *));
	if (decoder->APP  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->APP [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
	if (decoder->APP [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->APP[k] = decoder->APP[0] + k*nbCol;

	nbRow = N; nbCol = code->GF;
	/* likelihood [N][GF] */   /* VN modified */
	decoder->likelihood =calloc((size_t)nbRow,sizeof(softdata_t *));
	if (decoder->likelihood  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->likelihood [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
	if (decoder->likelihood [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->likelihood[k] = decoder->likelihood[0] + k*nbCol;

	/* likelihood_index [N][GF] */  /* VN modified */
	decoder->likelihood_index 	= calloc((size_t)nbRow,sizeof(int *));
	if (decoder->likelihood_index  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	decoder->likelihood_index [0] 	= calloc((size_t)nbRow*nbCol,sizeof(int));
	if (decoder->likelihood_index [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) decoder->likelihood_index[k] = decoder->likelihood_index[0] + k*nbCol;

	/*Initialisation vectors VtoCindex CtoVindex*/
	for (l=0; l<code->nbBranch; l++)
		for (k=0; k<decoder->nbMax; k++) {
			decoder->VtoCindex[l][k]=-1;
			decoder->CtoVindex[l][k]=-1;
		}
}

/**
 * Function name : ClearDecoder
 * Description   : Reinitialize all the edges in the decoder structure.
 * Inputs        :
 * 	- decoder_t decoder : Structure that describes the decoder (edges in the decoding graph)
 * 	- int GF  : order of the field
 * Outputs       :
 */
void ClearDecoder (decoder_t *decoder, int GF)
{
	int  k, l, n;
	/* VtoC [nbBranch][nbMax] */
	/* CtoV [nbBranch][nbMax] */
	/* VtoCindex [nbBranch][nbMax] */
	/* CtoVindex [nbBranch][nbMax] */
	for (k=0; k<decoder->nbBranch; k++) {
		for (l=0; l<decoder->nbMax; l++) {
			decoder->VtoC[k][l] = 0;
			decoder->CtoV[k][l] = 0;
			decoder->VtoCindex[k][l]=-1;
			decoder->CtoVindex[k][l]=-1;
		}
	}

	/* likelihood [N][nbMax] */
	/* likelihood_index [N][nbMax] */
	for (n=0; n<decoder->N; n++) {
		for (l=0; l<GF; l++) {   /* VN modified */
			decoder->likelihood[n][l] = 0;
			decoder->likelihood_index[n][l] = -1;
		}
		/* APP [N][GF] */
		for (l=0; l<GF; l++)
			decoder->APP[n][l] = 0;
	}
}

/**
 * Function name : FreeDecoder
 * Description   : Free pointers in a decoder_t structure
 * Inputs        :
 * 	- decoder_t decoder : Structure that describes the decoder (edges in the decoding graph)
 * Outputs       :
 */
void FreeDecoder (decoder_t *decoder)
{
	free(decoder->VtoC[0] );free(decoder->VtoC);
	free(decoder->CtoV[0] );free(decoder->CtoV);
	free(decoder->VtoCindex[0] );free(decoder->VtoCindex);
	free(decoder->CtoVindex[0] );free(decoder->CtoVindex);
        free(decoder->CV_arr[0] );free(decoder->CV_arr);
        free(decoder->VC_arr[0] );free(decoder->VC_arr);
        free(decoder->IndiceCV_arr[0] );free(decoder->IndiceCV_arr);
        free(decoder->IndiceVC_arr[0] );free(decoder->IndiceVC_arr);
	free(decoder->CV);
        free(decoder->VC);
	free(decoder->IndiceCV);
        free(decoder->IndiceVC);

        free(decoder->APP[0] );free(decoder->APP);
	free(decoder->likelihood[0] );free(decoder->likelihood);
	free(decoder->likelihood_index[0] );free(decoder->likelihood_index);

}


/**
 * Function name : LoadTables
 * Description   : Memory allocation for the tables and Initialization of the tables.
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * 	- int GF    : order of the field
 * 	- int logGF : logGF = log2(GF)
 * Outputs       :
 */
void LoadTables (table_t *table, int GF, int logGF)
{
	int nbRow, nbCol, g,k,l;

	if(GF!=64 && GF!=256)	{
		printf("The binary image of GF(%d) is not available in this version of the program. Please try GF(64) or GF(256)\n",GF);
		exit(EXIT_FAILURE);
	}

	nbRow = GF; nbCol = logGF;
	/* BINGF [GF][logGF] */
	table->BINGF =calloc((size_t)nbRow,sizeof(int *));
	if (table->BINGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	table->BINGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
	if (table->BINGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) table->BINGF[k] = table->BINGF[0] + k*nbCol;

	nbRow = GF; nbCol = GF;
	/* ADDGF [GF][logGF] */
	table->ADDGF =calloc((size_t)nbRow,sizeof(int *));
	if (table->ADDGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	table->ADDGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
	if (table->ADDGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) table->ADDGF[k] = table->ADDGF[0] + k*nbCol;

	/* MULGF [GF][logGF] */
	table->MULGF =calloc((size_t)nbRow,sizeof(int *));
	if (table->MULGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	table->MULGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
	if (table->MULGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) table->MULGF[k] = table->MULGF[0] + k*nbCol;

	/* DIVGF [GF][logGF] */
	table->DIVGF =calloc((size_t)nbRow,sizeof(int *));
	if (table->DIVGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	table->DIVGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
	if (table->DIVGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<nbRow; k++) table->DIVGF[k] = table->DIVGF[0] + k*nbCol;

	if(GF==256) {
		for(g=0;g<GF;g++)
			for(l=0;l<logGF;l++)
				table->BINGF[g][l] = BinGF_256[g][l];
		printf("Loading of the binary image of GF(64): Success\n");
		fflush(stdout);
	}

	/*
	 * Build the addition, multiplication and division tables (corresponding to GF[q])
	 */
	Table_Add_GF (table,GF,logGF);
	Table_Mul_GF (table->MULGF,GF);
	Table_Div_GF (table->DIVGF,GF);

}


/**
 * Function name : FreeTable
 * Description   : Free the memory in a table_t structure
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * Outputs       :
 */
void FreeTable(table_t *table)
{
	free(table->ADDGF[0]);
	free(table->MULGF[0]);
	free(table->DIVGF[0]);
	free(table->BINGF[0]);
	free(table->ADDGF);
	free(table->MULGF);
	free(table->DIVGF);
	free(table->BINGF);
}

/**
 * Function name : RandomBinaryGenerator
 * Description   : Uniform random binary generator (generate the information bits of the code (KBIN))
 * Description   : Uniform random binary generator (generate the information bits of the code (KBIN))
 * 		   The redundancy symbols are initialized to 0.
 * Inputs
 * 	- N 	: Code length
 * 	- M 	: Number of parity non-binary symbols
 * 	- GF    : Order of the field
 * 	- logGF : logGF = log2(GF)
 * Outputs       :
 *      - NBIN  : Binary representation of the codeword
 *      - KIN   : Binary representation of the information symbols
 *      - NSYMB : Non-binary symbols of the codeword
 *      - KSYMB : Information non-binary symbols
 */
void RandomBinaryGenerator (int N, int M,int GF,int logGF, int **NBIN,int **KBIN, int *NSYMB,int *KSYMB,int **BINGF)
{
	int k,q;

	/* Redundancy symbols */
	for (k=0;k<M;k++) {
		for (q=0;q<logGF;q++)
			NBIN[k][q]=0;
		NSYMB[k]=0;
	}

	/* Random and (bitwise) uniformly distributed information symbols */
	for (k=M;k<N;k++) {
		for (q=0;q<logGF;q++)
			NBIN[k][q]=floor(drand48()*2.0);;

		NSYMB[k]=Bin2GF(NBIN[k],GF,logGF,BINGF);
	}

	/* Copy of information symbols */
	for (k=0;k<N-M;k++)
		KSYMB[k]=NSYMB[k+M];
	/* Binary copy of information symbols */
	for (k=0;k<N-M;k++)
		for (q=0;q<logGF;q++)
			KBIN[k][q]=NBIN[k+M][q];
}

/**
 * Function name : GaussianElimination
 * Description   : Perform a Gaussian elimination on the parity-check matrix.
 * 		   The procedure stops when the
 * 		   The redundancy symbols are initialized to 0.
 * Inputs
 * 	- code->mat   : parity-check matrix
 * 	- table       : lookup tables for computations in GF(q)
 * Outputs       :
 *      - code->matUT : Upper triangular matrix used for encoding
 *      - code->Perm : Column permutation
 */
void GaussianElimination (code_t *code, table_t *table)
{
	const int N = code->N;
	const int M = code->M;
	int n,m,k,m1,ind, buf;

	code->matUT = calloc((size_t)M,sizeof(int *));
	if (code->matUT == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	code->matUT [0] = calloc((size_t)M*N,sizeof(int));
	if (code->matUT[0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
	for (k=1; k<M; k++) code->matUT[k] = code->matUT[0] + k*N;

	code->Perm 	= calloc(N,sizeof(int));
	if (code->Perm == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);

	for (n=0;n<N;n++) code->Perm[n] = n;

	for (m=0; m<M; m++) { for (k=0; k<code->rowDegree[m]; k++) { code->matUT[m][code->mat[m][k]] = code->matValue[m][k]; } }
	for (m=0; m<M; m++) {
           // for (k=0; k<8100; k++){if(code->matUT[62][k]!=0)printf("%d ",code->matUT[62][k]);}printf("\n");
		for (ind=m; ind<N; ind++) { if (code->matUT[m][ind]!=0) break; }
		if (ind==N) {printf("The matrix is not full rank (%d,%d)\n",m,ind); exit(EXIT_FAILURE); }
		buf = code->Perm[ind]; code->Perm[ind] = code->Perm[m]; code->Perm[m] = buf;
		for (m1=0; m1<M; m1++) {
			buf = code->matUT[m1][m];
			code->matUT[m1][m] = code->matUT[m1][ind];
			code->matUT[m1][ind] = buf;
		}
        //for(k=0;k<256;k++){for(m=0;m<256;m++){printf("%d ",table->ADDGF[m][k]);}printf("\n");}getchar();
		for (m1=m+1;m1<M;m1++) {
			if (code->matUT[m1][m]!=0) {
				buf=code->matUT[m1][m];
				for (n=m;n<N;n++) {
					if (code->matUT[m1][n]!=0)
						code->matUT[m1][n] = table->DIVGF[code->matUT[m1][n]][buf];
						//printf("%d ",code->matUT[m1][n]); code->matUT[m1][n] = table->DIVGF[code->matUT[m1][n]][buf];printf("%d \n",code->matUT[m1][n]);
				}//getchar();
				for (n=m;n<N;n++) {
					if (code->matUT[m1][n]!=0)
						code->matUT[m1][n] = table->MULGF[code->matUT[m1][n]][code->matUT[m][m]];
				}
				for (n=m;n<N;n++) {
					code->matUT[m1][n] = table->ADDGF[code->matUT[m1][n]][code->matUT[m][n]];//printf("%d ",table->ADDGF[code->matUT[m1][n]][code->matUT[m][n]]);
				}
			}
		}
	}
}


/**
 * Function name : Encoding
 * Description   : Encode the information bits into a codeword.
 * 		   matUT beeing upper triangular, the backsubstitution method is used.
 * 		   The M first symbols in NSYMB are redundancy symbols (before deinterleaving)
 * Inputs
 * 	- NSYMB  (N-M last symbols of NSYMB are information symbols)
 * Outputs
 *      - Codeword
 *      - NBIN : binary copy of the codeword
 */
int Encoding (code_t *code, table_t *table, int *CodeWord, int **NBIN, int *NSYMB)
{
	const int N = code->N;
	const int M = code->M;
	const int logGF = code->logGF;
	int n,m,q,buf;

	/* Backsubstitution */
	for (m=M-1; m>=0; m--) {
		buf=0;
		for (n=m+1;n<N;n++) {
			if (code->matUT[m][n]!=0)
				buf = table->ADDGF[buf][table->MULGF[code->matUT[m][n]][NSYMB[n]]];
		}
		/* Systematic codeword (interleaved) */
		NSYMB[m] = table->DIVGF[buf][code->matUT[m][m]];
	}

	/* De-interleaving */
	for (n=0; n<N; n++)
		CodeWord[code->Perm[n]] = NSYMB[n];

	/* Binary copy of the codeword: */
	for (n=0; n<N; n++) {
		for (q=0;q<logGF;q++)
			NBIN[n][q] = table->BINGF[CodeWord[n]][q];
	}

	return(0);
}

/**
 * Function name : ModelChannelAWGN
 * Description   : BPSK modulation + AWGN noise on the codeword.
 *                 The function computes the liklihoods corresponding to the noisy observations.
 * Inputs
 * 	- NBIN : Binary copy of the codeword
 * 	- EbN  : Signal to noise ratio in terms of Eb/No (dB)
 * Outputs
 *      - decoder->likelihood
 */
void ModelChannelAWGN (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN)
{
	const int N = code->N;
	const int logGF = code->logGF;
	const int nbMax = decoder->nbMax;
	int n,k,g,q;
	float w,u,v,sigma;
	float TMP[256];
    int som;
	float **NoisyBin = calloc(N,sizeof(float *));
	for (n=0;n<N;n++) NoisyBin[n] = calloc(2,sizeof(float));
/* Binary-input AWGN channel : */
//printf("je suis la ");
 //float R1=2.8284271247;
int i;

 float R=0.433860915637;

float TABLEAU[256][2]={-7.500000, -7.500000,
-7.500000, -6.500000,
-7.500000, -5.500000,
-7.500000, -4.500000,
-7.500000, -3.500000,
-7.500000, -2.500000,
-7.500000, -1.500000,
-7.500000, -0.500000,
-7.500000, 0.500000,
-7.500000, 1.500000,
-7.500000, 2.500000,
-7.500000, 3.500000,
-7.500000, 4.500000,
-7.500000, 5.500000,
-7.500000, 6.500000,
-7.500000, 7.500000,
-6.500000, -7.500000,
-6.500000, -6.500000,
-6.500000, -5.500000,
-6.500000, -4.500000,
-6.500000, -3.500000,
-6.500000, -2.500000,
-6.500000, -1.500000,
-6.500000, -0.500000,
-6.500000, 0.500000,
-6.500000, 1.500000,
-6.500000, 2.500000,
-6.500000, 3.500000,
-6.500000, 4.500000,
-6.500000, 5.500000,
-6.500000, 6.500000,
-6.500000, 7.500000,
-5.500000, -7.500000,
-5.500000, -6.500000,
-5.500000, -5.500000,
-5.500000, -4.500000,
-5.500000, -3.500000,
-5.500000, -2.500000,
-5.500000, -1.500000,
-5.500000, -0.500000,
-5.500000, 0.500000,
-5.500000, 1.500000,
-5.500000, 2.500000,
-5.500000, 3.500000,
-5.500000, 4.500000,
-5.500000, 5.500000,
-5.500000, 6.500000,
-5.500000, 7.500000,
-4.500000, -7.500000,
-4.500000, -6.500000,
-4.500000, -5.500000,
-4.500000, -4.500000,
-4.500000, -3.500000,
-4.500000, -2.500000,
-4.500000, -1.500000,
-4.500000, -0.500000,
-4.500000, 0.500000,
-4.500000, 1.500000,
-4.500000, 2.500000,
-4.500000, 3.500000,
-4.500000, 4.500000,
-4.500000, 5.500000,
-4.500000, 6.500000,
-4.500000, 7.500000,
-3.500000, -7.500000,
-3.500000, -6.500000,
-3.500000, -5.500000,
-3.500000, -4.500000,
-3.500000, -3.500000,
-3.500000, -2.500000,
-3.500000, -1.500000,
-3.500000, -0.500000,
-3.500000, 0.500000,
-3.500000, 1.500000,
-3.500000, 2.500000,
-3.500000, 3.500000,
-3.500000, 4.500000,
-3.500000, 5.500000,
-3.500000, 6.500000,
-3.500000, 7.500000,
-2.500000, -7.500000,
-2.500000, -6.500000,
-2.500000, -5.500000,
-2.500000, -4.500000,
-2.500000, -3.500000,
-2.500000, -2.500000,
-2.500000, -1.500000,
-2.500000, -0.500000,
-2.500000, 0.500000,
-2.500000, 1.500000,
-2.500000, 2.500000,
-2.500000, 3.500000,
-2.500000, 4.500000,
-2.500000, 5.500000,
-2.500000, 6.500000,
-2.500000, 7.500000,
-1.500000, -7.500000,
-1.500000, -6.500000,
-1.500000, -5.500000,
-1.500000, -4.500000,
-1.500000, -3.500000,
-1.500000, -2.500000,
-1.500000, -1.500000,
-1.500000, -0.500000,
-1.500000, 0.500000,
-1.500000, 1.500000,
-1.500000, 2.500000,
-1.500000, 3.500000,
-1.500000, 4.500000,
-1.500000, 5.500000,
-1.500000, 6.500000,
-1.500000, 7.500000,
-0.500000, -7.500000,
-0.500000, -6.500000,
-0.500000, -5.500000,
-0.500000, -4.500000,
-0.500000, -3.500000,
-0.500000, -2.500000,
-0.500000, -1.500000,
-0.500000, -0.500000,
-0.500000, 0.500000,
-0.500000, 1.500000,
-0.500000, 2.500000,
-0.500000, 3.500000,
-0.500000, 4.500000,
-0.500000, 5.500000,
-0.500000, 6.500000,
-0.500000, 7.500000,
0.500000, -7.500000,
0.500000, -6.500000,
0.500000, -5.500000,
0.500000, -4.500000,
0.500000, -3.500000,
0.500000, -2.500000,
0.500000, -1.500000,
0.500000, -0.500000,
0.500000, 0.500000,
0.500000, 1.500000,
0.500000, 2.500000,
0.500000, 3.500000,
0.500000, 4.500000,
0.500000, 5.500000,
0.500000, 6.500000,
0.500000, 7.500000,
1.500000, -7.500000,
1.500000, -6.500000,
1.500000, -5.500000,
1.500000, -4.500000,
1.500000, -3.500000,
1.500000, -2.500000,
1.500000, -1.500000,
1.500000, -0.500000,
1.500000, 0.500000,
1.500000, 1.500000,
1.500000, 2.500000,
1.500000, 3.500000,
1.500000, 4.500000,
1.500000, 5.500000,
1.500000, 6.500000,
1.500000, 7.500000,
2.500000, -7.500000,
2.500000, -6.500000,
2.500000, -5.500000,
2.500000, -4.500000,
2.500000, -3.500000,
2.500000, -2.500000,
2.500000, -1.500000,
2.500000, -0.500000,
2.500000, 0.500000,
2.500000, 1.500000,
2.500000, 2.500000,
2.500000, 3.500000,
2.500000, 4.500000,
2.500000, 5.500000,
2.500000, 6.500000,
2.500000, 7.500000,
3.500000, -7.500000,
3.500000, -6.500000,
3.500000, -5.500000,
3.500000, -4.500000,
3.500000, -3.500000,
3.500000, -2.500000,
3.500000, -1.500000,
3.500000, -0.500000,
3.500000, 0.500000,
3.500000, 1.500000,
3.500000, 2.500000,
3.500000, 3.500000,
3.500000, 4.500000,
3.500000, 5.500000,
3.500000, 6.500000,
3.500000, 7.500000,
4.500000, -7.500000,
4.500000, -6.500000,
4.500000, -5.500000,
4.500000, -4.500000,
4.500000, -3.500000,
4.500000, -2.500000,
4.500000, -1.500000,
4.500000, -0.500000,
4.500000, 0.500000,
4.500000, 1.500000,
4.500000, 2.500000,
4.500000, 3.500000,
4.500000, 4.500000,
4.500000, 5.500000,
4.500000, 6.500000,
4.500000, 7.500000,
5.500000, -7.500000,
5.500000, -6.500000,
5.500000, -5.500000,
5.500000, -4.500000,
5.500000, -3.500000,
5.500000, -2.500000,
5.500000, -1.500000,
5.500000, -0.500000,
5.500000, 0.500000,
5.500000, 1.500000,
5.500000, 2.500000,
5.500000, 3.500000,
5.500000, 4.500000,
5.500000, 5.500000,
5.500000, 6.500000,
5.500000, 7.500000,
6.500000, -7.500000,
6.500000, -6.500000,
6.500000, -5.500000,
6.500000, -4.500000,
6.500000, -3.500000,
6.500000, -2.500000,
6.500000, -1.500000,
6.500000, -0.500000,
6.500000, 0.500000,
6.500000, 1.500000,
6.500000, 2.500000,
6.500000, 3.500000,
6.500000, 4.500000,
6.500000, 5.500000,
6.500000, 6.500000,
6.500000, 7.500000,
7.500000, -7.500000,
7.500000, -6.500000,
7.500000, -5.500000,
7.500000, -4.500000,
7.500000, -3.500000,
7.500000, -2.500000,
7.500000, -1.500000,
7.500000, -0.500000,
7.500000, 0.500000,
7.500000, 1.500000,
7.500000, 2.500000,
7.500000, 3.500000,
7.500000, 4.500000,
7.500000, 5.500000,
7.500000, 6.500000,
7.500000, 7.500000, };


for(i=0; i<256; i++){TABLEAU[i][0]=R*TABLEAU[i][0];TABLEAU[i][1]=R*TABLEAU[i][1];}

//float somme1=0;
//for(i=0; i<256; i++){somme1=somme1+TABLEAU[i][0]*TABLEAU[i][0]+TABLEAU[i][1]*TABLEAU[i][1];}
//somme1=somme1/256;
//printf("%f",somme1); getchar();

float ATT[N];
	sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));
	for (n=0;n<N;n++) {
            som=0;

	    for (q=0;q<8;q++)
	    {som = som + NBIN[n][q]*pow(2,q);}
	   //printf("\n %d \n",som );
w=drand48();

float A=sqrt(-log(w));//printf(" %f %f\n",w,A);
ATT[n]=A;

		for (q=0; q<2; q++) {
			u=drand48();
			v=drand48();
			/* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
			NoisyBin[n][q] = TABLEAU[som][q]*A+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
		}
   //printf("%f %f %f %f \n",NoisyBin[n][0], TABLEAU[som][0] , NoisyBin[n][1] , TABLEAU[som][1]);getchar();
	}

	/* Compute the Log Likelihood Ratio messages */
	for (n=0; n<N; n++) {

        //printf("%d \n",code->GF);getchar();
		for(k=0; k<code->GF; k++) {
         som=0;
                 for (q=0;q<8;q++)
	    {som = som + BinGF_256[k][q]*pow(2,q);}

       if(ATT[n]>0 ) TMP[k] = -(SQR(NoisyBin[n][0]-TABLEAU[som][0]*ATT[n])/(2.0*SQR(sigma)))-SQR(NoisyBin[n][1]-ATT[n]*TABLEAU[som][1])/(2.0*SQR(sigma));
         else
        TMP[k]=0;
		//printf("%f \n", TMP[k]);
		}
		//getchar();
		for(k=0; k<code->GF; k++) {
			decoder->likelihood[n][k] = -1e5;
			decoder->likelihood_index[n][k] = -1;
			for (g=0; g<code->GF; g++)  {
				if (TMP[g] > decoder->likelihood[n][k]) {
					decoder->likelihood[n][k] = TMP[g];
					decoder->likelihood_index[n][k] = g;
				}
			}
			TMP[decoder->likelihood_index[n][k]] = -1e5;
		}

	}

	for (n=0;n<N;n++) free(NoisyBin[n]); free(NoisyBin);
}

/**
 * Function name : DataPassLog
 * Description   : Process the Variable to Check messages in the decoding graph.
 * Inputs
 * 	- decoder->CtoV
 *      - decoder->likelihood
 * Outputs
 * 	- decoder->VtoC
 */
void DataPassLog (code_t *code, decoder_t *decoder, float *TMP, int *TMPindex,float offset)
{

        const int N 	= code->N;
	const int GF 	= code->GF;
	const int nbMax = decoder->nbMax;
	int   g,k;
	float TMPP[256]; for (g=0; g<256; g++) TMPP[g] = 0;

     for (g=0; g<GF; g++) TMPP[TMPindex[g]]=TMP[g]-offset;

     for (k=0;k<(nbMax);k++)
       {
	if(decoder->IndiceCV[k]!=-1)
	  TMPP[decoder->IndiceCV[k]]=TMPP[decoder->IndiceCV[k]]+decoder->CV[k]+offset;
	else break;
       }


     /* Rearranging of the output messages */

     for(k=0;k<nbMax;k++)
       {
	decoder->VC[k]=-1e5; decoder->IndiceVC[k]=0;
	for (g=0;g<GF;g++)  {if (TMPP[g]>decoder->VC[k]) { decoder->VC[k]=TMPP[g]; decoder->IndiceVC[k]=g; } }
	TMPP[decoder->IndiceVC[k]]=-1e5;
       }

     /* Normalisation of the output messages */

     for(g=0;g<nbMax;g++) decoder->VC[g]=(decoder->VC[g]-decoder->VC[nbMax-1]);

}

int ElementaryStep(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nbMax,int nbOper, float offset){

   int nmU = nbMax; int nmV = nbMax; int nmS = nbMax;

   float loc_Output[nbMax];
   int loc_IndiceOut[nbMax];

    int u, i,j, s, ss, Indice_aux,Stp=0;
    int nl=4; // number of candidates (i.e. number of elements in the comparator)
    int pos;

   float tab_aux[nbMax][nbMax];
   float tab_comp[3][nl];
   int GFvalues_already_in_Out[GF];

    int debug = 0;


//    nl = (-1+ sqrt((float)(1+8*nmS)))/2;  // possible calculation


    //printf("nmU = %d; nmV = %d;nmS = %d;\n", nmU, nmV, nmS);
    if (debug){
    printf("nmU = %d; nmV = %d;nmS = %d;\n", nmU, nmV, nmS);
    printf("Input U = "); for (i=0; i<nmU; i++){        printf("%f ", Input1[i]);    }
    printf("\nInput V = ");    for (i=0; i<nmV; i++){        printf("%f ", Input2[i]);    }

    printf("Indice U = "); for (i=0; i<nmU; i++){        printf("%d ", IndiceInput1[i]);    }
    printf("\nIndice V = ");    for (i=0; i<nmV; i++){        printf("%d ", IndiceInput2[i]);    }
    printf("\n");
    }

    // init
    for (i=0; i<GF; i++){  GFvalues_already_in_Out[i]=-1;     }

    for (i=0; i<nbMax; i++){ // if nbOper < number of operations needed to completely fill de nbMax elements of Output, then -1000 and -1 are provided
        loc_Output[i] = -1000;
        loc_IndiceOut[i] = -1;
    }

    // tab_aux
    if (debug){   printf("\n*********** \n\n TAB_AUX \n*********** \n");}

    for (i=0; i<nmU; i++){
        for (j=0; j<nmV; j++){
        tab_aux[i][j] = Input1[i]+Input2[j];
        if (debug) printf("%.2f ", tab_aux[i][j]);
        }
        if (debug) printf("\n");
    }
    // tab_GFvalues


    if (debug){printf("\n\n*********** \n\n TAB_IndiceValues \n*********** \n");
        for (i=0; i<nmU; i++){
            for (j=0; j<nmV; j++){
                printf("%d ", ADDGF[IndiceInput1[i]][IndiceInput2[j]]);
            }
        printf("\n");
    }
}

    // tab_comp, comparator matrix, contains the "competitor elements" and their coordinates in "tab_aux"
    // Initialisation of tab_comp
    for (j=0; j<nl-1; j++){
            tab_comp[0][j]= tab_aux[j][0];
            //tab_aux[j][0]=-1;   /* for each element getting from tab_aux to tab_comp, its position in tab_aux becomes -1
                                         // This way, we control the problem described in *modified 12/02/2008*/
            tab_comp[1][j]=j;
            tab_comp[2][j]=0;
    }

	j=nl-1;
	tab_comp[0][j]= tab_aux[j-1][1];
    tab_comp[1][j]=j-1;
    tab_comp[2][j]=1;


    // debug : initial tab_comp
    if (debug){     printf("\n*********** \n initial TAB_COMP \n*********** \n");
        for (i=0; i<3; i++){
            for (j=0; j<nl; j++){
                printf("%.2f ", tab_comp[i][j]);
            }
            printf("\n");
        }
    }


	// filling Out and IndiceOut

    s = 0;

    for(ss=0; ss < nbOper; ss++){

//        printf("s=%d; ss= %d \n", s, ss);
        pos = maximum(tab_comp[0], nl);

      if ((IndiceInput1[(int)(tab_comp[1][pos])]==-1) || (IndiceInput2[(int)(tab_comp[2][pos])]==-1)){
            printf("\n \n out of bounds : IndiceInput1 = %d , IndiceInput1 = %d \n", IndiceInput1[(int)(tab_comp[1][pos])], IndiceInput2[(int)(tab_comp[2][pos])]);
          break;
      }

      Indice_aux = ADDGF[IndiceInput1[(int)(tab_comp[1][pos])]][IndiceInput2[(int)(tab_comp[2][pos])]];

        // control redundancy in the output list
        // transfer from tab_comp to output Out, IndiceOut
        if (GFvalues_already_in_Out[Indice_aux] == -1){
            loc_Output[s] = tab_comp[0][pos];
            loc_IndiceOut[s] = Indice_aux;
            GFvalues_already_in_Out[Indice_aux] = 1;
            s++;
        }

        if (s==nmS) break;

        if (0) printf("\n s=%d; pos = %d", s, pos);

        //control limits of tab_aux
        if ((tab_comp[1][pos]>= nmU-1)||(tab_comp[2][pos]>=nmV-1)) {
            break;
            printf("\n\n qqqq  out of bounds tab_aux \n");
        }

        // update tab_comp with next value from tab_aux


		u = (pos>>1);

		tab_comp[1][pos]=tab_comp[1][pos] + u;

		tab_comp[2][pos]=tab_comp[2][pos] + 1 - u;

		tab_comp[0][pos]= tab_aux[(int)tab_comp[1][pos]][(int)tab_comp[2][pos]];

		if (0)
                 { // : evolution of tab_comp (for each element getting to Out[]
        printf("\n*********** \n TAB_COMP ##%d\n*********** \n", ss);
            for (i=0; i<3; i++)
                          {
                for (j=0; j<nl; j++){ printf("% f ", tab_comp[i][j]); } printf("\n");
              }
         }


}

        /* Calcul de la longuer effective de vecteur de sortie rt Normalisation*/

          if (s==nbMax) Stp=nbMax; else Stp=s;
          for(i=0;i<Stp;i++) loc_Output[i]=(loc_Output[i]-loc_Output[Stp-1]);

          for(i=0;i<nbMax;i++) {Output[i]=loc_Output[i]; IndiceOut[i]=loc_IndiceOut[i];}

}



int maximum(float *Y, int nl){

	float aux = 0.0;
	int pos = 0;
	int i;

	for (i=0; i<nl ; i++){
		if (Y[i] > aux){
			aux = Y[i];
			pos = i;
		}
	}

return pos;

}



/**
 * Function name : CheckPassLogEMS
 * Description   : Process the Check to Variable messages in the decoding graph.
 * Inputs
 * 	- decoder->VtoC
 * Outputs
 * 	- decoder->CtoV
 */
void CheckPassLogEMS (int node,decoder_t *decoder, code_t *code, table_t *table, float offset,int NbOper)
{
  const int nbMax = decoder->nbMax;

  int        m,t,k,g,kk,S,c,k1,Stp;
  int        **MatriceInterIndice,*OutForwardIndice,*OutBackwardIndice,*OutForwardIndice1,*OutBackwardIndice1;
  float      **MatriceInter,*OutForward,*OutBackward,*OutForward1,*OutBackward1;



  //Temporary buffers (for the F/B recursivity)

  OutForward=(float *)calloc(nbMax,sizeof(float));
  OutForwardIndice=(int *)calloc(nbMax,sizeof(int));
  OutBackward=(float *)calloc(nbMax,sizeof(float));
  OutBackwardIndice=(int *)calloc(nbMax,sizeof(int));

  OutForward1=(float *)calloc(nbMax,sizeof(float));
  OutForwardIndice1=(int *)calloc(nbMax,sizeof(int));
  OutBackward1=(float *)calloc(nbMax,sizeof(float));
  OutBackwardIndice1=(int *)calloc(nbMax,sizeof(int));


      // Number of steps F/B = S
      S=2*(code->rowDegree[node]-2);

// MatriceInter - store the results (reeal values) of the F/B recursivity
      MatriceInter=(float **)calloc(S,sizeof(float*));
      for (k=0;k<S;k++) {MatriceInter[k]=(float *)calloc(nbMax,sizeof(float)); for (k1=0;k1<nbMax;k1++) {MatriceInter[k][k1]=-1e5;}}

//MatriceInterIndice - that store the results (GF(q) values) of the F/B recursivity
      MatriceInterIndice=(int **)calloc(S,sizeof(int*));
      for (k=0;k<S;k++) {MatriceInterIndice[k]=(int *)calloc(nbMax,sizeof(int)); for (k1=0;k1<nbMax;k1++) {MatriceInterIndice[k][k1]=-1;}}

//Initialization of the temporary buffers
   for (k=0;k<nbMax;k++) {OutForward[k]=-1e5;OutBackward[k]=-1e5;OutForwardIndice[k]=-1;OutBackwardIndice[k]=-1;}


//Rotation par la valeur non nulle dans la matrice pour VtoC
for(t=0;t<code->rowDegree[node];t++)
   {
     for(k=0;k<nbMax;k++)
        {
         if (decoder->IndiceVC_arr[t][k]!=-1)
           {
            c=decoder->IndiceVC_arr[t][k];
            decoder->IndiceVC_arr[t][k]=table->MULGF[c][code->matValue[node][t]];
           }
       else printf("IndiceVC_arr[%d][%d]=%d\n",t,k,decoder->IndiceVC_arr[t][k]);
        }
   }



//Initialisation algorithm
      for(k=0;k<nbMax;k++)
	{
	  OutForward[k]=decoder->VC_arr[0][k];
	  OutForwardIndice[k]=decoder->IndiceVC_arr[0][k];
	  OutBackward[k]=decoder->VC_arr[code->rowDegree[node]-1][k];
	  OutBackwardIndice[k]=decoder->IndiceVC_arr[code->rowDegree[node]-1][k];
	}

// Start of the recusivity
      for(kk=1;kk<(code->rowDegree[node]-1);kk++)
	{

           for(k=0;k<nbMax;k++)
		{
              //forward step
		OutForward1[k]=decoder->VC_arr[kk][k];
		OutForwardIndice1[k]=decoder->IndiceVC_arr[kk][k];
              //backward step
		OutBackward1[k]=decoder->VC_arr[code->rowDegree[node]-kk-1][k];
		OutBackwardIndice1[k]=decoder-> IndiceVC_arr[code->rowDegree[node]-kk-1][k];
		}

	  for(k=0;k<nbMax;k++)
	    {
	      // Storage of the intermediate vectors
	      MatriceInter[kk-1][k]=OutForward[k];
	      MatriceInterIndice[kk-1][k]=OutForwardIndice[k];
	      MatriceInter[2*(code->rowDegree[node]-2)-kk][k]=OutBackward[k];
	      MatriceInterIndice[2*(code->rowDegree[node]-2)-kk][k]=OutBackwardIndice[k];
	    }

	  //forward step
	if(kk<(code->rowDegree[node]-1))
		 ElementaryStep(OutForward,OutForward1,OutForwardIndice,OutForwardIndice1,OutForward,OutForwardIndice,table->ADDGF,code->GF,nbMax,NbOper,offset);


	  //backward step
	if(kk<(code->rowDegree[node]-1))
	      ElementaryStep(OutBackward,OutBackward1,OutBackwardIndice,OutBackwardIndice1,OutBackward,OutBackwardIndice,table->ADDGF,code->GF,nbMax,NbOper,offset);

       }


      //Update of vectors CtoV (first and last)
      for(k=0;k<nbMax;k++)
	{
	  // last vector CV_arr
	  decoder->CV_arr[code->rowDegree[node]-1][k]=OutForward[k];
	  decoder->IndiceCV_arr[code->rowDegree[node]-1][k]=OutForwardIndice[k];
	  // first vector CV_arr
	  decoder->CV_arr[0][k]=OutBackward[k];
	  decoder->IndiceCV_arr[0][k]=OutBackwardIndice[k];

	}

      //Update of the others vectors CtoV
      for(k=0;k<(code->rowDegree[node]-2);k++)
	{

         ElementaryStep(MatriceInter[k],MatriceInter[(code->rowDegree[node]-2)+k],MatriceInterIndice[k],MatriceInterIndice[(code->rowDegree[node]-2)  +k],OutForward,OutForwardIndice,table->ADDGF,code->GF,nbMax,NbOper,offset);
	  for(g=0;g<nbMax;g++) {decoder->CV_arr[k+1][g]=OutForward[g]; decoder->IndiceCV_arr[k+1][g]=OutForwardIndice[g];}
	}


//Rotation par la valeur non nulle dans la matrice pour CtoV
for(t=0;t<code->rowDegree[node];t++)
	  {
 	    for(k=0;k<nbMax;k++)
              {
                if (decoder->IndiceCV_arr[t][k]== -1)
                  {
                    Stp=k;
                    break;
                  }
                else
                    Stp=nbMax;
              }

	   for(k=0;k<Stp;k++) {decoder->IndiceCV_arr[t][k]=table->DIVGF[decoder->IndiceCV_arr[t][k]][code->matValue[node][t]]; }

	  }


      for (k=0;k<S;k++) free(MatriceInterIndice[k]); free(MatriceInterIndice);
      for (k=0;k<S;k++) free(MatriceInter[k]); free(MatriceInter);

  free(OutForward);free(OutForwardIndice);free(OutBackward);free(OutBackwardIndice);
  free(OutForward1);free(OutForwardIndice1);free(OutBackward1);free(OutBackwardIndice1);

}

/**
 * Function name : ComputeAPPLog
 * Description   : Compute the A Posteriori Probabilities
 * Inputs
 * 	- decoder->CtoV
 * 	- offset correction factor
 * Outputs
 * 	- decoder->APP
 */
void ComputeAPPLog (code_t *code, decoder_t *decoder, float offset)
{
	const int N = code->N;
	const int GF = code->GF;

	int  numB,n,g,t;
	float  **TMP;

	TMP=(float **)calloc(decoder->nbBranch,sizeof(float *));
	for (g=0;g<decoder->nbBranch;g++) 	TMP[g]=(float *)calloc(GF,sizeof(float));

	for (n=0;n<decoder->nbBranch;n++) {
		for(g=0;g<GF;g++)
			TMP[n][g]=-offset;
		for(g=0;g<decoder->nbMax;g++)
			if (decoder->CtoVindex[n][g]!=-1)
				TMP[n][decoder->CtoVindex[n][g]] = decoder->CtoV[n][g];
	}

	numB=0;
	for (n=0; n<N; n++) {
		//for(g=0;g<GF;g++)
		//	decoder->APP[n][g]=-offset;  //c'est une idiotie
		for (g=0;g<GF;g++)
			decoder->APP[n][decoder->likelihood_index[n][g]] = decoder->likelihood[n][g];

		for(g=0;g<GF;g++) {
			for (t=numB; t<numB+code->columnDegree[n]; t++)
				decoder->APP[n][g] = decoder->APP[n][g] + TMP[code->interleaver[t]][g];
		}
		numB = numB + code->columnDegree[n];
	}

	for (g=0;g<decoder->nbBranch;g++) free(TMP[g]);free(TMP);
}

/**
 * Function name : Syndrom
 * Description   : Compute the syndom of a message
 * Inputs
 * 	- code structure code_t
 * 	- table_t tableGF : lookup table
 * 	- message
 * Outputs
 * 	- synd is 0 iff. the decided message is a codeword
 * 	(the value of synd is not meaningful if synd != 0 )
 */
int Syndrom (code_t *code, int *decide, table_t *tableGF)
{
	int k,l;
	int synd;

	synd = 0;
	for (k=0; k<code->M; k++) {
		for (l=0; l<code->rowDegree[k]; l++)
			synd = tableGF->ADDGF[synd][tableGF->MULGF[code->matValue[k][l]][decide[code->mat[k][l]]]];
		if (synd != 0)
			break;
	}

	return (synd);
}


/**
 * Function name : Decision
 * Description   : Make a hard decision given the APP
 * Inputs
 * 	- APP
 * 	- N
 * 	- GF
 * Outputs
 * 	- decision
 */
void Decision( int *decision,float **APP,int N,int GF)
{
	int n,g,ind;
	float max;

	for (n=0;n<N;n++) {
		max = APP[n][0];
		ind=0;
		for (g=1;g<GF;g++)
			if (APP[n][g]>max) {
				max=APP[n][g]; ind=g;
			}
		decision[n]=ind;
	}
}

/*
 * Main program
 */
int main(int argc, char * argv[])
{
    int nb_erreurrrrrr_1 =0;
    int nb_erreurrrrrr_2 =0;
    int nb_erreurrrrrr_3 =0;
    int nb_erreurrrrrr_4 =0;
    int nb_erreurrrrrr_5 =0;

    int nbErroneousFrames_1 = 0;
    int nbErroneousFrames_2 = 0;
    int nbErroneousFrames_3 = 0;
    int nbErroneousFrames_4 = 0;
    int nbErroneousFrames_5 = 0;

	int 		t,k,l,n,iter,i;
	int 		**KBIN,*KSYMB,**NBIN,*NSYMB;
	int 		*decide,*CodeWord, *Pe;
	int 		nb,NbIterMax;
	float 		EbN;
	int 		NbOper,NbMonteCarlo;
	float 	offset;
	char *FileName,*FileMatrix,*name;
	int synd=0, nbErrors_1, nbErrors_2, nbErrors_3, nbErrors_4, nbErrors_5, nbUndetectedErrors = 0;
	int BUFF[8];
        float *TMP ;
        int *TMPindex,*Flag;


	code_t code;
	table_t table;
	decoder_t decoder;

        int node,numB,numBint ;




	/*
	 * Command line arguments
	 */
	if (argc < 8) {
		printf("Syntax:\n %s\n ",argv[0]);
		printf("\n\t NbMonteCarlo     : # simulated frames ");
		printf("\n\t NbIterMax        : # of maximum decoding iterations");
		printf("\n\t FileMatrix       : File name of the parity-check matrix");
		printf("\n\t EbN              : Eb/No (dB) ");
		printf("\n\t NbMax            : size of truncated messages ");
		printf("\n\t Offset           : offset correction factor (0.4 -- 1) ");
		printf("\n\t NbOper           : Maximum number of operations for sorting \n");
		return (EXIT_FAILURE);
	}
	FileName 	= malloc(STR_MAXSIZE);
	FileMatrix 	= malloc(STR_MAXSIZE);
	name 		= malloc(STR_MAXSIZE);


	NbMonteCarlo 		= atoi(argv[1]);
	NbIterMax 		= atoi(argv[2]);
	strcpy(FileMatrix,argv[3]);
	EbN 			= atof(argv[4]);
	decoder.nbMax 		= atoi(argv[5]);
	offset  		= atof(argv[6]);
	NbOper  		= atoi(argv[7]);
	printf("\n\t NbMonteCarlo     : %d", NbMonteCarlo);
	printf("\n\t NbIterMax        : %d", NbIterMax);
	printf("\n\t FileMatrix       : %s", FileMatrix);
	printf("\n\t Eb/No (dB)       : %g", EbN);
	printf("\n\t NbMax            : %d", decoder.nbMax);
	printf("\n\t Offset           : %g", offset);
	printf("\n\t NbOper           : %d\n",NbOper);


	Pe	= calloc(NbIterMax,sizeof(int));
	LoadCode (FileMatrix, &code);
	LoadTables (&table, code.GF, code.logGF);
	AllocateDecoder (&code, &decoder);
	GaussianElimination (&code, &table);



	/*
	 * Memory  allocation
	 */
	NBIN=(int **)calloc(code.N,sizeof(int *));
	for (n=0;n<code.N;n++)  	NBIN[n]=(int *)calloc(code.logGF,sizeof(int));
	KBIN=(int **)calloc(code.K,sizeof(int *));
	for (k=0;k<code.K;k++) 	KBIN[k]=(int *)calloc(code.logGF,sizeof(int));

	NSYMB=(int *)calloc(code.N,sizeof(int));
	KSYMB=(int *)calloc(code.K,sizeof(int));
	CodeWord=(int *)calloc(code.N,sizeof(int));

	decide=(int *)calloc(code.N,sizeof(int));

        Flag=(int *)calloc(code.N,sizeof(int));
        TMP = (float *) calloc (code.GF,sizeof(float ));
        TMPindex = (int *) calloc (code.GF,sizeof(int ));

      	for (nb=1; nb<=NbMonteCarlo; nb++) {
		/* Decoder re-initialization */
		ClearDecoder (&decoder, code.GF);

		/* Generate uniformly distributed information bits (KBIN)) */
		RandomBinaryGenerator (code.N, code.M, code.GF, code.logGF, NBIN, KBIN, NSYMB, KSYMB, table.BINGF);

		/* Encode the information bits KBIN to a (non binary) codeword NSYMB */
		Encoding (&code, &table, CodeWord, NBIN, NSYMB);

		/* Noisy channel (AWGN)*/
		ModelChannelAWGN (&code, &decoder, &table, NBIN, EbN);


      numB=0; /* numB is the edge number */


	/***********************************************/
	/* Implementation of the horizontal scheduling */

	/* initialisation step : intrinsic likelihoods from the channel go directly to check nodes before starting the first iteration */
     for (n=0;n<code.N;n++){
        for (t=0;t<code.columnDegree[n];t++){

            numBint=code.interleaver[numB+t];

            for(k=0;k<decoder.nbMax;k++){
		            decoder.VtoC[numBint][k]=decoder.likelihood[n][k]; decoder.VtoCindex[numBint][k]=decoder.likelihood_index[n][k];
            }

		   }
		  numB=numB+code.columnDegree[n];
	  }

        nbErrors_1 = 0;
        nbErrors_2 = 0;
        nbErrors_3 = 0;
        nbErrors_4 = 0;
        nbErrors_5 = 0;

		/* Decoding iterations*/
		for (iter=0; iter<=NbIterMax; iter++)
		{
		    numB=0;
            for (node=0;node<code.M;node++){ /* Loop for the M Check nodes */
                       /*******************************************************************************************************************/
                       /*******************************************************************************************************************/
				for (i=0;i<code.rowDegree[node];i++){
					for (k=0;k<decoder.nbMax;k++){
						/* Filling of the VC_arr (RowDegree[n] by NbMax array) with the corresponding VtoC */
						decoder.VC_arr[i][k]=decoder.VtoC[numB+i][k];
						decoder.IndiceVC_arr[i][k]=decoder.VtoCindex[numB+i][k];
	             }
             }

             CheckPassLogEMS (node,&decoder, &code, &table, offset, NbOper);

             for (i=0;i<code.rowDegree[node];i++){
     	          for (k=0;k<decoder.nbMax;k++){
							/* Filling of the CV_arr (RowDegree[n] by NbMax1 array) in the corresponding CtoV */
                     decoder.CtoV[numB+i][k]=decoder.CV_arr[i][k];
                     decoder.CtoVindex[numB+i][k]=decoder.IndiceCV_arr[i][k];;
                }
              }
                      /*******************************************************************************************************************/
                      /*******************************************************************************************************************/

              for (i=0;i<code.rowDegree[node];i++){
                for (k=0;k<code.GF;k++){
                    TMP[k]=decoder.likelihood[code.mat[node][i]][k];
                    TMPindex[k]=decoder.likelihood_index[code.mat[node][i]][k];
                }
                for (k=0;k<decoder.nbMax;k++){
                    decoder.CV[k]=decoder.CtoV[i+numB][k];/* Filling of the CV vector with the corresponding CtoV */
       	            decoder.IndiceCV[k]=decoder.CtoVindex[i+numB][k];
		        }

                DataPassLog (&code, &decoder,TMP,TMPindex,offset);

                Flag[code.mat[node][i]]=1^Flag[code.mat[node][i]]; /*Flag:is used to select the appropriate edge at the VN. It can be used only for dv=2*/

			     for (k=0;k<decoder.nbMax;k++){
		          decoder.VtoC[code.NtoB[code.mat[node][i]][Flag[code.mat[node][i]]]][k]=decoder.VC[k]; /* Filling of the VC in the corresponding VtoC */
		          decoder.VtoCindex[code.NtoB[code.mat[node][i]][Flag[code.mat[node][i]]]][k]=decoder.IndiceVC[k];
		         }
            }


            numB=numB+code.rowDegree[node];


            /*******************************************************************************************************************/
            /*******************************************************************************************************************/

            } /* End of the node update*/

       if(iter == 9)
		{

				/* Compute the Bit Error Rate (BER)*/
		for (k=0; k<code.K; k++) {
			for (l=0; l<code.logGF; l++)
				BUFF[l] = table.BINGF[decide[k]][l];
			for (l=0; l<code.logGF; l++)
				if (BUFF[l] != NBIN[k][l])
					nbErrors_1 ++ ;
		}

		}

		 if(iter == 19)
		{

				/* Compute the Bit Error Rate (BER)*/
		nbErrors_2 = 0;
		for (k=0; k<code.K; k++) {
			for (l=0; l<code.logGF; l++)
				BUFF[l] = table.BINGF[decide[k]][l];
			for (l=0; l<code.logGF; l++)
				if (BUFF[l] != NBIN[k][l])
					nbErrors_2 ++ ;
		}

		}

		 if(iter == 29)
		{

				/* Compute the Bit Error Rate (BER)*/
		nbErrors_3 = 0;
		for (k=0; k<code.K; k++) {
			for (l=0; l<code.logGF; l++)
				BUFF[l] = table.BINGF[decide[k]][l];
			for (l=0; l<code.logGF; l++)
				if (BUFF[l] != NBIN[k][l])
					nbErrors_3 ++ ;
		}

		}

		 if(iter == 39)
		{

				/* Compute the Bit Error Rate (BER)*/
		nbErrors_4 = 0;
		for (k=0; k<code.K; k++) {
			for (l=0; l<code.logGF; l++)
				BUFF[l] = table.BINGF[decide[k]][l];
			for (l=0; l<code.logGF; l++)
				if (BUFF[l] != NBIN[k][l])
					nbErrors_4 ++ ;
		}

		}

		 if(iter == 49)
		{

				/* Compute the Bit Error Rate (BER)*/
		nbErrors_5 = 0;
		for (k=0; k<code.K; k++) {
			for (l=0; l<code.logGF; l++)
				BUFF[l] = table.BINGF[decide[k]][l];
			for (l=0; l<code.logGF; l++)
				if (BUFF[l] != NBIN[k][l])
					nbErrors_5 ++ ;
		}}
            ComputeAPPLog (&code, &decoder, offset);
			Decision (decide, decoder.APP, code.N, code.GF);



			synd = Syndrom (&code, decide, &table);
			if (synd == 0)
				break;
		// fin des iteration
		}

		if (nbErrors_1 > 12) {
			nbErroneousFrames_1 ++;
			if (synd == 0)
				nbUndetectedErrors ++;
		}

		if (nbErrors_2 >12 ) {
			nbErroneousFrames_2 ++;
			if (synd == 0)
				nbUndetectedErrors ++;
		}

		if (nbErrors_3 > 12) {
			nbErroneousFrames_3 ++;
			if (synd == 0)
				nbUndetectedErrors ++;
		}

		if (nbErrors_4 > 12) {
			nbErroneousFrames_4 ++;
			if (synd == 0)
				nbUndetectedErrors ++;
		}

		if (nbErrors_5 > 12) {
			nbErroneousFrames_5 ++;
			if (synd == 0)
				nbUndetectedErrors ++;
		}


		 if (nbErrors_1 > 12) {nb_erreurrrrrr_1 = nb_erreurrrrrr_1+nbErrors_1;}
		 if (nbErrors_2 > 12) {nb_erreurrrrrr_2 = nb_erreurrrrrr_2+nbErrors_2;}
		 if (nbErrors_3 > 12) {nb_erreurrrrrr_3 = nb_erreurrrrrr_3+nbErrors_3;}
		 if (nbErrors_4 > 12) {nb_erreurrrrrr_4 = nb_erreurrrrrr_4+nbErrors_4;}
		 if (nbErrors_5 > 12) {nb_erreurrrrrr_5 = nb_erreurrrrrr_5+nbErrors_5;}

        printf("\r<%d> FER= %d/ %d = %g BER= %d/ %d = %g \n",
				nbUndetectedErrors, nbErroneousFrames_1, nb,(double)nbErroneousFrames_1/nb, nb_erreurrrrrr_1 , nb,(double)nb_erreurrrrrr_1/nb/code.N/code.logGF);
fflush(stdout);
		 printf("\r<%d> FER= %d/ %d = %g BER= %d/ %d = %g \n",
				nbUndetectedErrors, nbErroneousFrames_2, nb,(double)nbErroneousFrames_2/nb, nb_erreurrrrrr_2 , nb,(double)nb_erreurrrrrr_2/nb/code.N/code.logGF);
fflush(stdout);
		 printf("\r<%d> FER= %d/ %d = %g BER= %d/ %d = %g \n",
				nbUndetectedErrors, nbErroneousFrames_3, nb,(double)nbErroneousFrames_3/nb, nb_erreurrrrrr_3 , nb,(double)nb_erreurrrrrr_3/nb/code.N/code.logGF);
fflush(stdout);
		 printf("\r<%d> FER= %d/ %d = %g BER= %d/ %d = %g \n",
				nbUndetectedErrors, nbErroneousFrames_4, nb,(double)nbErroneousFrames_4/nb, nb_erreurrrrrr_4 , nb,(double)nb_erreurrrrrr_4/nb/code.N/code.logGF);
fflush(stdout);
		 printf("\r<%d> FER= %d/ %d = %g BER= %d/ %d = %g \n",
				nbUndetectedErrors, nbErroneousFrames_5, nb,(double)nbErroneousFrames_5/nb, nb_erreurrrrrr_5 , nb,(double)nb_erreurrrrrr_5/nb/code.N/code.logGF);
fflush(stdout);

		if (nbErroneousFrames_5==100)
			break;
	}
	printf("\n");

	free(FileName); free(FileMatrix); free(name);
	free(Pe); free(decide); free(CodeWord);
	free(KSYMB); free(NSYMB); free(TMP); free(TMPindex);

	for (n=0;n<code.N;n++) free(NBIN[n]);free(NBIN);
	for (k=0;k<code.K;k++) free(KBIN[k]); free(KBIN);

	FreeCode(&code); FreeDecoder(&decoder); FreeTable(&table);
	return (EXIT_SUCCESS);
}
