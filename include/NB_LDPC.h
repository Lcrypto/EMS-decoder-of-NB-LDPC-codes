#ifndef NB_LDPC_H_INCLUDED
#define NB_LDPC_H_INCLUDED

/*!
 * \file NB_LDPC.h
 */





#define STR_MAXSIZE 350    //la taille maximale permise pour les chaines de caract�res. Utilis�e dans la fonction calloc().


#include "struct.h"
#include "init.h"
#include "tools.h"
#include "channel.h"
#include "syndrome_decoder.h"
#include "bubble_decoder.h"



static char usage_txt[] =
"\n\
	Syntax:\n\
	\n\t NbMonteCarlo     : # simulated frames\
	\n\t NbIterMax        : # of maximum decoding iterations\
	\n\t FileMatrix       : File name of the parity-check matrix\
	\n\t EbN              : Eb/No (dB)\
	\n\t NbMax            : size of truncated messages\
	\n\t Offset           : offset correction factor (0.4 -- 1)\
	\n\t NbOper           : Maximum number of operations for sorting \n\
";






#endif // NB_LDPC_H_INCLUDED
