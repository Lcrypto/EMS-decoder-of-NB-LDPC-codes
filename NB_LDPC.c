/*!
 * \file NB_LDPC.c
 * \brief Non-binary LDPC reduced complexity decoder with horizontal scheduling
 * \author C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015
 * \details

   Extended Min-Sum Decoder
  Horizontal Scheduling
    Layered decoder
    syndrome based architecture

 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "./include/NB_LDPC.h"





/*!
 * \fn int main(int argc, char * argv[])
 * \brief main program

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
int main(int argc, char * argv[])
{
    int 		k,l,n,iter,i,g;
    int 		**KBIN,*KSYMB,**NBIN,*NSYMB;
    int 		*decide,*CodeWord;
    int 		nb,NbIterMax;
    float 		EbN;
    int 		NbOper,NbMonteCarlo;
    float 	offset;
    char *FileName,*FileMatrix,*name;
    int synd=0, nbErrors, nbErroneousFrames = 0, nbUndetectedErrors = 0;
    int total_errors =0;







    code_t code;
    table_t table;
    decoder_t decoder;

    int node,numB ;

    int Idum=-1; // initialization of random generator
    srand(5);

    /*
     * Command line arguments
     */
    if (argc < 8)
    {
        printf("File:\n %s\n ",argv[0]);
        printf(usage_txt);
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
    printf(" Monte-Carlo simulation of Non-Binary LDPC decoder \n\n");
    printf("Simulation parameters:\n");
    printf("\n\t NbMonteCarlo     : %d", NbMonteCarlo);
    printf("\n\t NbIterMax        : %d", NbIterMax);
    printf("\n\t FileMatrix       : %s", FileMatrix);
    printf("\n\t Eb/No (dB)       : %g", EbN);
    printf("\n\t NbMax            : %d", decoder.nbMax);
    printf("\n\t Offset           : %g", offset);
    printf("\n\t NbOper           : %d\n",NbOper);

    LoadCode (FileMatrix, &code);
    LoadTables (&table, code.GF, code.logGF);
    AllocateDecoder (&code, &decoder);
    GaussianElimination (&code, &table);

// output results in a file
    FILE *opfile;
    char note[40]="FB30";
    printf("\n\t Note             : %s\n",note);
    char file_name [70];
    time_t current_time;
    char* c_time_string;


    sprintf (file_name,"./data/results_N%d_CR%0.2f_GF%d_IT%d_Offset%0.1f_nm%d_%s.txt",code.N,code.rate,code.GF,NbIterMax, offset,decoder.nbMax,note);
//sprintf (file_name,"./data/results_N%d_GF%d_IT%d_nm%d_%s.txt",code.N,code.GF,NbIterMax,decoder.nbMax,note);


    current_time = time(NULL);
    c_time_string = ctime(&current_time);
    printf("Simulation started at time: %s \n", c_time_string);


    /*
     * Memory  allocation
     */
    NBIN=(int **)calloc(code.N,sizeof(int *));
    for (n=0; n<code.N; n++)  	NBIN[n]=(int *)calloc(code.logGF,sizeof(int));
    KBIN=(int **)calloc(code.K,sizeof(int *));
    for (k=0; k<code.K; k++) 	KBIN[k]=(int *)calloc(code.logGF,sizeof(int));

    NSYMB=(int *)calloc(code.N,sizeof(int));
    KSYMB=(int *)calloc(code.K,sizeof(int));
    CodeWord=(int *)calloc(code.N,sizeof(int));

    decide=(int *)calloc(code.N,sizeof(int));

int stat_on = 20; // o if no stat , else ECN size, (nbMax)
int stat_bubble[ 3 * (code.rowDegree[0]-2) * stat_on * stat_on ];
for (i=0; i<3 * (code.rowDegree[0]-2) * stat_on * stat_on; i++)
{
    stat_bubble[i]=0;
}

    // check that dc is constant
    int dc_min=100;
    int dc_max=0;
    for (node=0; node<code.M; node++) /* Loop for the M Check nodes */
    {
        if (dc_max < code.rowDegree[node])
        {
            dc_max = code.rowDegree[node];
        }
        if (dc_min > code.rowDegree[node])
        {
            dc_min = code.rowDegree[node];
        }
    }
    if (dc_min != dc_max)
    {
        printf("d_c is not constant: dc_min= %d ; dc_max=%d !!!!!! \n", dc_min,dc_max);
    }


    int sum_it;


    int n_cv = NbOper;

////generate config_table
//    int d_1=40;
//    int d_2=15;
//    int d_3=5;
//    int **config_table;
//    int config_table_size;
//
//    config_table = build_config_table(&config_table_size,dc_max,d_1, d_2, d_3);
//    sort_config_table( config_table, config_table_size, dc_max);
//
//    config_table_size = 1000;
//
//    int j;
//    k=0;
//    l=0;
//    g=0;
//    int config_table_size_2d[dc_max];
//    for (j=0; j<dc_max; j++)
//    {
//        config_table_size_2d[j]=0;
//    }
//
//    for (i=0; i<config_table_size; i++)
//    {
////    printf(" \n \n");
////    printf("%d :\t",i);
////    for(j=0; j<dc_max; j++)
////        printf(" %d \t",config_table[i][j]);
//        if (config_table[i][0]==0) k++;
//        if (config_table[i][1]==0) g++;
//        if (config_table[i][dc_max-1]==0) l++;
//
//        for (j=0; j<dc_max; j++)
//        {
//            if (config_table[i][j]==0) config_table_size_2d[j]++;
//        }
//
//    }
//    printf("\n tablesize:%d, decorrelated size1=%d size2=%d size3=%d \n",config_table_size,k,g,l);
//    for (j=0; j<dc_max; j++)
//    {
//        printf(" %d \t",config_table_size_2d[j]);
//    }
//    getchar();





//    getchar();

    softdata_t Mvc_temp[dc_max][code.GF];
    softdata_t Mvc_temp2[dc_max][code.GF];
    softdata_t Mcv_temp[dc_max][code.GF+1];
    sum_it=0;




    for (nb=1; nb<=NbMonteCarlo; nb++)
    {
        /* Decoder re-initialization */

        /* Generate uniformly distributed information bits (KBIN)) */
        RandomBinaryGenerator (code.N, code.M, code.GF, code.logGF, KBIN, KSYMB, table.BINGF,&Idum);

        /* Encode the information bits KBIN to a (non binary) codeword NSYMB */
        Encoding (&code, &table, CodeWord, NBIN, KSYMB);

        /* Noisy channel (AWGN)*/
        ModelChannel_AWGN_BPSK (&code, &decoder, &table,  NBIN, EbN,&Idum);
        //ModelChannel_AWGN_64 (&code, &decoder, NBIN, EbN,&Idum);
        //ModelChannel(&code, &decoder,  NBIN, EbN,&Idum);
        //ModelChannel_AWGN_256QAM_4D (&code, &decoder, NBIN, EbN,&Idum);

        numB=0; /* numB is the edge number */


        /***********************************************/
        /* Implementation of the horizontal scheduling */

        /* initialisation step : init all Mcv with 0 */
        for (numB=0; numB<code.nbBranch; numB++)
        {
            for(k=0; k<code.GF; k++)
            {
                decoder.CtoV[numB][k]=0;
            }
        }
        // init APP with soft input
        for (n=0; n<code.N; n++)
        {
            for (k=0; k<code.GF; k++)
            {
                decoder.APP[n][decoder.intrinsic_GF[n][k]]=decoder.intrinsic_LLR[n][k];
            }

        }


//                    printf(" \n LLR \n",i);
//                    for (k=0; k< code.GF; k++)
//                    {
//                        printf(" %f \t %d \n ",decoder.intrinsic_LLR[0][k],decoder.intrinsic_GF[0][k] );
//
//                    }
//                getchar();
//
//
//
//                    printf(" \n APP \n",i);
//                    for (k=0; k< code.GF; k++)
//                    {
//                        printf("%d \t %f \n ",k,decoder.APP[0][k] );
//
//                    }
//                getchar();





        /* Decoding iterations*/
        for (iter=0; iter < NbIterMax - 1; iter++)
        {

//if (iter>10){stat_on=1;}else{stat_on=0;}

            numB=0;
            for (node=0; node<code.M; node++) /* Loop for the M Check nodes */
            {
                /*******************************************************************************************************************/
                /*******************************************************************************************************************/





                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for (k=0; k < code.GF; k++)
                    {

                        Mvc_temp[i][k] = decoder.APP[code.mat[node][i]][k] - decoder.CtoV[numB+i][k];
                        Mvc_temp2[i][k] = Mvc_temp[i][k];
                    }
                }


//                for (i=0; i < code.rowDegree[node]; i++)
//                {
//                    printf(" \n APP read at address %d, Mvc = SO - Mcv \n ",code.mat[node][i]);
//                    for(k=0; k < code.GF; k++)
//                    {
//                     printf(" %f -  %f = %f \n",decoder.APP[code.mat[node][i]][k],decoder.CtoV[numB+i][k],Mvc_temp[i][k]);
//
//                    }
//                    printf(" \n ");
//                }
//                    getchar();


                // sorting Mvc values
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for(k=0; k<decoder.nbMax; k++)
                    {
                        decoder.M_VtoC_LLR[i][k]=+1e5;
                        decoder.M_VtoC_GF[i][k]=0;
                        for (g=0; g<code.GF; g++)
                        {
                            if (Mvc_temp2[i][g] < decoder.M_VtoC_LLR[i][k])
                            {
                                decoder.M_VtoC_LLR[i][k]=Mvc_temp2[i][g];
                                decoder.M_VtoC_GF[i][k]=g;
                            }
                        }
                        Mvc_temp2[i][decoder.M_VtoC_GF[i][k]]=+1e5;
                    }

                    // Normalisation  */
                    for(g=1; g<decoder.nbMax; g++) decoder.M_VtoC_LLR[i][g]= decoder.M_VtoC_LLR[i][g]-decoder.M_VtoC_LLR[i][0];
                    decoder.M_VtoC_LLR[i][0]=0.0;
                }

//                for (i=0; i < code.rowDegree[node]; i++)
//                {
//                    printf(" \n sorted Mvc %d\n ",i);
//                    for(k=0; k<decoder.nbMax; k++)
//                    {
//                     printf(" %d \t %f \t %d \n",k,decoder.M_VtoC_LLR[i][k],decoder.M_VtoC_GF[i][k]);
//
//                    }
//                }
//                    getchar();


                //    syndrome_ems(node, &decoder, &code, &table, config_table, config_table_size, dc_max, offset, n_cv);// n_cv=NbOper

//syndrome_ems3(node, &decoder, &code, &table, config_table, config_table_size, dc_max, offset);
//CheckPassLogEMS_bayes (node, &decoder, &code, &table, offset);
                CheckPassLogEMS (node, &decoder, &code, &table,NbOper,offset);

                //CheckPassLogEMS_presorting (node, &decoder, &code, &table, NbOper,offset,stat_bubble,stat_on);
                //CheckPassLogEMS_presorting_synd(node, &decoder, &code, &table, NbOper,offset,stat_bubble,stat_on);
                //CheckPassLogEMS_presorting_synd_bayes(node, &decoder, &code, &table, NbOper,offset,stat_bubble,stat_on);
                //CheckPassLogEMS_tree(node, &decoder, &code, &table, offset,stat_bubble);



//                //for (i=0; i<code.rowDegree[node]; i++)
//                //{
//                i=0;
//                    printf(" Mcv%d \t node:%d \t iter=%d \n",i, node, iter);
//                    for (k=0; k<code.GF; k++)
//                    {
//                        printf("%d \t %f \t %d \n ",k, decoder.M_CtoV_LLR[i][k],decoder.M_CtoV_GF[i][k] );
//
//                    }
//                //}
//                getchar();


                // compute Mcv_temp
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for (k=0; k < code.GF; k++)
                    {
                        Mcv_temp[i][decoder.M_CtoV_GF[i][k]] = decoder.M_CtoV_LLR[i][k];
                    }
                }

//                    printf(" \n Mcv_temp \n ");
//                    for(k=0; k<code.GF; k++)
//                    {
//                     printf(" %d \t %f \n",k,Mcv_temp[0][k]);
//
//                    }
//                    getchar();



                // save in Mcv FIFO
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for (k=0; k < code.GF; k++)
                    {
                        decoder.CtoV[numB+i][k]=Mcv_temp[i][k];
                    }
                }


                //compute SO
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for (k=0; k < code.GF; k++)
                    {
                        decoder.APP[code.mat[node][i]][k] = Mcv_temp[i][k] + Mvc_temp[i][k];
                    }
                }

//                    printf(" \n SO = Mcv + Mvc \n ");
//                    for(k=0; k<code.GF; k++)
//                    {
//                     printf("  %f \t + \t %f \t = \t %f \n",Mcv_temp[0][k], Mvc_temp[0][k],decoder.APP[code.mat[node][0]][k] );
//
//                    }
//                    getchar();

                numB=numB+code.rowDegree[node];


                /*******************************************************************************************************************/
                /*******************************************************************************************************************/

            } /* End of the node update*/

            Decision (decide, decoder.APP, code.N, code.GF);
            synd = Syndrom (&code, decide, &table);
            if (synd == 0)
                break;
        }

        sum_it= sum_it+ iter +1;



        /* Compute the Bit Error Rate (BER)*/
        nbErrors = 0;
        for (k=0; k<code.K; k++)
        {
            for (l=0; l<code.logGF; l++)
                if (table.BINGF[decide[k]][l] != NBIN[k][l])
                    nbErrors ++ ;
        }

//        //fake BCH with 12errors correction capacity
 //       if (nbErrors<13) nbErrors=0;


        total_errors = total_errors + nbErrors;
        if (nbErrors != 0)
        {
            nbErroneousFrames ++;
            if (synd == 0)
                nbUndetectedErrors ++;
        }
        printf("\r<%d> FER= %d / %d = %f BER= %d / x = %f  avr_it=%.2f",
               nbUndetectedErrors, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb,total_errors,(double)total_errors/(nb*code.K*code.logGF),(double)(sum_it)/nb);
        fflush(stdout);

//                printf("\r<%d> FER=  %d / %d = %f ",nbUndetectedErrors, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb);
//                fflush(stdout);


        if (nbErroneousFrames == 40)
            break;



    }
    printf(" \n results are printed in file %s \n",file_name);




//opfile=fopen("./data/bubble_stat.txt","w");
//    for (k=0; k <(code.rowDegree[0]-2)*3; k++ )
//    {
//
//        printf( "bubble check number: %d \n",k );
//
//        for (i=0; i<decoder.nbMax; i++)
//        {
//            for (l=0; l<decoder.nbMax; l++)
//            {
//    printf(" %d ", stat_bubble[i + decoder.nbMax * l + k * decoder.nbMax * decoder.nbMax ]);
//    fprintf(opfile," %d ", stat_bubble[i + decoder.nbMax * l + k * decoder.nbMax * decoder.nbMax ]);
//            }
//            printf(" \n ");
//            fprintf(opfile," \n ");
//        }
//
//    }

//// backward stat
//opfile=fopen("./data/bubble_stat.txt","w");
//    for (k=0; k <(code.rowDegree[0]-1); k++ ) //backward
//    //for (k=0; k <2; k++ )// tree
//    {
//
//        printf( "bubble check number: %d \n",k );
//
//        for (i=0; i<stat_on; i++)
//        {
//            for (l=0; l<stat_on; l++)
//            {
//    printf(" %d ", stat_bubble[i + stat_on * l + k * stat_on * stat_on ]);
//
//    fprintf(opfile," %d ", stat_bubble[i + stat_on * l + k * stat_on * stat_on ]);
//            }
//            printf(" \n ");
//            fprintf(opfile," \n ");
//        }
//
//    }








    current_time = time(NULL);
    c_time_string = ctime(&current_time);

    opfile=fopen(file_name,"a");
//opfile=fopen("./data/results.txt","a");
    if ((opfile)==NULL)
    {
        printf(" \n !! file not found \n ");
    }
    else
    {
        fprintf(opfile," SNR:%.2f: \t FER= %d / %d = %f ", EbN, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb );
        fprintf(opfile," \t BER= %d / x = \t %f  avr_it= \t %.2f \t time: %s",total_errors, (double)total_errors/(double)(nb*code.K*code.logGF),(double)(sum_it)/nb , c_time_string );
    }
    fclose(opfile);
    printf(" \n results printed \n ");


    printf("\n");
    printf("Simulation complete at time: %s", c_time_string);
    getchar();

    free(FileName);
    free(FileMatrix);
    free(name);
    free(decide);
    free(CodeWord);
    free(KSYMB);
    free(NSYMB);

    for (n=0; n<code.N; n++) free(NBIN[n]);
    free(NBIN);
    for (k=0; k<code.K; k++) free(KBIN[k]);
    free(KBIN);

    FreeCode(&code);
    FreeDecoder(&decoder);
    FreeTable(&table);
    return (EXIT_SUCCESS);
}
