/*!
 * \file syndrome_decoder.c
 * \brief NB_LDPC_decoder based on syndrome
 * \author C.Marchand
 * \copyright BSD copyright
 * \date 03/03/2015
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
#include "./include/syndrome_decoder.h"


/*!
* \fn syndrome_ems
* \brief syndrome architecture as described in "A new architecture for high speed, low latency NB-LDPC check node processing"
*
*/
int syndrome_ems(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max, float offset, int n_cv)
{
    int i,j,t,k,temp;
    syndrome_type* syndrome_set;
    syndrome_set =(syndrome_type *)calloc(config_table_size,sizeof(syndrome_type));
    int dc;
    int decorrelate_size;

    GF_element *Mcv_temp;
    Mcv_temp = calloc(code->GF,sizeof(GF_element));

    GF_element decorrelated_syndrome[config_table_size];


    /// Mvc_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<decoder->nbMax; k++)
        {
            temp=decoder->M_VtoC_GF[t][k];
            decoder->M_VtoC_GF[t][k]=table->MULGF[temp][code->matValue[node][t]];
        }
    }

    int updated_gf[code->GF];
    float sat;
    float temp_mvc_llr[dc_max][code->GF];
    int temp_mvc_gf[dc_max][code->GF];

    int index_order[dc_max];
    int border = 4;
    //float array_float[code->GF];
    presorting_mvc(decoder,code,dc_max,border,index_order);
    //for (dc=0; dc<dc_max; dc++) { index_order[dc]=dc ; } // no presorting



    /// syndrome calculation //
    for (i=0 ; i< config_table_size; i++)
    {
        syndrome_set[i].LLR = 0;
        syndrome_set[i].GF =  0;
        syndrome_set[i].config = i;

        for (j=0; j<dc_max; j++)
        {
            syndrome_set[i].LLR = syndrome_set[i].LLR +  decoder->M_VtoC_LLR[j][config_table[i][j]];  // min-sum algo
            //if (syndrome_set[i].LLR < decoder->M_VtoC_LLR[j][config_table[i][j]] ) {syndrome_set[i].LLR = decoder->M_VtoC_LLR[j][config_table[i][j]];} //min-max algo

            syndrome_set[i].GF =  table->ADDGF [syndrome_set[i].GF] [ decoder->M_VtoC_GF[j][config_table[i][j]]];
        }
    }


 //   /// sorting syndrome
    sorting(syndrome_set, config_table_size);


//    printf(" after sorting \n");
//        for ( j=0 ; j < config_table_size; j++ )
//        {
//         printf("%d \t LLR:%.2f \t GF:%d  \t index:%d \n",j, syndrome_set[j].LLR, syndrome_set[j].GF, syndrome_set[j].config) ;
//        }
//    getchar();


    /// decorrelator //
    for (dc=0; dc<dc_max; dc++)
    {
        i=0;
        for (j=0; j < config_table_size ; j++)
        {

            if(config_table[syndrome_set[j].config][dc] == 0)
            {
                decorrelated_syndrome[i].LLR = syndrome_set[j].LLR;
                decorrelated_syndrome[i].GF = table->ADDGF[syndrome_set[j].GF][decoder->M_VtoC_GF[dc][0]];
                i++;
            }

        }
        //printf("dc:%d \t config_size:%d \t decorrelate_size: %d \n",dc, config_table_size,i);getchar();
        decorrelate_size=i;

//        printf(" after decorrelation \n");
//            for ( j=0 ; j < decorrelate_size; j++ )
//            {
//             printf("%d \t LLR:%.2f \t GF:%d \n",j, decorrelated_syndrome[j].LLR, decorrelated_syndrome[j].GF) ;
//            }
//        getchar();



        for(j=0; j<code->GF; j++)
        {
            updated_gf[j]=0;
        }


        /// remove redondancy, add offset and route Mcv

        // init
        for(j=0; j<code->GF; j++)
        {
            // decoder->M_CtoV_LLR[dc][j] = decorrelated_syndrome[n_cv - 1].LLR + offset;
            decoder->M_CtoV_LLR[dc][j] = 1500.0;
            decoder->M_CtoV_GF[dc][j]= j;
        }

//        max=0;
//sum=0;

        for(j=0; j<decorrelate_size; j++)
        {

            if (updated_gf[decorrelated_syndrome[j].GF] ==1)
            {
                //decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF] =decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF] *0.9375;
                decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF] = bayes(  decorrelated_syndrome[j].LLR, decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF])     ;

//if (decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF]> max )
//max = decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF];


            }
            else
            {
                decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF] = decorrelated_syndrome[j].LLR;
                updated_gf[decorrelated_syndrome[j].GF]=1;

//if (decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF]> max )
//max = decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF];
//
//sum = sum + decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF];

            }


        }


/// saturation ///
        //init
//        for(j=0; j<code->GF; j++)
//        {
//         array_float[j]=decoder->M_CtoV_LLR[dc][j] ;
//        }
//
////        //pre-saturation
//        for(j=0; j<code->GF; j++)
//        {
//         if ((array_float[j] > max) && (updated_gf[j]==1)) max=array_float[j];
//
//        }
//        for(j=0; j<code->GF; j++)
//        {
//         if (updated_gf[j]==0) array_float[j]=max;
//        }
//
//
//        for(j=0; j<code->GF; j++)
//        {
//         printf(" %f  %d \n",array_float[j],updated_gf[j]);
//
//        }
//        getchar();
//          sat = select_sort(array_float, code->GF  , n_cv-5);

        // sat=max/4;
        sat = decorrelated_syndrome[n_cv-1+ 3*dc ].LLR;
        //printf("max:%f \t sat=%f \t sat2=%f\n",max/2, sat,decorrelated_syndrome[n_cv - 1].LLR );getchar();


        for(j=0; j<code->GF; j++)
        {

            if (decoder->M_CtoV_LLR[dc][j] > sat)
                decoder->M_CtoV_LLR[dc][j] = sat + offset;


//            if (decoder->M_CtoV_LLR[dc][j] > max)
//            decoder->M_CtoV_LLR[dc][j] = max + offset;

        }


//printf (" max: %f , mcv: %f ",max,decorrelated_syndrome[n_cv - 1].LLR );
//getchar();
//
//
//i=0;
//for(j=0; j<code->GF; j++)
//{
// printf(" %d ",updated_gf[j] );
// i=i+updated_gf[j];
//}
//printf(" \n total=%d \n ", i);
//getchar();




    }// end dc loop




/// reorder mcv
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < code->GF; j++)
        {
            temp_mvc_gf[index_order[i]][j] = decoder->M_CtoV_GF[i][j];
            temp_mvc_llr[index_order[i]][j] = decoder->M_CtoV_LLR[i][j];
        }
    }
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < code->GF; j++)
        {
            decoder->M_CtoV_GF[i][j] =temp_mvc_gf[i][j];
            decoder->M_CtoV_LLR[i][j] = temp_mvc_llr[i][j];
        }
    }







    /// Mcv_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<code->GF; k++)
        {
            decoder->M_CtoV_GF[t][k]=table->DIVGF[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }
    }


//    for (i=0; i<dc_max; i++)
//    {
//            printf("decorrelation results \n");
//            printf("Mcv%d \n",i);
//            for (k=0; k<code->GF; k++)
//            {
//                printf(" LLR:%f \t GF:%d \n",decoder->M_CtoV_LLR[i][k],decoder->M_CtoV_GF[i][k]);
//            }
//            getchar();
//    }

    free(Mcv_temp);
    free(syndrome_set);
    return(0);

}




int presorting_mvc(decoder_t *decoder,const code_t *code,int dc_max,int border, int *index_order)
{

    int i,j;
    int deviation;
    float temp_llr[dc_max];
    int index_min;
    float temp_mvc_llr[dc_max][code->GF];
    int temp_mvc_gf[dc_max][code->GF];



    float min;


//// check result
//printf("mvc before sorting \n");
//for(i = 0; i < dc_max; i++)
//{
//    for(j = 0; j < decoder->nbMax; j++)
//    {
//     printf(" %f ",decoder->M_VtoC_LLR[i][j] );
//    }
//    printf(" \n");
//}
//getchar();


/// edge ordering using mvc[1]
    deviation = 1;
    for (i=0; i< dc_max; i++)
    {
        temp_llr[i]=decoder->M_VtoC_LLR[i][deviation];
    }
/// sorting and compute index order

    for(i = 0; i < dc_max; i++)
    {
        min=10000.0;
        for(j = 0; j < dc_max; j++)
        {
            if (temp_llr[j]<min)
            {
                min = temp_llr[j];
                index_min = j;
            }
        }
        temp_llr[index_min]=15000.0;
        index_order[i]=index_min;

    }
/// reorder mvc
    for(i = 0; i < dc_max; i++)
    {
        for(j = 0; j < decoder->nbMax; j++)
        {
            temp_mvc_gf[i][j] = decoder->M_VtoC_GF[index_order[i]][j];
            temp_mvc_llr[i][j] = decoder->M_VtoC_LLR[index_order[i]][j];
        }
    }
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < decoder->nbMax; j++)
        {
            decoder->M_VtoC_GF[i][j] =temp_mvc_gf[i][j];
            decoder->M_VtoC_LLR[i][j] = temp_mvc_llr[i][j];
        }

    }


//// check result
//printf("\n result after first first sort \n");
//for(i = 0; i < dc_max; i++)
//{
//    for(j = 0; j < decoder->nbMax; j++)
//    {
//     printf(" %f ",decoder->M_VtoC_LLR[i][j] );
//    }
//    printf(" \n");
//}
//getchar();




/// sort at mvc[2]
    deviation=2;
    for (i=0; i< border; i++)
    {
        temp_llr[i]=decoder->M_VtoC_LLR[i][deviation];
    }


//// find min index
//    min=10000.0;
//    for(j = 0; j < border; j++)
//    {
//       if (temp_llr[j]<min)
//        {
//            min = temp_llr[j];
//            index_min = j;
//        }
//    }

/// sorting and compute index order
    int index_order_temp[border];
    int index_order_temp2[border];
    for(i = 0; i < border; i++)
    {
        min=10000.0;
        for(j = 0; j < border; j++)
        {
            if (temp_llr[j]<min)
            {
                min = temp_llr[j];
                index_min = j;
            }
        }
        temp_llr[index_min]=15000.0;
        index_order_temp[i]=index_min;
        // printf(" \n index_min=%d \n",index_min);
    }



//// reorder mvc
//    for(j = 0; j < decoder->nbMax; j++)
//    {
//     temp_mvc_gf[0][j] = decoder->M_VtoC_GF[0][j] ;
//     temp_mvc_llr[0][j] = decoder->M_VtoC_LLR[0][j] ;
//
//     decoder->M_VtoC_GF[0][j]=decoder->M_VtoC_GF[index_min][j]   ;
//     decoder->M_VtoC_LLR[0][j]=decoder->M_VtoC_LLR[index_min][j]  ;
//
//     decoder->M_VtoC_GF[index_min][j]=temp_mvc_gf[0][j]   ;
//     decoder->M_VtoC_LLR[index_min][j]=temp_mvc_llr[0][j]  ;
//    }





/// reorder mvc
    for(i = 0; i < border; i++)
    {
        for(j = 0; j < decoder->nbMax; j++)
        {
            temp_mvc_gf[i][j] = decoder->M_VtoC_GF[index_order_temp[i]][j];
            temp_mvc_llr[i][j] = decoder->M_VtoC_LLR[index_order_temp[i]][j];
        }
    }
    for(i = 0; i < border; i++)
    {

        for(j = 0; j < decoder->nbMax; j++)
        {
            decoder->M_VtoC_GF[i][j] =temp_mvc_gf[i][j];
            decoder->M_VtoC_LLR[i][j] = temp_mvc_llr[i][j];
        }

    }


    // update index

//temp=index_order[0];
//index_order[0]=index_order[index_min];
//index_order[index_min]=temp;

//update index

    for(i = 0; i < border; i++)
    {
        index_order_temp2[i] = index_order[ index_order_temp[i]];
    }
    for(i = 0; i < border; i++)
    {
        index_order[i] = index_order_temp2[i];
    }


//// check result
//printf("\n result after second sort \n");
//for(i = 0; i < dc_max; i++)
//{
//    for(j = 1; j < decoder->nbMax; j++)
//    {
//     printf("%.2f\t",decoder->M_VtoC_LLR[i][j] );
//    }
//    printf(" \n");
//}
//getchar();
//// check index
//for(i = 0; i < dc_max; i++)
//{
//    printf(" %d \t", index_order[i]);
//}





return(0);


}




/*!
* \fn syndrome_ems
* \brief syndrome architecture as described in "A new architecture for high speed, low latency NB-LDPC check node processing"
*
*/
int syndrome_ems_median(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max)
{
    int i,j,t,k,temp;
    syndrome_type* syndrome_set;
    syndrome_set =(syndrome_type *)calloc(config_table_size,sizeof(syndrome_type));
    int dc;
    int decorrelate_size;

    GF_element *Mcv_temp;
    Mcv_temp = calloc(code->GF,sizeof(GF_element));

    GF_element decorrelated_syndrome[config_table_size];

    float offset= 0.0;
    float sat2;

//int n_cv=decoder->nbMax;



    /// Mvc_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<decoder->nbMax; k++)
        {
            temp=decoder->M_VtoC_GF[t][k];
            decoder->M_VtoC_GF[t][k]=table->MULGF[temp][code->matValue[node][t]];
        }
    }


    /// syndrome calculation //
    for (i=0 ; i< config_table_size; i++)
    {
        syndrome_set[i].LLR = 0;
        syndrome_set[i].GF =  0;
        syndrome_set[i].config = i;

        for (j=0; j<dc_max; j++)
        {
            syndrome_set[i].LLR = syndrome_set[i].LLR +  decoder->M_VtoC_LLR[j][config_table[i][j]];  // min-sum algo
            //if (syndrome_set[i].LLR < decoder->M_VtoC_LLR[j][config_table[i][j]] ) {syndrome_set[i].LLR = decoder->M_VtoC_LLR[j][config_table[i][j]];} //min-max algo

            syndrome_set[i].GF =  table->ADDGF [syndrome_set[i].GF] [ decoder->M_VtoC_GF[j][config_table[i][j]]];
        }
    }


//    /// sorting
//    sorting(syndrome_set, config_table_size);


//    printf(" after sorting \n");
//        for ( j=0 ; j < config_table_size; j++ )
//        {
//         printf("%d \t LLR:%.2f \t GF:%d  \t index:%d \n",j, syndrome_set[j].LLR, syndrome_set[j].GF, syndrome_set[j].config) ;
//        }
//    getchar();


    float temp_llr[config_table_size];

    /// decorrelator //
    for (dc=0; dc<dc_max; dc++)
    {
        i=0;
        for (j=0; j < config_table_size ; j++)
        {

            if(config_table[syndrome_set[j].config][dc] == 0)
            {
                decorrelated_syndrome[i].LLR = syndrome_set[j].LLR;
                decorrelated_syndrome[i].GF = table->ADDGF[syndrome_set[j].GF][decoder->M_VtoC_GF[dc][0]];
                temp_llr[i] = syndrome_set[j].LLR;

                i++;
            }

        }
        decorrelate_size=i;

//        printf(" after decorrelation \n");
//            for ( j=0 ; j < decorrelate_size; j++ )
//            {
//             printf("%d \t LLR:%.2f \t GF:%d \n",j, decorrelated_syndrome[j].LLR, decorrelated_syndrome[j].GF) ;
//            }
//        getchar();

//sat = decorrelated_syndrome[n_cv - 1].LLR;




        sat2= median_median_x(temp_llr ,decorrelate_size,16,6);

//printf("sat1: %f \t sat2: %f \n",sat,sat2); getchar();


        /// remove redundancy, add offset and route Mcv
        for(j=0; j<code->GF; j++)
        {
            decoder->M_CtoV_LLR[dc][j] = sat2 + offset;
            decoder->M_CtoV_GF[dc][j]= j;
        }
        for(j=0; j<decorrelate_size; j++)
        {
            if (decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF] > decorrelated_syndrome[j].LLR)
            {
                decoder->M_CtoV_LLR[dc][decorrelated_syndrome[j].GF] = decorrelated_syndrome[j].LLR;
            }

        }

    }// end dc loop




    /// Mcv_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<code->GF; k++)
        {
            decoder->M_CtoV_GF[t][k]=table->DIVGF[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }
    }


//    for (i=0; i<dc_max; i++)
//    {
//            printf("decorrelation results \n");
//            printf("Mcv%d \n",i);
//            for (k=0; k<decoder->nbMax; k++)
//            {
//                printf(" LLR:%f \t GF:%d \n",decoder->M_CtoV_LLR[i][k],decoder->M_CtoV_GF[i][k]);
//            }
//            getchar();
//    }

    free(Mcv_temp);
    free(syndrome_set);
    return(0);

}





/*!
* \fn syndrome_ems2
* \brief modification of syndrome architecture, case of a maximum of 2 deviations
*/
int syndrome_ems2(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max, float offset)
{
    int i,j,t,k,temp;
    GF_element* syndrome_set;
    syndrome_set =(GF_element *)calloc(config_table_size,sizeof(GF_element));
    int dc;

    GF_element GF_temp;
    GF_element *Mcv_temp;
    Mcv_temp = calloc(code->GF,sizeof(GF_element));

    float sat=10;
    float sat2;

int result[dc_max];


    /// Mvc_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<decoder->nbMax; k++)
        {
            temp=decoder->M_VtoC_GF[t][k];
            decoder->M_VtoC_GF[t][k]=table->MULGF[temp][code->matValue[node][t]];
        }
    }

    float temp_mvc_llr[dc_max][code->GF];
    int temp_mvc_gf[dc_max][code->GF];
    int index_order[dc_max];
    int border = 6;
    presorting_mvc(decoder,code,dc_max,border,index_order);
    //for (dc=0; dc<dc_max; dc++) { index_order[dc]=dc ; } // no presorting



    /// syndrome calculation //
    for (i=0 ; i< config_table_size; i++)
    {
        syndrome_set[i].LLR = 0;
        syndrome_set[i].GF =  0;

        for (j=0; j<dc_max; j++)
        {
            syndrome_set[i].LLR = syndrome_set[i].LLR +  decoder->M_VtoC_LLR[j][config_table[i][j]];  // min-sum algo
            //if (syndrome_set[i].LLR < decoder->M_VtoC_LLR[j][config_table[i][j]] ) {syndrome_set[i].LLR = decoder->M_VtoC_LLR[j][config_table[i][j]];} //min-max algo

            syndrome_set[i].GF =  table->ADDGF [syndrome_set[i].GF] [ decoder->M_VtoC_GF[j][config_table[i][j]]];
        }
    }


/// remove GF redundancy by keeping the 3 minimum and index, storing in GF order
    GF_syndrome_type GF_syndrome[code->GF];

// init
    for ( i=0 ; i < code->GF; i++ )
    {
        GF_syndrome[i].min1=sat;
        GF_syndrome[i].min2=sat;
        GF_syndrome[i].min3=sat;
        GF_syndrome[i].index1=0;
        GF_syndrome[i].index2=0;
        GF_syndrome[i].index3=0;

    }

    GF_syndrome_type synd_temp;

    for (i=0 ; i< config_table_size; i++)
    {
        GF_temp = syndrome_set[i];
        synd_temp = GF_syndrome[GF_temp.GF];
        if (synd_temp.min1 > GF_temp.LLR)
        {
//            synd_temp.min3 = synd_temp.min2;
//            synd_temp.index3 = synd_temp.index2;
            synd_temp.min2 = synd_temp.min1;
            synd_temp.index2 = synd_temp.index1;
            synd_temp.min1 = GF_temp.LLR;
            synd_temp.index1 = i;

        }
        else if(synd_temp.min2 > GF_temp.LLR )
        {
//            synd_temp.min3 = synd_temp.min2;
//            synd_temp.index3 = synd_temp.index2;
            synd_temp.min2 = GF_temp.LLR;
            synd_temp.index2 = i;
        }
//        else if(synd_temp.min3 > GF_temp.LLR )
//        {
//            synd_temp.min3 = synd_temp.min2;
//            synd_temp.index3 = synd_temp.index2;
//        }

        GF_syndrome[GF_temp.GF]=synd_temp;
    }


//    for ( i=0 ; i < code->GF; i++ )
//    {
//     printf("%d \t min1:%.2f \t min2:%.2f \t min3:%.2f \t index1:%d \t index2:%d \n ",i, GF_syndrome[i].min1, GF_syndrome[i].min2, GF_syndrome[i].min3, GF_syndrome[i].index1, GF_syndrome[i].index2 ) ;
//    }
//getchar();

    //second step
    for (i=0 ; i< config_table_size; i++)
    {
        GF_temp = syndrome_set[i];
        synd_temp = GF_syndrome[GF_temp.GF];


        if(  (synd_temp.min3 > GF_temp.LLR)  &&  (synd_temp.min2 < GF_temp.LLR) )
        {
            temp = check_deviation(config_table,dc_max, synd_temp.index1, synd_temp.index2,i, result);
            if (temp == 0)
            {
            synd_temp.min3 = GF_temp.LLR;
            synd_temp.index3 =i;


//                printf(" \n index1:%d  LLR=%f \n", synd_temp.index1, synd_temp.min1);
//                for (k=0;k<dc_max;k++)
//                {
//                    printf(" %d ", config_table[synd_temp.index1][k]);
//                }
//                printf(" \n index2:%d  LLR=%f \n", synd_temp.index2, synd_temp.min2);
//                for (k=0;k<dc_max;k++)
//                {
//                    printf(" %d ", config_table[synd_temp.index2][k]);
//                }
//                printf(" \n  index3:%d LLR3=%f \n",synd_temp.index3, synd_temp.min3);
//                for (k=0;k<dc_max;k++)
//                {
//                    printf(" %d ", config_table[synd_temp.index3][k]);
//                }
//
//
//getchar();

            }
        }

        GF_syndrome[GF_temp.GF]=synd_temp;
    }




//    for ( i=0 ; i < code->GF; i++ )
//    {
//     printf("%d \t min1:%.2f \t min2:%.2f \t min3:%.2f \t index1:%d \t index2:%d  \t index3:%d\n ",i, GF_syndrome[i].min1, GF_syndrome[i].min2, GF_syndrome[i].min3, GF_syndrome[i].index1, GF_syndrome[i].index2,GF_syndrome[i].index3 ) ;
//    }
//getchar();






     /// decorrelator 3min//
     int dev1,dev2;
     double min1,min2,min3;
    for (dc=0; dc<dc_max; dc++)
    {
        for (j=0; j < code->GF ; j++)
        {

            dev1=config_table[GF_syndrome[j].index1][dc];
            dev2=config_table[GF_syndrome[j].index2][dc];
            min1= GF_syndrome[j].min1;
            min2= GF_syndrome[j].min2;
            min3= GF_syndrome[j].min3;



//printf("llr1:%f  llr2:%f  llr3:%f \n ",GF_syndrome[j].min1,GF_syndrome[j].min2,GF_syndrome[j].min3 );

//            // normal
//            if( dev1 == 0)
//            {
//                Mcv_temp[j].LLR = min1;
//            }
//            else if ( dev2 == 0)
//            {
//                Mcv_temp[j].LLR = min2;
//            }
//            else
//            {
//                Mcv_temp[j].LLR =min3;
//            }
//            //printf("llr1:%f ",Mcv_temp[j].LLR );



                // using bayes
                if (min1==sat)
                {
                    Mcv_temp[j].LLR =sat;
                }
            else if( (dev1 == 0) && (dev2 == 0))
            {
                Mcv_temp[j].LLR = bayes(min2 , min1  );
            }
            else if ( dev1 == 0)
            {
                Mcv_temp[j].LLR = bayes(min3,min1);
            }
            else if (dev2 == 0)
            {
                Mcv_temp[j].LLR = bayes(min3,min2);
            }
            else
            {
                Mcv_temp[j].LLR =min3;
            }

//printf("llr2:%f \n",Mcv_temp[j].LLR );
//getchar();
        }







//printf("mcv_temp \n");
//    for ( j=0 ; j < code->GF; j++ )
//    {
//     printf("%d \t LLR:%.2f \t GF:%d \n",j, Mcv_temp[j].LLR, Mcv_temp[j].GF) ;
//    }
//getchar();


        for(j=0; j<code->GF; j++)
        {
            decoder->M_CtoV_LLR[dc][j] = Mcv_temp[j].LLR;
            decoder->M_CtoV_GF[dc][j]= table->ADDGF [j] [ decoder->M_VtoC_GF[dc][0]];;
        }


//sorting(Mcv_temp, code->GF);
//     sat2=Mcv_temp[18].LLR;

//     printf(" \n sorted Mcv \n");
//    for ( j=0 ; j < code->GF; j++ )
//    {
//     printf("%d \t LLR:%.2f \t GF:%d \n",j, Mcv_temp[j].LLR, Mcv_temp[j].GF) ;
//    }
//     printf(" \n saturation:%.2f \n", sat2); getchar();



        float array_test[code->GF];
        for (i=0; i<code->GF; i++)
        {
            array_test[i]=Mcv_temp[i].LLR;
        }



        //sat2= median_median256(array_test);
        sat2= median_median64(array_test);



////for (i=0; i<64; i++)
////{
////    printf(" %.2f  ",array_test[i]);
////}
//printf(" \n sat2=%.2f \t median=%.2f ",sat2, median);
//getchar();




        //offset and sat
        for(k=0; k<code->GF; k++)
        {
            if (decoder->M_CtoV_LLR[dc][k] > sat2) decoder->M_CtoV_LLR[dc][k]= sat2 + offset;

        }
//printf("Mcv: \n");
//    for ( j=0 ; j < code->GF; j++ )
//    {
//     printf("%d \t LLR:%.2f \n",j, decoder->M_CtoV_LLR[dc][j]) ;
//    }
//getchar();


    }// end dc loop



/// reorder mcv
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < code->GF; j++)
        {
            temp_mvc_gf[index_order[i]][j] = decoder->M_CtoV_GF[i][j];
            temp_mvc_llr[index_order[i]][j] = decoder->M_CtoV_LLR[i][j];
        }
    }
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < code->GF; j++)
        {
            decoder->M_CtoV_GF[i][j] =temp_mvc_gf[i][j];
            decoder->M_CtoV_LLR[i][j] = temp_mvc_llr[i][j];
        }
    }






    /// Mcv_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<code->GF; k++)
        {
            decoder->M_CtoV_GF[t][k]=table->DIVGF[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }
    }


//    for (i=0; i<dc_max; i++)
//    {
//            printf("decorrelation results \n");
//            printf("Mcv%d \n",i);
//            for (k=0; k<decoder->nbMax; k++)
//            {
//                printf(" LLR:%f \t GF:%d \n",decoder->M_CtoV_LLR[i][k],decoder->M_CtoV_GF[i][k]);
//            }
//            getchar();
//    }

    free(Mcv_temp);
    free(syndrome_set);
    return(0);

}





/*!
* \fn syndrome_ems3
* \brief modification of syndrome architecture, case of a maximum of 3 deviations
*/
int syndrome_ems3(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max, float offset)
{
    int i,j,t,k,temp;
    GF_element* syndrome_set;
    syndrome_set =(GF_element *)calloc(config_table_size,sizeof(GF_element));
    int dc;

    GF_element GF_temp;
    GF_element *Mcv_temp;
    Mcv_temp = calloc(code->GF,sizeof(GF_element));

    float sat=20;
    float sat2;
    //float offset= 0.5;






    /// Mvc_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<decoder->nbMax; k++)
        {
            temp=decoder->M_VtoC_GF[t][k];
            decoder->M_VtoC_GF[t][k]=table->MULGF[temp][code->matValue[node][t]];
        }
    }


    float temp_mvc_llr[dc_max][code->GF];
    int temp_mvc_gf[dc_max][code->GF];
    int index_order[dc_max];
    int border = 6;
    presorting_mvc(decoder,code,dc_max,border,index_order);
    //for (dc=0; dc<dc_max; dc++) { index_order[dc]=dc ; } // no presorting




    /// syndrome calculation //
    for (i=0 ; i< config_table_size; i++)
    {
        syndrome_set[i].LLR = 0;
        syndrome_set[i].GF =  0;

        for (j=0; j<dc_max; j++)
        {
            syndrome_set[i].LLR = syndrome_set[i].LLR +  decoder->M_VtoC_LLR[j][config_table[i][j]];  // min-sum algo
            //if (syndrome_set[i].LLR < decoder->M_VtoC_LLR[j][config_table[i][j]] ) {syndrome_set[i].LLR = decoder->M_VtoC_LLR[j][config_table[i][j]];} //min-max algo

            syndrome_set[i].GF =  table->ADDGF [syndrome_set[i].GF] [ decoder->M_VtoC_GF[j][config_table[i][j]]];
        }
    }


/// remove GF redundancy by keeping the 3 minimum and index, storing in GF order
    GF_syndrome_type GF_syndrome[code->GF];

// init
    for ( i=0 ; i < code->GF; i++ )
    {
        GF_syndrome[i].min1=sat;
        GF_syndrome[i].min2=sat;
        GF_syndrome[i].min3=sat;
        GF_syndrome[i].min4=sat;
        GF_syndrome[i].index1=0;
        GF_syndrome[i].index2=0;
        GF_syndrome[i].index2=0;

    }

    GF_syndrome_type synd_temp;

    for (i=0 ; i< config_table_size; i++)
    {
        GF_temp = syndrome_set[i];
        synd_temp = GF_syndrome[GF_temp.GF];
        if (synd_temp.min1 > GF_temp.LLR)
        {
            synd_temp.min4 = synd_temp.min3;
            synd_temp.min3 = synd_temp.min2;
            synd_temp.index3 = synd_temp.index2;
            synd_temp.min2 = synd_temp.min1;
            synd_temp.index2 = synd_temp.index1;
            synd_temp.min1 = GF_temp.LLR;
            synd_temp.index1 = i;

        }
        else if(synd_temp.min2 > GF_temp.LLR )
        {
            synd_temp.min4 = synd_temp.min3;
            synd_temp.min3 = synd_temp.min2;
            synd_temp.index3 = synd_temp.index2;
            synd_temp.min3 = synd_temp.min2;
            synd_temp.min2 = GF_temp.LLR;
            synd_temp.index2 = i;
        }
        else if(synd_temp.min3 > GF_temp.LLR )
        {
            synd_temp.min4 = synd_temp.min3;
            synd_temp.min3 = GF_temp.LLR;
            synd_temp.index3 = i;
        }
        else if(synd_temp.min4 > GF_temp.LLR )
        {
            synd_temp.min4 = GF_temp.LLR;
        }

        GF_syndrome[GF_temp.GF]=synd_temp;

    }



//    for ( i=0 ; i < code->GF; i++ )
//    {
//     printf("%d \t min1:%.2f \t min2:%.2f \t min3:%.2f \t index1:%d \t index2:%d \n ",i, GF_syndrome[i].min1, GF_syndrome[i].min2, GF_syndrome[i].min3, GF_syndrome[i].index1, GF_syndrome[i].index2 ) ;
//    }
//getchar();



    /// decorrelator //
    for (dc=0; dc<dc_max; dc++)
    {

        // decorrelate

        for (j=0; j < code->GF ; j++)
        {

            if(config_table[GF_syndrome[j].index1][dc] == 0)
            {
                Mcv_temp[j].LLR = GF_syndrome[j].min1;
            }
            else if (config_table[GF_syndrome[j].index2][dc] == 0)
            {
                Mcv_temp[j].LLR = GF_syndrome[j].min2;
            }
            else if (config_table[GF_syndrome[j].index3][dc] == 0)
            {
                Mcv_temp[j].LLR = GF_syndrome[j].min3;
            }
            else
            {
                Mcv_temp[j].LLR = GF_syndrome[j].min4;
//                printf(" \n min4 is required! \n");
//                printf("dc=%d GF=%d \n ", dc, j);
//                printf(" \n index1:%d  LLR=%f \n", GF_syndrome[j].index1, GF_syndrome[j].min1);
//                for (k=0;k<dc_max;k++)
//                {
//                    printf(" %d ", config_table[GF_syndrome[j].index1][k]);
//                }
//                printf(" \n index2:%d  LLR=%f \n", GF_syndrome[j].index2, GF_syndrome[j].min2);
//                for (k=0;k<dc_max;k++)
//                {
//                    printf(" %d ", config_table[GF_syndrome[j].index2][k]);
//                }
//                printf(" \n index3:%d  LLR=%f \n", GF_syndrome[j].index3, GF_syndrome[j].min3);
//                for (k=0;k<dc_max;k++)
//                {
//                    printf(" %d ", config_table[GF_syndrome[j].index3][k]);
//                }
//getchar();




            }

        }

//printf("mcv_temp \n");
//    for ( j=0 ; j < code->GF; j++ )
//    {
//     printf("%d \t LLR:%.2f \t GF:%d \n",j, Mcv_temp[j].LLR, Mcv_temp[j].GF) ;
//    }
//getchar();


        for(j=0; j<code->GF; j++)
        {
            decoder->M_CtoV_LLR[dc][j] = Mcv_temp[j].LLR;
            decoder->M_CtoV_GF[dc][j]= table->ADDGF [j] [ decoder->M_VtoC_GF[dc][0]];;
        }


//sorting(Mcv_temp, code->GF);
//     sat2=Mcv_temp[18].LLR;

//     printf(" \n sorted Mcv \n");
//    for ( j=0 ; j < code->GF; j++ )
//    {
//     printf("%d \t LLR:%.2f \t GF:%d \n",j, Mcv_temp[j].LLR, Mcv_temp[j].GF) ;
//    }
//     printf(" \n saturation:%.2f \n", sat2); getchar();



        float array_test[code->GF];
        for (i=0; i<code->GF; i++)
        {
            array_test[i]=Mcv_temp[i].LLR;
        }



        //sat2= median_median256(array_test);
sat2= median_median64(array_test);



////for (i=0; i<64; i++)
////{
////    printf(" %.2f  ",array_test[i]);
////}
//printf(" \n sat2=%.2f \t median=%.2f ",sat2, median);
//getchar();




        //offset and sat
        for(k=0; k<code->GF; k++)
        {
            if (decoder->M_CtoV_LLR[dc][k] > sat2) decoder->M_CtoV_LLR[dc][k]= sat2 + offset;

        }
//printf("Mcv: \n");
//    for ( j=0 ; j < code->GF; j++ )
//    {
//     printf("%d \t LLR:%.2f \n",j, decoder->M_CtoV_LLR[dc][j]) ;
//    }
//getchar();


    }// end dc loop


/// reorder mcv
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < code->GF; j++)
        {
            temp_mvc_gf[index_order[i]][j] = decoder->M_CtoV_GF[i][j];
            temp_mvc_llr[index_order[i]][j] = decoder->M_CtoV_LLR[i][j];
        }
    }
    for(i = 0; i < dc_max; i++)
    {

        for(j = 0; j < code->GF; j++)
        {
            decoder->M_CtoV_GF[i][j] =temp_mvc_gf[i][j];
            decoder->M_CtoV_LLR[i][j] = temp_mvc_llr[i][j];
        }
    }



    /// Mcv_rotation //
    for(t=0; t<dc_max; t++)
    {
        for(k=0; k<code->GF; k++)
        {
            decoder->M_CtoV_GF[t][k]=table->DIVGF[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }
    }


//    for (i=0; i<dc_max; i++)
//    {
//            printf("decorrelation results \n");
//            printf("Mcv%d \n",i);
//            for (k=0; k<decoder->nbMax; k++)
//            {
//                printf(" LLR:%f \t GF:%d \n",decoder->M_CtoV_LLR[i][k],decoder->M_CtoV_GF[i][k]);
//            }
//            getchar();
//    }

    free(Mcv_temp);
    free(syndrome_set);
    return(0);

}






/*!
* \fn sorting
* \brief sorting function
*
*/
int sorting(syndrome_type* syndrome_set, int syndrome_table_size)
{
    int i,j;
    syndrome_type syndrome_temp;

    // insertion sort
    for(j = 1; j < syndrome_table_size; j++)    // Start with 1 (not 0)
    {
        syndrome_temp = syndrome_set[j];
        i=j-1;
        while( (i >=0) && ( syndrome_set[i].LLR > syndrome_temp.LLR ) )
        {
            syndrome_set[i+1] = syndrome_set[i]; //Higher values move up
            i=i-1;

        }
        syndrome_set[i+1] = syndrome_temp;    //Put key into its proper location
    }
    return(0);
}

/*!
* \fn select_sort
* \brief selection algorithm: return the order smallest element
*
*/
float select_sort(float* array_float,int size_array, int order)
{
    int i,j;
    float temp[size_array];
    float sorted_array[size_array];
    float min;
    int index_min;

    for(i=0; i<size_array; i++)
    {
        temp[i]=array_float[i];
    }

    // sorting order first values
    for(i=0; i<size_array; i++)
    {
        min=10000.0;
        for(j=0; j<size_array; j++)
        {
            if (min > temp[j])
            {
                min=temp[j];
                index_min=j;
            }

        }
        temp[index_min]=15000.0;
        sorted_array[i]=min;
    }
    return(sorted_array[order]);

}

/*!
* \fn median_median64
* \brief compute median of median on 64 values
*
*/
float median_median64(float* array_float)
{
    // test  median of 8 median of aray size 8
    int i,j;
    float medians[8];
    float temp[8];
    float median_median;

    for( i=0 ; i<8; i++)
    {
        for (j=0; j<8; j++)
        {
            temp[j]=array_float[i*8+j];
        }
        medians[i] = select_sort(temp,8,4);

    }
    median_median= select_sort(medians,8,4);


//// print kth order statistic
//j=0;
//for (i=0 ; i<64; i++)
//{
//    if (median_median > array_float[i])
//    {
//        j++;
//    }
//}
//printf("\n order is: %d \n",j);
//getchar();


    return(median_median);

}



/*!
* \fn median_median256
* \brief compute median of median on 256 values
*
*/
float median_median256(float* array_float)
{
    // test  median of 16 median of aray size 16
    int i,j;
    float medians[16];
    float temp[16];
    float median_median;

    for( i=0 ; i<16; i++)
    {
        for (j=0; j<16; j++)
        {
            temp[j]=array_float[i*16+j];
        }
        medians[i] = select_sort(temp,16,4);

    }
    median_median= select_sort(medians,16,4);


//// print kth order statistic
//j=0;
//for (i=0 ; i<256; i++)
//{
//    if (median_median > array_float[i])
//    {
//        j++;
//    }
//}
//printf("\n order is: %d \n",j);
//getchar();


    return(median_median);

}



/*!
* \fn median_median_x
* \brief compute median of median on x values
*
*/
float median_median_x(float* array_float, int array_float_size, int median_size, int median_select)
{
    // test  median of 16 median of aray size 16
    int i,j;
    int nb_medians;
    nb_medians=array_float_size/median_size;
    float medians[nb_medians];
    float temp[median_size];
    float median_median;

    for( i=0 ; i<nb_medians; i++)
    {
        for (j=0; j<median_size; j++)
        {
            temp[j]=array_float[i*median_size+j];
        }
        medians[i] = select_sort(temp,median_size,median_select);

    }
    median_median= select_sort(medians,nb_medians,median_select);


// print kth order statistic
    j=0;
    for (i=0 ; i<array_float_size; i++)
    {
        if (median_median > array_float[i])
        {
            j++;
        }
    }
    printf("\n order is: %d \n",j);
    getchar();


    return(median_median);

}







int** AllocateArray(int line, int column)
{
    int i;
    int** pArray = (int**)calloc(line,sizeof(int*));
    for ( i = 0; i < line; i++ )
    {
        pArray[i] = (int*)calloc(column,sizeof(int));
    }
    return pArray;
}

int compute_config_table_size(int dc_max, int d_1, int d_2, int d_3)
{

    int config_table_size;
    int card_d_1, card_d_2, card_d_3;

    card_d_1=dc_max*d_1;
    card_d_2=combin(dc_max,2)*d_2*d_2;
    card_d_3=combin(dc_max,3)*d_3*d_3*d_3;

    config_table_size = 1 + card_d_1 +card_d_2 + card_d_3;

    printf( "computed config size:%d  (%d + %d + %d + 1) ",config_table_size, card_d_1,card_d_2,card_d_3);

    return(config_table_size);


}

int** build_config_table(int* config_table_size_p,int dc_max,int d_1,int d_2, int d_3)
{
    int i;
    int config_table_size_d=0;
    int** config_table = (int**)calloc(dc_max,sizeof(int*));
    //int config_table_size;
    *config_table_size_p = compute_config_table_size(dc_max, d_1, d_2, d_3);
    //*config_table_size_p = 1000;
    config_table = AllocateArray(*config_table_size_p, dc_max);
//   *config_table_size_p = gen_config_table(config_table,dc_max,d_1, d_2, d_3);
    *config_table_size_p = gen_config_table2(config_table,dc_max,d_1,d_2,d_3);
 //   *config_table_size_p = gen_config_table3(config_table,dc_max,d_1,d_2,d_3);
 //   *config_table_size_p = gen_config_table4(config_table,dc_max,d_1,d_2,d_3);

    // count number of config after decorrelation
    for( i=0; i< *config_table_size_p ; i++)
    {
        if (config_table[i][0] == 0) config_table_size_d++;
    }


//    for( i=0; i< *config_table_size_p ; i++)
//    {
//        for(j=0;j<dc_max;j++)
//        {printf(" %d ",config_table[i][j]);}
//        printf(" \t \t %d \n ",i);
//    }
    printf("d_1=%d, d_2=%d, d_3=%d, config table size: %d , size afer decorrelation: %d \n",d_1, d_2, d_3, *config_table_size_p, config_table_size_d);
    //getchar();


    return (config_table);
}


/*!
 * \fn gen_config_table
 * \brief generate configuration table
 * Inputs
 * 	- d_c, d_1, d_2, d_3
 * Outputs
 * 	- config table
 *  - size of config table
 *
 * \author C. Marchand
 */
int gen_config_table(int ** config_table,int dc_max,int d_1,int d_2, int d_3)
{
    int config;
    int i,j,k,l,m,n;

    config=1;

    // d_1
    for (i=0; i<dc_max; i++)
    {
        for (j=0; j<d_1 ; j++)
        {
            config_table[config][i]=j+1;
            config++;
        }
    }

    //d_2
    for (i=0; i<dc_max-1; i++)
    {
        for (j=i+1; j<dc_max ; j++)
        {
            for (k=0; k<d_2; k++ )
            {
                for (l=0; l<d_2; l++)
                {
                    config_table[config][i]=k+1;
                    config_table[config][j]=l+1;
                    config++;
                }
            }
        }
    }

    //d_3
    for (i=0; i<dc_max-2; i++)
    {
        for (j=i+1; j<dc_max-1 ; j++)
        {
            for (k=j+1; k<dc_max; k++ )
            {
                for (l=0; l<d_3; l++)
                {
                    for (m=0; m<d_3; m++)
                    {
                        for (n=0; n<d_3; n++)
                        {
                            config_table[config][i]=l+1;
                            config_table[config][j]=m+1;
                            config_table[config][k]=n+1;
                            config++;
                        }
                    }
                }
            }
        }
    }


    return (config);
}

/*!
 * \fn gen_config_table2
 * \brief generate configuration table with trapeze shape
 * Inputs
 * 	- d_c, d_1, d_2, d_3
 * Outputs
 * 	- config table
 *  - size of config table
 *
 * \author C. Marchand
 */
int gen_config_table2(int ** config_table,int dc_max,int d_1,int d_2, int d_3)
{

    int config;
    int i,j,k,l,m,n,o,p;
    int d_4=2;

    config=1;

    // d_1
    for (i=0; i<dc_max; i++)
    {
        for (j=0; j<d_1 ; j++)
        {
            config_table[config][i]=j+1;
            config++;
        }
    }

    //d_2
    for (i=0; i<dc_max-1; i++)
    {
        for (j=i+1; j<dc_max ; j++)
        {
            for (k=0; k<d_2; k++ )
            {
                for (l=0; l<d_2; l++)
                {
                    if ( k+l < d_2)
                    {
                        config_table[config][i]=k+1;
                        config_table[config][j]=l+1;
                        config++;
                    }

                }
            }
        }
    }

    //d_3
    for (i=0; i<dc_max-2; i++)
    {
        for (j=i+1; j<dc_max-1 ; j++)
        {
            for (k=j+1; k<dc_max; k++ )
            {
                for (l=0; l<d_3; l++)
                {
                    for (m=0; m<d_3; m++)
                    {
                        for (n=0; n<d_3; n++)
                        {
                            if (l+m+n<d_3)
                            {
                                config_table[config][i]=l+1;
                                config_table[config][j]=m+1;
                                config_table[config][k]=n+1;
                                config++;
                            }
                        }
                    }
                }
            }
        }
    }


//    //d_4
for (o=0; o<dc_max-3; o++)
{
    for (i=o+1; i<dc_max-2; i++)
    {
        for (j=i+1; j<dc_max-1 ; j++)
        {
            for (k=j+1; k<dc_max; k++ )
            {
                for (l=0; l<d_4; l++)
                {
                    for (m=0; m<d_4; m++)
                    {
                        for (n=0; n<d_4; n++)
                        {
                            for (p=0; p<d_4; p++)
                            {
                                if (l+m+n<d_4)
                                {
                                    config_table[config][i]=l+1;
                                    config_table[config][j]=m+1;
                                    config_table[config][k]=n+1;
                                    config_table[config][o]=p+1;
                                    config++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}




    return (config);
}





/*!
 * \fn gen_config_table3
 * \brief generate configuration table with trapeze shape
 * Inputs
 * 	- d_c, d_1, d_2, d_3
 * Outputs
 * 	- config table
 *  - size of config table
 *
 * \author C. Marchand
 */
int gen_config_table3(int ** config_table,int d_c,int d_1)
{
    int config;
    int i,j;

    config=1;

    // d_1
    for (i=0; i<d_c; i++)
    {
        for (j=0; j<d_1 ; j++)
        {
            config_table[config][i]=j+1;
            config++;
        }
    }

    //d_2
    for (i=0; i<d_c-1; i++)
    {
        for (j=i+1; j<d_c ; j++)
        {
            config_table[config][i]=1;
            config_table[config][j]=1;
            config++;
            config_table[config][i]=1;
            config_table[config][j]=2;
            config++;
            config_table[config][i]=2;
            config_table[config][j]=1;
            config++;

        }
    }



    return (config);
}




/*!
 * \fn gen_config_table2
 * \brief generate configuration table with trapeze shape, not regular
 * Inputs
 * 	- d_c, d_1, d_2, d_3
 * Outputs
 * 	- config table
 *  - size of config table
 *
 * \author C. Marchand
 */
int gen_config_table4(int ** config_table,int dc_max,int d_1,int d_2, int d_3)
{

    int config;
    int i,j,k,l;

    config=1;
    int border=6;
    int border0=dc_max-3;

    // d_1 part1
    for (i=0; i<border; i++)
    {
        for (j=0; j<d_1 ; j++)
        {
            config_table[config][i]=j+1;
            config++;
        }
    }
    //d_1 part2
    for (i=border; i<border0; i++)
    {
        for (j=0; j<d_2 ; j++)
        {
            config_table[config][i]=j+1;
            config++;
        }
    }
    //d_1 part3
    for (i=border0; i<dc_max; i++)
    {
            config_table[config][i]=1;
            config++;
    }



    //d_2 part1
    for (i=0; i<border-1; i++)
    {
        for (j=i+1; j<border ; j++)
        {
            for (k=0; k<d_2; k++ )
            {
                for (l=0; l<d_2; l++)
                {
                    if ( k+l < d_2)
                    {
                        config_table[config][i]=k+1;
                        config_table[config][j]=l+1;
                        config++;
                    }

                }
            }
        }
    }
    //d_2 part2
    for (i=border0-1; i > border-1 ; i--)
    {
        for (j=border-1; j>-1 ; j--)
        {
            config_table[config][i]=1;
            config_table[config][j]=1;
            config++;
        }
    }
    //d_2 part3
    for (i=border; i < border0 ; i++)
    {
        config_table[config][0]=2;
        config_table[config][i]=1;
        config++;
    }


////d_2 part4
//for (i=1; i<border; i++)
//{
//    for (j=d_2; j<d_1; j++)
//    {
//        config_table[config][0]=2;
//        config_table[config][i]=j;
//        config++;
//    }
//}



//    //d_2 part4
//    for (i=border; i < dc_max ; i++)
//    {
//                        config_table[config][1]=2;
//                        config_table[config][i]=1;
//                        config++;
//    }

//    //d_2 part5
//    for (i=border; i < dc_max ; i++)
//    {
//                        config_table[config][2]=2;
//                        config_table[config][i]=1;
//                        config++;
//    }






//    //d_3 part1 trapeze shape
//    for (i=0; i<border-2; i++)
//    {
//        for (j=i+1; j<border-1 ; j++)
//        {
//            for (m=j+1; m<border ; m++)
//            {
//                for (k=0; k<d_2; k++ )
//                {
//                    for (l=0; l<d_2; l++)
//                    {
//                        for (n=0; n<d_2; n++)
//                        {
//                            if ( k+l+n < d_2)
//                            {
//                                config_table[config][i]=k+1;
//                                config_table[config][j]=l+1;
//                                config_table[config][m]=n+1;
//                                config++;
//                            }
//                        }
//
//                    }
//                }
//            }
//        }
//    }

    //d_3 part1
    for (i=0; i<border-2; i++)
    {
        for (j=i+1; j<border-1 ; j++)
        {
            for (k=j+1; k<border; k++ )
            {

                config_table[config][i]=1;
                config_table[config][j]=1;
                config_table[config][k]=1;
                config++;
            }
        }
    }
    //d_3 part2
    for (j=1; j<border-1 ; j++)
    {
        for (k=j+1; k<border; k++ )
        {

            config_table[config][0]=2;
            config_table[config][j]=1;
            config_table[config][k]=1;
            config++;
        }
    }
//            //d_3 part3
//    for (j=1; j<border-1 ; j++)
//        {
//            for (k=j+1; k<border; k++ )
//            {
//
//                config_table[config][0]=3;
//                config_table[config][j]=1;
//                config_table[config][k]=1;
//                config++;
//            }
//        }


    //d_4
    for (l=0; l<border-3; l++)
    {
        for (i=l+1; i<border-2; i++)
        {
            for (j=i+1; j<border-1 ; j++)
            {
                for (k=j+1; k<border; k++ )
                {
                    config_table[config][i]=1;
                    config_table[config][j]=1;
                    config_table[config][k]=1;
                    config_table[config][l]=1;
                    config++;
                }
            }
        }
    }
    //d_4 part 2
    for (i=1; i<border-2; i++)
    {
        for (j=i+1; j<border-1 ; j++)
        {
            for (k=j+1; k<border; k++ )
            {

                config_table[config][0]=2;
                config_table[config][i]=1;
                config_table[config][j]=1;
                config_table[config][k]=1;
                config++;
            }
        }
    }

////    //d_5
//
//    for (n=0; n<border-4; n++)
//    {
//        for (l=n+1; l<border-3; l++)
//        {
//            for (i=l+1; i<border-2; i++)
//            {
//                for (j=i+1; j<border-1 ; j++)
//                {
//                    for (k=j+1; k<border; k++ )
//                    {
//
//                        config_table[config][i]=1;
//                        config_table[config][j]=1;
//                        config_table[config][k]=1;
//                        config_table[config][l]=1;
//                        config_table[config][n]=1;
//                        config++;
//                    }
//                }
//            }
//        }
//    }
////    //d_5 part2
//    for (l=1; l<border-3; l++)
//    {
//        for (i=l+1; i<border-2; i++)
//        {
//            for (j=i+1; j<border-1 ; j++)
//            {
//                for (k=j+1; k<border; k++ )
//                {
//                    config_table[config][0]=2;
//                    config_table[config][i]=1;
//                    config_table[config][j]=1;
//                    config_table[config][k]=1;
//                    config_table[config][l]=1;
//                    config++;
//                }
//            }
//        }
//    }
////
////
////
//// d_6
//config_table[config][0]=1;
//config_table[config][1]=1;
//config_table[config][2]=1;
//config_table[config][3]=1;
//config_table[config][4]=1;
//config_table[config][5]=1;
//config++;


    return (config);
}







/*!
 * \fn factorial
 * \brief return factorial
 * \author C. Marchand
 */
int factorial( int p)
{
    int facts = 1;
    int i;
    for( i = 1; i<= p; i++)
        facts = facts * i;
    return( facts);
}

/*!
 * \fn combin
 * \brief return number of combination C(n,k)
 * \author C. Marchand
 */
int combin( int n, int r)
{
    return( factorial( n) / (factorial( r) * factorial(n- r) ) ) ;
}


double bayes(double M1 , double M2  )// M1 and M2 must be absolute values
{
   // return ( f(f(M1)+f(M2) ) );



    float dif;
    float min;

    // min-sum
    if (M1 < M2 )
    {
        min = M1;
        dif = M2 - M1;
    }
    else
    {
        min = M2;
        dif = M1 - M2;
    }

//////////modified bayes for syndrome
//    if (dif<0.05)
//    {
//        min=0.5*min;
//        //min=min-0.5;
//    }
//    else if (dif<0.2)
//
//    {
//        min=0.75*min;
//        //min=min-0.2;
//    }
//    else if (dif<1)
//    {
//        min=0.825*min;
//        //min=min-0.05;
//    }
//    else if (dif<2)
//    {
//        min=0.9375*min;
//        //min=min-0.01;
//    }
////////if (min<0) {min=0.0;}



//////modified bayes for bubble
    if (dif<0.1)
    {
        min=0.5*min;
    }
    else if (dif<0.2)

    {
        min=0.75*min;
    }
    else if (dif<1)
    {
        min=0.825*min;
    }
    else if (dif<2)
    {
        min=0.9375*min;
    }

    return(min);


}

/*****************************************************************/
double f(double x)        //f(x)=-ln (tanh(x/2))
{
    double y = x;
    if(y>16.6335)
    {
        y = 16.6335;
    }
    else
    {
        if(y<1.19452e-7)
        {
            y = 1.19452e-7;
        }
    }
    double res = log(1+exp(y))- log(-1+exp(y));
    return (res);
}

int check_deviation(int ** config_table,int dc_max,int index1,int index2,int index3,int * result)
{

    int dc;
    int sum=0;
    int temp;

//                printf("\n check deviation \n");
//                 printf("\n index1: %d \n ",index1);
//                for (dc=0;dc<dc_max;dc++)
//                {
//                    printf(" %d ", config_table[index1][dc]);
//                }
//                printf(" \n index2:%d \n",index2);
//                for (dc=0;dc<dc_max;dc++)
//                {
//                    printf(" %d ", config_table[index2][dc]);
//                }
//                printf(" \n index3:%d \n",index3);
//                for (dc=0;dc<dc_max;dc++)
//                {
//                    printf(" %d ", config_table[index3][dc]);
//                }



                for (dc=0;dc<dc_max;dc++)
                {

                    result[dc]=config_table[index1][dc]*config_table[index2][dc] ;
                }

                for (dc=0;dc<dc_max;dc++)
                {
                    temp = result[dc]*config_table[index3][dc] ;
                    result[dc]=temp;
                    sum=sum+temp;
                }


//                printf(" \n result \n");
//                for (dc=0;dc<dc_max;dc++)
//                {
//                    printf(" %d ", result[dc]);
//                }
//                printf("sum=%d",sum);

                return(sum);



}

int sort_config_table(int ** config_table,int config_table_size,int dc_max)
{

    //sorting config table

    //allocate cost to each config
    float config_cost[config_table_size];
    int config,dc;
    for (config=0; config<config_table_size; config++)
    {
        config_cost[config]=0;
        for (dc=0; dc<dc_max; dc++)
        {
            if (config_table[config][dc]>0)
            {
                config_cost[config]=config_cost[config] + config_table[config][dc]+ 3.0*dc;
            }

        }
    }

    //// check result
    //for (config=0; config<config_table_size; config++)
    //{
    //      printf(" \n \n");
    //    printf("%d :\t",config);
    //    for(dc=0; dc<dc_max; dc++)
    //        printf(" %d \t",config_table[config][dc]);
    //    printf(" \t = %f", config_cost[config]);
    //    getchar();
    //}




    //// sort cost
    syndrome_type syndrome_set[config_table_size];

    for (config=0; config<config_table_size; config++)
    {
     syndrome_set[config].LLR=config_cost[config];
     syndrome_set[config].config = config;
    }


     sorting(syndrome_set, config_table_size);

    //// check result
    //for (config=0; config<config_table_size; config++)
    //{
    //      printf(" \n \n");
    //    printf("%d :\t",config);
    //    for(dc=0; dc<dc_max; dc++)
    //        printf(" %d \t",config_table[syndrome_set[config].config][dc]);
    //    printf(" \t = %f", syndrome_set[config].LLR);
    //    getchar();
    //}


    //build new config_table
        int **config_table_temp;
    config_table_temp = AllocateArray(config_table_size, dc_max);
    for (config=0; config<config_table_size; config++)
    {
        for(dc=0; dc<dc_max; dc++)
            config_table_temp[config][dc] =  config_table[syndrome_set[config].config][dc];
    }

    for (config=0; config<config_table_size; config++)
    {
        for(dc=0; dc<dc_max; dc++)
            config_table[config][dc]=  config_table_temp[config][dc];
    }

    //// check result
    //for (config=0; config<config_table_size; config++)
    //{
    //      printf(" \n \n");
    //    printf("%d :\t",config);
    //    for(dc=0; dc<dc_max; dc++)
    //        printf(" %d \t",config_table[config][dc]);
    //    printf(" \t = %f", syndrome_set[config].LLR);
    //    getchar();
    //}

    return(0);
}


