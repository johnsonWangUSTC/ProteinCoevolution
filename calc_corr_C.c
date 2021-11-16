#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NaN 255
#define UNDEF -1

static int McLachlan[20][20]=   {{8, 1, 3, 4, 1, 3, 3, 2, 3, 2, 3, 3, 4, 3, 2, 4, 3, 3, 1, 1},
                                {1, 9, 1, 0, 0, 1, 3, 1, 0, 0, 3, 1, 0, 0, 1, 2, 2, 1, 2, 1},
                                {3, 1, 8, 5, 1, 3, 4, 0, 3, 1, 2, 5, 3, 4, 1, 3, 3, 1, 0, 1},
                                {4, 0, 5, 8, 0, 3, 2, 1, 4, 1, 1, 4, 4, 5, 3, 4, 4, 2, 1, 2},
                                {1, 0, 1, 0, 9, 0, 4, 3, 0, 5, 5, 0, 1, 0, 1, 2, 1, 3, 6, 6},
                                {3, 1, 3, 3, 0, 8, 2, 1, 3, 1, 1, 3, 3, 2, 3, 3, 2, 2, 1, 0},
                                {3, 3, 4, 2, 4, 2, 8, 2, 4, 2, 3, 4, 3, 4, 5, 3, 4, 2, 3, 4},
                                {2, 1, 0, 1, 3, 1, 2, 8, 1, 5, 5, 1, 1, 0, 1, 2, 3, 5, 3, 3},
                                {3, 0, 3, 4, 0, 3, 4, 1, 8, 2, 1, 4, 3, 4, 5, 3, 3, 2, 1, 1},
                                {2, 0, 1, 1, 5, 1, 2, 5, 2, 8, 6, 1, 1, 3, 2, 2, 3, 5, 3, 3},
                                {3, 3, 2, 1, 5, 1, 3, 5, 1, 6, 8, 2, 1, 3, 1, 2, 3, 4, 1, 2},
                                {3, 1, 5, 4, 0, 3, 4, 1, 4, 1, 2, 8, 1, 4, 3, 5, 3, 1, 0, 2},
                                {4, 0, 3, 4, 1, 3, 3, 1, 3, 1, 1, 1, 8, 3, 3, 3, 3, 2, 0, 0},
                                {3, 0, 4, 5, 0, 2, 4, 0, 4, 3, 3, 4, 3, 8, 5, 4, 3, 2, 2, 1},
                                {2, 1, 1, 3, 1, 3, 5, 1, 5, 2, 1, 3, 3, 5, 8, 4, 3, 2, 3, 2},
                                {4, 2, 3, 4, 2, 3, 3, 2, 3, 2, 2, 5, 3, 4, 4, 8, 5, 2, 3, 3},
                                {3, 2, 3, 4, 1, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 5, 8, 3, 2, 1},
                                {3, 1, 1, 2, 3, 2, 2, 5, 2, 5, 4, 1, 2, 2, 2, 2, 3, 8, 2, 3},
                                {1, 2, 0, 1, 6, 1, 3, 3, 1, 3, 1, 0, 0, 2, 3, 3, 2, 2, 9, 6},
                                {1, 1, 1, 2, 6, 0, 4, 3, 1, 3, 2, 2, 0, 1, 2, 3, 1, 3, 6, 9}};

static int BLOSUM62[24][24] = {{4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4}, 
                            {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4}, 
                            {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4} ,
                            {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4} ,
                            {0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4} ,
                            {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
                            {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4}, 
                            {0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4} ,
                            {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
                            {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 ,-3, -3, -1, -4}, 
                            {-1, -2, -3, -4, -1 ,-2, -3, -4, -3,  2,  4, -2 , 2 , 0 ,-3, -2, -1 ,-2, -1 , 1 ,-4, -3 ,-1, -4},
                            {-1,  2,  0, -1, -3 , 1,  1, -2, -1, -3 ,-2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4}, 
                            {-1, -1, -2, -3, -1,  0, -2 ,-3 ,-2,  1 , 2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4}, 
                            {-2, -3, -3, -3, -2, -3, -3, -3, -1 , 0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1 ,-4}, 
                            {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1 ,-4 ,-3, -2, -2, -1, -2, -4}, 
                            {1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2 , 0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4}, 
                            {0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4}, 
                            {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4}, 
                            {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 ,-3, -2, -1, -4},
                            {0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4}, 
                            {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4}, 
                            {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4}, 
                            {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4}, 
                            {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}}; 

static int EQUIV[24][24] = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}};

float corr_C(float* a, float* b, int M, float diag){
    int i, V=0;
    int *idx;
    float ma=0, mb=0, sa=0, sb=0, sab=0;
    idx = (int*)malloc(M*sizeof(int));
    for(i=0; i<M; i++){
        if(a[i]!=NaN && b[i]!=NaN){
            idx[i]=1;
            ma += a[i];
            mb += b[i]; 
            V++;       
        }  
        else{
            idx[i]=0;
        }
    }
    if(V==0){
        printf("All undefined\n");
        return(0.0);
    }
    ma /= V;
    mb /= V;
    for(i=0; i<M; i++){
        if(idx[i]!=0){
            sa += (a[i]-ma) * (a[i]-ma);
            sb += (b[i]-mb) * (b[i]-mb);
            sab += (a[i]-ma) * (b[i]-mb);
        }
    }
    free(idx);
    if(sa == 0 || sb == 0){
        printf("Zero Deviation\n");
        return(0.0);
    }
    return(sab/sqrt(sa)/sqrt(sb));
}

float ** calc_corr_C(int * MSA, int N, int L, int type){
    int * mat;
    int p;
    if(type==0){
        p = 20;
        mat = &McLachlan[0][0];
        printf("McLachlan\n");
    }
    else if(type==1){
        p = 24;
        mat = &BLOSUM62[0][0];
        printf("BLOSUM62\n");
    }
    else if(type==2){
        p = 24;
        mat = &EQUIV[0][0];
        printf("Equivalence\n");
    }
    else{
        int (*mat)[20];
        printf("WRONG TYPE!\n");
        return(NULL);
    }
    float ** DAT;
    float diag, res;
    int i, j, k, ai, aj;
    DAT = (float**)malloc(L*sizeof(float*));
    for(i=0; i<L; i++){
        DAT[i] = (float*)malloc(N*N*sizeof(float));
    }
    for(k=0; k<L; k++){
        printf("pos:%d\n",k);
        for(i=0; i<N; i++){
            for(j=0; j<=i; j++){
                ai = MSA[i*L+k];
                aj = MSA[j*L+k];
                if(ai!=UNDEF && aj!=UNDEF){
                    DAT[k][i*N+j] = DAT[k][j*N+i] = (float)(*(mat+ai*p+aj));
                }
                else{
                    DAT[k][i*N+j] = DAT[k][j*N+i] = NaN;
                }
            }
        }
    }
    //return(DAT);
    printf("Converted To Events!\n");
    float ** corr;
    corr = (float**)malloc(L*sizeof(float*));
    for(i=0;i<L;i++){
        corr[i] = (float*)malloc(L*sizeof(float));
    }
    for(i=0; i<L; i++){
        for(j=0; j<=i; j++){
            if(i==j){diag = 1;}
            else{diag = 0;}
            printf("%d,%d:",i,j);
            res = corr_C(DAT[i],DAT[j],N*N,diag);
            corr[i][j] = corr[j][i] = res;
            printf("%.2f\n",res);
        }
    }
    return(corr);
}

float ** calc_corr_C_v2(int * MSA, int N, int L, int type){
    int * mat;
    int p;
    if(type==0){
        p = 20;
        mat = &McLachlan[0][0];
        printf("McLachlan\n");
    }
    else if(type==1){
        p = 24;
        mat = &BLOSUM62[0][0];
        printf("BLOSUM62\n");
    }
    else if(type==2){
        p = 24;
        mat = &EQUIV[0][0];
        printf("Equivalence\n");
    }
    else{
        int (*mat)[20];
        printf("WRONG TYPE!\n");
        return(NULL);
    }
    float ** DAT;
    float diag, res;
    int i, j, k, ai, aj;
    DAT = (float**)malloc(L*sizeof(float*));
    for(i=0; i<L; i++){
        DAT[i] = (float*)malloc(N*(N-1)/2*sizeof(float));
    }
    for(k=0; k<L; k++){
        printf("pos:%d\n",k);
        for(i=1; i<N; i++){
            for(j=0; j<i; j++){
                ai = MSA[i*L+k];
                aj = MSA[j*L+k];
                if(ai!=UNDEF && aj!=UNDEF){
                    //if(i==N-1 && j==i-1){printf("index:%d,%d\n",i,j);}
                    DAT[k][i*(i-1)/2+j] = (float)(*(mat+ai*p+aj));
                }
                else{
                    //if(i==N-1 && j==i-1){printf("index:%d,%d\n",i,j);}
                    DAT[k][i*(i-1)/2+j] = NaN;
                }
            }
        }
    }
    //return(DAT);
    printf("Converted To Events!\n");
    float ** corr;
    printf("*\n");
    corr = (float**)malloc(L*sizeof(float*));
    for(i=0;i<L;i++){
        corr[i] = (float*)malloc(L*sizeof(float));
    }
    for(i=0; i<L; i++){
        for(j=0; j<=i; j++){
            if(i==j){diag = 1;}
            else{diag = 0;}
            printf("%d,%d:",i,j);
            res = corr_C(DAT[i],DAT[j],N*(N-1)/2,diag);
            corr[i][j] = corr[j][i] = res;
            printf("%.2e\n",res);
        }
    }
    return(corr);
}


float cov_C(float* a, float* b, int M, float diag){
    int i, V=0;
    int *idx;
    float ma=0, mb=0, sa=0, sb=0, sab=0;
    idx = (int*)malloc(M*sizeof(int));
    for(i=0; i<M; i++){
        if(a[i]!=NaN && b[i]!=NaN){
            idx[i]=1;
            ma += a[i];
            mb += b[i]; 
            V++;       
        }  
        else{
            idx[i]=0;
        }
    }
    if(V<=1){
        printf("Less than 2 samples\n");
        return((float)diag);
    }
    ma /= V;
    mb /= V;
    for(i=0; i<M; i++){
        if(idx[i]!=0){
            sab += (a[i]-ma) * (b[i]-mb);
        }
    }
    free(idx);
    return(sab/(V-1));
}


float ** calc_cov_C(int * MSA, int N, int L, int type){
    int * mat;
    int p;
    if(type==0){
        p = 20;
        mat = &McLachlan[0][0];
        printf("McLachlan\n");
    }
    else if(type==1){
        p = 24;
        mat = &BLOSUM62[0][0];
        printf("BLOSUM62\n");
    }
    else if(type==2){
        p = 24;
        mat = &EQUIV[0][0];
        printf("Equivalence\n");
    }
    else{
        int (*mat)[20];
        printf("WRONG TYPE!\n");
        return(NULL);
    }
    float ** DAT;
    float diag, res;
    int i, j, k, ai, aj;
    DAT = (float**)malloc(L*sizeof(float*));
    for(i=0; i<L; i++){
        DAT[i] = (float*)malloc(N*N*sizeof(float));
    }
    for(k=0; k<L; k++){
        printf("pos:%d\n",k);
        for(i=0; i<N; i++){
            for(j=0; j<=i; j++){
                ai = MSA[i*L+k];
                aj = MSA[j*L+k];
                if(ai!=UNDEF && aj!=UNDEF){
                    DAT[k][i*N+j] = DAT[k][j*N+i] = (float)(*(mat+ai*p+aj));
                }
                else{
                    DAT[k][i*N+j] = DAT[k][j*N+i] = NaN;
                }
            }
        }
    }
    //return(DAT);
    printf("Converted To Events!\n");
    float ** cov;
    cov = (float**)malloc(L*sizeof(float*));
    for(i=0;i<L;i++){
        cov[i] = (float*)malloc(L*sizeof(float));
    }
    for(i=0; i<L; i++){
        for(j=0; j<=i; j++){
            if(i==j){diag = 1;}
            else{diag = 0;}
            printf("%d,%d:",i,j);
            res = cov_C(DAT[i],DAT[j],N*N,diag);
            cov[i][j] = cov[j][i] = res;
            printf("%.2f\n",res);
        }
    }
    return(cov);
}