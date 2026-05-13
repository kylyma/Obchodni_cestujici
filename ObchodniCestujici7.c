#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define pocet 6
//Last version

//gcc -o ObchodniCestujici7 ObchodniCestujici7.c -O3 -lm -Wall

void CreateMol(int N, char aaa[],int maxCycle,int cycPrint,float T, float k){

    FILE *molf;

     molf = fopen("trajectory.mol","a");
     int i;

    fprintf(molf,"! Soubor vygenerovan v ramci reseni problemu obchodniho cestujiciho.\n");
    fprintf(molf,"!Number of cities=%d, Steps=%d, Printed every %d steps, Initial Temperature= %f, k=%f\n",N+1,maxCycle,cycPrint,T,k);
    fprintf(molf,"\n");
    //fprintf(molf,"parameter_set = sea\n");
    fprintf(molf,"number_of_atoms = %d\n",N+1);
    fprintf(molf,"\n");
    fprintf(molf,"\n");
    fprintf(molf,"atoms\n");
    fprintf(molf,"! i  atom-id  a-type charge chir nbonds bound_atoms\n");


    fprintf(molf,"%3d %3d-H \t H      0    0      2     %d %d\n",0,0,1,N);

    for (i=1;i<N;i++)
    {
    fprintf(molf,"%3d %3d-H \t H      0    0      2     %d %d\n",i,i,i-1,i+1);
    }
    fprintf(molf,"%3d %3d-H \t H      0    0      2     %d %d\n",N,N,N-1,0);

    fclose(molf);
}

int isFileExists(const char *path){
    // Try to open file
    FILE *fptr = fopen(path, "r");

    // If file does not exists
    if (fptr == NULL)
        return 0;

    // File exists hence close file and return true.
    fclose(fptr);

    return 1;
}

void VypisKonfig(int N, char aaa[] ,float *X, float *Y){
    FILE *conf;
    conf = fopen(aaa, "wt+");
    for(int i=0;i<=N;i++){
    fprintf(conf,"%f \t %f\n",X[i],Y[i]);
    }
    fclose(conf);
}

void VypisSetup(int N, char aaa[] ,unsigned long int maxCycle, int cycPrint,float T0,float  k,int variant,float x,float y ){
    FILE *conf;
    conf = fopen(aaa, "wt+");
    fprintf(conf,"Cityes= %d\n",N);
    fprintf(conf,"Cyc= %lu \n",maxCycle);
    fprintf(conf,"Pr= %d\n",cycPrint);
    fprintf(conf,"T0= %f\n",T0);
    fprintf(conf,"k= %f\n",k);
    fprintf(conf,"var= %d\n",variant);
    fprintf(conf,"x= %f\n",x);
    fprintf(conf,"y= %f\n",y);
    fprintf(conf,"mix= 1\n");
    fclose(conf);
}

void VypisCas(char aaa[] ,double  seconds ){
    FILE *trip;
    int minuty, hodiny,dny;

    if(seconds<60.0f){
    trip = fopen(aaa, "a");
    fprintf(trip,"!Simulated annealing run time was %f sec\n", seconds);
    fclose(trip);
    printf("Simulated annealing run time was %f sec\n", seconds);
    }
    else if(seconds<3600.0f){
    minuty=seconds/60;
    seconds=seconds-minuty*60;

    trip = fopen(aaa, "a");
    fprintf(trip,"!Simulated annealing run time was %d minuts %f sec\n", minuty, seconds);
    fclose(trip);
    printf("Simulated annealing run time was %d minuts %f sec\n", minuty, seconds);

    }
    else if(seconds<86400.0f){
    hodiny=seconds/3600;
    seconds=seconds-hodiny*3600;
    minuty=seconds/60;
    seconds=seconds-minuty*60;

    trip = fopen(aaa, "a");
    fprintf(trip,"!Simulated annealing run time was %d hours %d minits %f sec\n",hodiny ,minuty, seconds);
    fclose(trip);
    printf("Simulated annealing run time was %d hours %d minits %f sec\n",hodiny ,minuty, seconds);
    }
    else{
    dny=seconds/86400;
    seconds=seconds-dny*86400;
    hodiny=seconds/3600;
    seconds=seconds-hodiny*3600;
    minuty=seconds/60;
    seconds=seconds-minuty*60;

    trip = fopen(aaa, "a");
    fprintf(trip,"!Simulated annealing run time was %d days %d hours %d minits %f sec\n",dny,hodiny ,minuty, seconds);
    fclose(trip);
    printf("Simulated annealing run time was %d days %d hours %d minits %f sec\n",dny,hodiny ,minuty, seconds);

    }


}

void VypisCasOpravy(char aaa[] ,double  seconds ){
    FILE *trip;
    int minuty, hodiny,dny;

    if(seconds<60.0f){
    trip = fopen(aaa, "a");
    fprintf(trip,"!My correction run time was %f sec\n", seconds);
    fprintf(trip,"!Trip length after my correction folows:\n");
    fclose(trip);
    printf("My correction run time was %f sec\n", seconds);
    //printf("Trip length after my correction folows:\n");
    }
    else if(seconds<3600.0f){
    minuty=seconds/60;
    seconds=seconds-minuty*60;

    trip = fopen(aaa, "a");
    fprintf(trip,"!My correction run time was %d minuts %f sec\n", minuty, seconds);
    //fprintf(trip,"!Trip length after my correction folows:\n");
    fclose(trip);
    printf("My correction run time was %d minuts %f sec\n", minuty, seconds);
    }
    else if(seconds<86400.0f){
    hodiny=seconds/3600;
    seconds=seconds-hodiny*3600;
    minuty=seconds/60;
    seconds=seconds-minuty*60;

    trip = fopen(aaa, "a");
    fprintf(trip,"!My correction run time was %d hours %d minits %f sec\n",hodiny ,minuty, seconds);
    fprintf(trip,"!Trip length after my correction folows:\n");
    fclose(trip);
    printf("My correction run time was %d hours %d minits %f sec\n",hodiny ,minuty, seconds);
    }
    else{
    dny=seconds/86400;
    seconds=seconds-dny*86400;
    hodiny=seconds/3600;
    seconds=seconds-hodiny*3600;
    minuty=seconds/60;
    seconds=seconds-minuty*60;

    trip = fopen(aaa, "a");
    fprintf(trip,"!My correction run time was %d days %d hours %d minits %f sec\n",dny,hodiny ,minuty, seconds);
    fprintf(trip,"!Trip length after my correction folows:\n");
    fclose(trip);
    printf("My correction run time was %d days %d hours %d minits %f sec\n",dny,hodiny ,minuty, seconds);

    }

}

void NactiSetup(char aaa[],int *N,unsigned long int *maxCycle,int *cycPrint, float *T0,float *k,int *variant,float *x,float *y,int *mix){
    FILE *conf;

     if (! isFileExists(aaa))
        {
            printf("File does not exists\n");
            exit(0);
        }

    conf = fopen(aaa, "r");
    if(fscanf(conf, "Cityes= %d\n",N)!=1) {
        printf("Can't read Cityes!!\n");
        exit(0);
    }
    //fscanf(conf, "Cyc= %llu\n",maxCycle);

    if(fscanf(conf, "Cyc= %lu \n",maxCycle)!=1){
        printf("Can't read Cycles!!\n");
        exit(0);
    }

    if(fscanf(conf, "Pr= %d\n",cycPrint)!=1){
        printf("Can't read Pr!!\n");
        exit(0);
    }

    if(fscanf(conf, "T0= %f\n",T0)!=1){
        printf("Can't read T0!!\n");
        exit(0);
    }

    if(fscanf(conf, "k= %f\n",k)!=1){
        printf("Can't read k!!\n");
        exit(0);
        }

    if(fscanf(conf, "var= %d\n",variant)!=1){
        printf("Can't read var!!\n");
        exit(0);
    }
    if(fscanf(conf, "x= %f \n",x)!=1){
        printf("Can't read x!!\n");
        exit(0);
    }
    if(fscanf(conf, "y= %f \n",y)!=1){
        printf("Can't read y !!\n");
        exit(0);
    }
    if(fscanf(conf, "mix= %d \n",mix)!=1){
        printf("Can't read mix !!\n");
        exit(0);
    }
    fclose(conf);
}

void NactiKonfig(int N, char aaa[] ,float *X,float  *Y ){
    FILE *conf;
    int i=0;
    strcat(aaa,".conf");

    if (! isFileExists(aaa)){
        printf("File does not exists\n");
        exit(0);
    }

    conf = fopen(aaa, "r");
    while(fscanf(conf,"%f \t %f", &X[i], &Y[i]) == 2){
        //printf("%f \t %f\n", X[i], Y[i]);
        if(i==N){break;}
        ++i;
    }
    fclose(conf);
}

void CreatePlb(char aaa[],float *hlavicka){
    FILE *plb;

    plb = fopen(strcat(aaa,".plb"), "wb");

    fwrite(hlavicka,1,sizeof(hlavicka),plb);
    fclose(plb);
}

void WriteTrip(char aaa[],float size,unsigned long  int step,int cycPrint ,float T0){
    FILE *trip;

    trip = fopen(aaa, "a");
    //fprintf(trip,"%15llu     %f     %f \n",step,T0,fabs(size));
    fprintf(trip,"%15ld     %lu      %f     %f \n",step,step/cycPrint,T0,fabs(size));
    fclose(trip);
}

float NahR(){
    return (float)rand()/(float)(RAND_MAX);
}

float NahRinterval(float zacatek, float konec){
    return (float)rand()/(float)(RAND_MAX)*(konec-zacatek)+zacatek;
}

int NahN(int N){
    return rand() % (N + 1 - 0) + 0;
}

float Cesta(float *X,float *Y,int NC,int N){
    if(NC==N){
        return sqrt((X[0]-X[NC])*(X[0]-X[NC])+(Y[0]-Y[NC])*(Y[0]-Y[NC]));
    }
    else{
        return sqrt((X[NC+1]-X[NC])*(X[NC+1]-X[NC])+(Y[NC+1]-Y[NC])*(Y[NC+1]-Y[NC]));
    }

}

float SumS(float *X,float *Y ,int N){
    float SS=0.0;
    for(int i=0;i<=N;i++){
    SS=SS+Cesta(X,Y,i,N);
    }
    return SS;
}

void PrintXY(float *X,float *Y,int N){
    printf("N  X  Y\n");
    for(int i=0;i<=N;i++){
        printf("%d  %f  %f\n",i,X[i],Y[i]);
    }

}

void mixVectors(int N, float *X, float *Y){
    float *X0,*Y0;
    int size = N;
    int *elements = malloc(sizeof(int)*size);

    // inizialize
    for (int i = 0; i < size; ++i)
        elements[i] = i;

    for (int i = size - 1; i > 0; --i) {
  // generate random index
        int w = rand()%i;
  // swap items
        int t = elements[i];
        elements[i] = elements[w];
        elements[w] = t;
    }

    X0 = malloc( N*sizeof(float*));
    Y0 = malloc( N*sizeof(float*));

    for(int i=0;i<N;i++){
        X0[i]=X[i];
    Y0[i]=Y[i];
    }

    for(int i=0;i<N;i++){
        X[i]=X0[elements[i]];
        Y[i]=Y0[elements[i]];
    }

    free(X0);
    free(Y0);
}

float MeziMesty(float *X,float *Y ,int Prvni,int Druhe){
    return sqrt((X[Prvni]-X[Druhe])*(X[Prvni]-X[Druhe])+(Y[Prvni]-Y[Druhe])*(Y[Prvni]-Y[Druhe]));
}

float RozdilDelek(float *X,float *Y ,int CB1,int CB2){
    return (MeziMesty(X,Y ,CB1-1,CB2)+MeziMesty(X,Y ,CB1,CB2+1))-(MeziMesty(X,Y ,CB1-1,CB1)+MeziMesty(X,Y ,CB2,CB2+1));
}

/*void NtaPermutace(int x//pocet clenu rady,int porady//,int fact, int *a ){
    //vraci a upravenene na permutaci cislo poradi-cislovan od nuly
    int i, j;
    int y=0;
    int c;
    while (y != fact) {
        if(y==poradi){
          break;
          }

        i=x-2;
        while(a[i] > a[i+1] ) i--;
        j=x-1;
        while(a[j] < a[i] ) j--;
      c=a[j];
      a[j]=a[i];
      a[i]=c;
    i++;
    for (j = x-1; j > i; i++, j--) {
  c = a[i];
  a[i] = a[j];
  a[j] = c;
      }
    y++;
   }

}*/

int factorial(int n){
    if (n == 1)
        return 1;
    else
        return n * factorial(n - 1);
}

void CreateA(int N,int mesto, int *a){
    //N je 49 pro 50 mest

    for(int i=0;i<pocet;i++){
        a[i]=i+mesto;
        if(a[i]>N){
            a[i]=a[i]-N-1;
        }
        if(a[i]<0){
            a[i]=a[i]+N+1;
        }
    }
}

int minimum(int *a,int n){
    int min=a[0];
    for(int i=0;i<n;i++){
        if(a[i]<min){
            min=a[i];
        }
    }
    return min;
}

int maximum(int *a,int n){
    int max=a[0];
    for(int i=0;i<n;i++){
        if(a[i]>max){
            max=a[i];
        }
    }
    return max;
}

float SpoctiCestu(int *aray, float *X,float *Y ,int mesto,int N){
    int seznam[pocet+2];
    float suma=0.0f;

    seznam[pocet+1]=pocet+mesto;
    if(seznam[pocet+1]>N){
        seznam[pocet+1]=seznam[pocet+1]-N-1;
    }
    if(seznam[pocet+1]<0){
        seznam[pocet+1]=seznam[pocet+1]+N+1;
    }
    seznam[0]=mesto-1;
    if(seznam[0]>N){
        seznam[0]=seznam[0]-N-1;
    }
    if(seznam[0]<0){
        seznam[0]=seznam[0]+N+1;
    }
    for(int i=1;i<=pocet;i++){
        seznam[i]=aray[i-1];
    }
    /*printf("Seznam: \n");
    for(int i=0;i<(pocet+2);i++){
        printf(" %d ",seznam[i]);
    }
    printf("\n");*/
    for(int i=0;i<(pocet+1);i++){
        suma=suma+MeziMesty(X,Y ,seznam[i],seznam[i+1]);
    }

    return suma;
}

void PrintVect(int *a,int n){
    for(int i=0;i<n;i++){
        printf("%d",a[i]);
    }
    printf("\n");
    printf("----------------------------\n");
}

void NtaPermutace(int N/*pocet mest*/,int poradi, int *D,int mesto ){ //N je 49 pro 50 mest
    //poradi-cislovan0 od nuly do fact-1
    unsigned int i, j;
    unsigned int y=0;
    unsigned int c,a[pocet],b[pocet];
    unsigned int x=pocet;
    int fact=factorial(pocet);

    CreateA( N, mesto, D);

    if (poradi>(fact-1)){
        printf("Funkce NtaPermutace:\n");
        printf("Zadano prilis velke poradi.\n Vystup jede znovu dokola!!\n");
        return;
    }

    for(int i=0;i<pocet;i++){
        b[i]=D[i];
    }

    for(int i=0;i<pocet;i++){
        a[i]=i;
    }

    while (y != fact) {
        if(y==poradi){
            for(int i=0;i<pocet;i++){
                D[a[i]]=b[i];
            }
    /*for(int i=0;i<pocet;i++){
         printf("%d",a[i]);
        }
    printf("----\n");*/
            break;
        }

        i=x-2;
        //while(a[i] > a[i+1] )
        while(a[i] > a[i+1] )
            i--;
            j=x-1;
        //printf("i%d,j%d\n",i,j);
        //while(a[j] < a[i] ) j--;
        while(a[j] < a[i] ) j--;
      c=a[j];
      a[j]=a[i];
      a[i]=c;
    i++;
    for (j = x-1; j > i; i++, j--) {
        c = a[i];
         a[i] = a[j];
        a[j] = c;
    }
    y++;
   }

}

int VratPoradi(int cislo, int N){

    if (cislo>=0){
        if(N>=cislo){
            return cislo;
        }
        else{
            return (cislo -N-1);
        }
    }
    else{
        return (N+cislo+1);
    }
}

int VratIndex(int *pomocny, int hledany){
    for(int i=0;i<pocet;i++){
       if(pomocny[i]==hledany){
           return i;
       }
    }

}

void OpravaMesto(float *X,float *Y ,int mesto,int N,int variant){//N je 49 pro 50 mest
    int poradi[pocet],nejKonf,pomocny[pocet];
    float delka,nejDelka,prestupX[pocet],prestupY[pocet];

    for(int i=0;i<factorial(pocet);i++){
        NtaPermutace(N,i,poradi,mesto);
        delka=variant*SpoctiCestu(poradi,X,Y,mesto,N);
        if(i==0){
            nejDelka=delka;
            nejKonf=i;
        }
        else if(delka<=nejDelka){
            nejDelka=delka;
            nejKonf=i;
        }
    }
    NtaPermutace(N,nejKonf,poradi,mesto);

    for(int i=0;i<pocet;i++){
        pomocny[i]=VratPoradi(i+mesto,N);
    }


    for(int i=0;i<pocet;i++){
        prestupX[i]=X[VratPoradi(i+mesto,N)];
        prestupY[i]=Y[VratPoradi(i+mesto,N)];
    }

    for(int i=0;i<pocet;i++){

       X[VratPoradi(i+mesto,N)]=prestupX[VratIndex(pomocny, poradi[i])];
       Y[VratPoradi(i+mesto,N)]=prestupY[VratIndex(pomocny, poradi[i])];
       /*X[i+mesto]=prestupX[poradi[i]-mesto];
       Y[i+mesto]=prestupY[poradi[i]-mesto];*/
    }
}

double MojeOprava(float *X,float *Y,int N,int variant){//N je 49 pro 50 mest

    clock_t start, finish;

    start = clock();
    for(int i=0;i<=N;i++){
    OpravaMesto(X,Y,i,N,variant);
    }
    finish = clock();
    return (((double)(finish-start))/CLOCKS_PER_SEC);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

    char NAME[255],NAME2[255],NAME3[255],NAME4[255],NAME5[255],NAME6[255],NAME7[255];
    int N,CB1,CB2,cycPrint=100,variant=1,mix=1;//maxCycle=5000,cycle=1;
    unsigned long  int maxCycle=5000,cycle=1;
    FILE *plb,*trip;
    float BoxSize[3],hlavicka[2],*X,*Y,z=0.0f,x,y,T,A,T0,k,finalSum,korr;//data-4castice 3 souradnice x y z
    double time_taken,totalTime=0.0f;
    srand(time(NULL));
    clock_t start, finish;
    time_t t;
    srand((unsigned) time(&t));

    BoxSize[0]=1.0f;// x box size
    BoxSize[1]=1.0f;//y box size
    BoxSize[2]=z;
    hlavicka[1]=-3.0;

    if((argc != 9) && (argc != 10) && (argc != 2)){
      printf("WRONG number of arguments.\n");
      printf("Program for solving Travelling salesman problem by simulated anealing.\n");
      printf("Number of cities. Must be higher than 4\n");
      printf("Program must be started as folows:\n");
      printf("  ObchodniCestujici7 N NAME Cyc Pr T0 k xsize ysize\n");
      printf("If other one argument is added, distance will be MAXIMALIZED, not minimalized.\n");
      printf("  N is number of cities.\n");
      printf("  NAME is name of the job = name of output files.\n");
      printf("  Cyc is number of cycles.\n");
      printf("  Pr is writing period in cycles- says how often data will be writed into output files\n");
      printf("  T0 is initial temperature.\n");
      printf("       Tempetadure is decreasing during simulation from Tinit by T=maxCycle*k*T0/(cycle+(maxCycle*k-1));\n");
      printf("  k is small constant-cca 0.0001- says how fast T decreeses- Bigger k=>slower\n");
      printf("  xsize and ysize are sizes of simulation box\n");
      printf("\n");
      printf("Example:\n");
      printf("ObchodniCestujici7 50 Vystup 5000000 20000 200 0.0001 1 1\n");
      printf("Output files NAME.mol NAME.plb and NAME.trip will be created.\n");
      printf("Other way to start is folowing:\n");
      printf("ObchodniCestujici7 NAME\n");
      printf("Files NAME.conf and NAME.set must exist.\n");
      printf("\n");
      printf("WARNING:\n");
      printf("    NAME.mol NAME.plb and NAME.trip will be owerwrited if exists!\n");
      printf("\n");
      printf("Compilation:\n");
      printf("gcc -o ObchodniCestujici7 ObchodniCestujici7.c -O3 -lm -Wall \n");
      exit(0);
    }

    if(argc==2){
        strcpy(NAME,argv[1]);
        strcpy(NAME2,argv[1]);
        strcpy(NAME3,argv[1]);
        strcpy(NAME4,argv[1]);
        strcpy(NAME5,argv[1]);
        strcpy(NAME6,argv[1]);
        strcpy(NAME7,argv[1]);

        strcat(NAME3,".trip");
        strcat(NAME5,".set");
        strcat(NAME4,".conf");
        strcat(NAME7,".fin");

        NactiSetup(NAME5,&N,&maxCycle,&cycPrint,&T0,&k,&variant,&x,&y,&mix);
        BoxSize[0]=x;// x box size
        BoxSize[1]=y;//y box size
        BoxSize[2]=z;
        hlavicka[1]=-3.0;
        hlavicka[0]=(float)(N);
        T=T0;
        X = malloc( N*sizeof(float*));
        Y = malloc( N*sizeof(float*));

        NactiKonfig(N-1, NAME6 , X , Y );

        if(mix!=0){
            mixVectors(N,X,Y);                      //Zde lze vypnout promichani vstupnich vektoru zakomentovanoim tohoto radku
        }

        N=N-1;
        CreatePlb(NAME2,hlavicka);
        CreateMol(N,NAME,maxCycle,cycPrint,T0,k);
        trip = fopen(NAME3, "w");
        //fprintf(trip,"!Number of cities=%d, Steps=%llu, Printed every %d steps, Initial Temperature= %f, k=%f\n",N+1,maxCycle,cycPrint,T0,k);
        fprintf(trip,"!Number of cities=%d, Steps=%lu , Printed every %d steps, Initial Temperature= %f, k=%f\n",N+1,maxCycle,cycPrint,T0,k);
        fprintf(trip,"!Steps   Frames   Temperature   Length of trip\n");
        fclose(trip);
    }

    else{
        sscanf(argv[1],"%d", &N);//atof
        //sscanf(argv[3],"%llu", &maxCycle);//atof
        sscanf(argv[3],"%lu ", &maxCycle);//atof
        sscanf(argv[4],"%d", &cycPrint);//atof
        sscanf(argv[5],"%f", &T0);//atof
        sscanf(argv[6],"%f", &k);//atof
        sscanf(argv[7],"%f", &x);
        sscanf(argv[8],"%f", &y);

        if(N<=4){
            printf("ERROR!!!\n");
            printf("N=%d\n",N);
            printf("Number of cities must be higher than 4\n");
            exit(0);
        }

        if(argc==10){
            variant=-1;
        }

        strcpy(NAME,argv[2]);
        strcpy(NAME2,argv[2]);
        strcpy(NAME3,argv[2]);
        strcpy(NAME4,argv[2]);
        strcpy(NAME5,argv[2]);
        strcpy(NAME7,argv[2]);

        strcat(NAME3,".trip");
        strcat(NAME4,".conf");
        strcat(NAME5,".set");
        strcat(NAME7,".fin");

        BoxSize[0]=x;// x box size
        BoxSize[1]=y;//y box size
        BoxSize[2]=z;
        hlavicka[1]=-3.0;
        hlavicka[0]=(float)(N);
        T=T0;

        X = malloc( N*sizeof(float*));
        Y = malloc( N*sizeof(float*));
        N=N-1;

        CreatePlb(NAME2,hlavicka);
        CreateMol(N,NAME,maxCycle,cycPrint,T0,k);
        trip = fopen(NAME3, "w");
        //fprintf(trip,"!Number of cities=%d, Steps=%llu, Printed every %d steps, Initial Temperature= %f, k=%f\n",N+1,maxCycle,cycPrint,T0,k);
        fprintf(trip,"!Number of cities=%d, Steps=%lu , Printed every %d steps, Initial Temperature= %f, k=%f\n",N+1,maxCycle,cycPrint,T0,k);
        fprintf(trip,"!Steps \t frames \t Temperature \t Length of trip\n");
        fclose(trip);

        //Generate initial configuration
        for(int i=0;i<=N;i++){
            X[i]=NahRinterval( 0.0f,x);
            Y[i]=NahRinterval( 0.0f,y);
        }
    }

    VypisKonfig(N, NAME4 ,X, Y);
    VypisSetup(N+1, NAME5, maxCycle, cycPrint, T0, k,variant, x, y);
    //Write initial configuration
    plb = fopen(NAME2, "ab");
    fwrite(BoxSize,1,sizeof(BoxSize),plb);

    for(int i=0;i<=N;i++){
        fwrite(&X[i],1,sizeof(float),plb);
        fwrite(&Y[i],1,sizeof(float),plb);
        fwrite(&z,1,sizeof(float),plb);
    }

    fclose(plb);
    WriteTrip(NAME3,variant*SumS(X,Y,N),cycle,cycPrint,T);

    //
    //mam poc konfiguraci i zapsanou
    //

    korr=sqrt(x*y);
    //korr=sqrt(x*x+y*y);
    start = clock();


    while(cycle<=maxCycle){

        T=maxCycle*k*T0/(cycle+(maxCycle*k-1));

        //zisk CB1 a CB2
        CB1=NahN(N);
        CB2=CB1;
        while(abs(CB1-CB2)<2 || abs(CB1-CB2)==N){
            CB2=NahN(N);
        }

        if(CB2<CB1){
            A=CB1;
            CB1=CB2;
            CB2=A;
        }

        if ((RozdilDelek(X,Y ,CB1,CB2)*variant)<0.0 || exp(-(RozdilDelek(X,Y ,CB1,CB2)*variant)/(T*korr))>NahR()){
        //zmena cesty dle CB1 a CB2 v X1Y1
            for(int i=0;i<=(CB2-CB1)/2;i++){
                A=X[CB1+i];
                X[CB1+i]=X[CB2-i];
                X[CB2-i]=A;
                A=Y[CB1+i];
                Y[CB1+i]=Y[CB2-i];
                Y[CB2-i]=A;
            }

        }
        //vypis cesty
        if(cycle % cycPrint==0){
            plb = fopen(NAME2, "ab");
            fwrite(BoxSize,1,sizeof(BoxSize),plb);
            for(int i=0;i<=N;i++){
                fwrite(&X[i],1,sizeof(float),plb);
                fwrite(&Y[i],1,sizeof(float),plb);
                fwrite(&z,1,sizeof(float),plb);
            }
            fclose(plb);
            WriteTrip(NAME3,variant*SumS(X,Y,N),cycle,cycPrint,T);
        }
        cycle++;
    }

    if(maxCycle%cycPrint!=0){
        plb = fopen(NAME2, "ab");
        fwrite(BoxSize,1,sizeof(BoxSize),plb);
        for(int i=0;i<=N;i++){
            fwrite(&X[i],1,sizeof(float),plb);
            fwrite(&Y[i],1,sizeof(float),plb);
            fwrite(&z,1,sizeof(float),plb);
        }
        fclose(plb);
        WriteTrip(NAME3,variant*SumS(X,Y,N),maxCycle,cycPrint,T);
    }

    finish = clock();
    time_taken = ((double)(finish-start))/CLOCKS_PER_SEC;
    /*trip = fopen(NAME3, "a");
    fprintf(trip,"!Simulated annealing run time was %f sec\n", time_taken);
    fclose(trip);
    printf("Simulated annealing run time was %f sec\n", time_taken);*/
    VypisCas(NAME3, time_taken );
    printf("Final trip length from simulated anealing is: %f\n",SumS(X,Y,N));

    //Moje Oprava
    //if(variant==1){
    for(int i=0;i<=3;i++){
        totalTime=MojeOprava(X,Y,N,variant)+totalTime;
        finalSum=SumS(X,Y,N);
        plb = fopen(NAME2, "ab");
        fwrite(BoxSize,1,sizeof(BoxSize),plb);
        for(int i=0;i<=N;i++){
            fwrite(&X[i],1,sizeof(float),plb);
            fwrite(&Y[i],1,sizeof(float),plb);
            fwrite(&z,1,sizeof(float),plb);
        }
        fclose(plb);
        WriteTrip(NAME3,finalSum,cycle,cycPrint,T);
    }
    VypisCasOpravy(NAME3,totalTime);
    printf("Final length after My correction is: %f\n",finalSum);
    printf("");
    VypisKonfig(N, NAME7 ,X, Y);
    //}

    free(X);
    free(Y);
    return 0;
}
