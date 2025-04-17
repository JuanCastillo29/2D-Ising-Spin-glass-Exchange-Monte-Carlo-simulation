#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
#define NormRANu (2.3283063671E-10F)
#define Pi 3.14159265
#define MAXDIM  3
#define MAXREPLICAS 2


typedef struct{
    int value;
    unsigned int vecinos[2*MAXDIM];
    double dE;
    float prob;
}Espin;


typedef struct{
    Espin *lattice;
    int Dim;
    unsigned int Size;
    double Beta;
    double Energy;
}EspinLattice;


typedef struct{
    FILE *f;
    int Systemindex[MAXREPLICAS];
    double T;
} TimeEvolStorage;

void Input( int *input, double *TMax, double *TMin);

void InicializarSistema(EspinLattice *sistema, int L, double T);
void FreeSystemMemory(EspinLattice *sistemas, TimeEvolStorage *fichero);

#endif // MAIN_H_INCLUDED
