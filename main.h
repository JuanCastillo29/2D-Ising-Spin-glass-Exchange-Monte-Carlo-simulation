#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
#define NormRANu (2.3283063671E-10F)
#define Pi 3.14159265


typedef struct{
    int value;
    unsigned int vecinos[4];
    double dE;
    float prob;
}Espin;


typedef struct{
    Espin *lattice;
    unsigned int Size;
    double Beta;
    double Energy;
}EspinLattice;

void Input( int *input, double *TMax, double *TMin);

void InicializarSistema(EspinLattice *sistema, int L);
void FreeSystemMemory(EspinLattice *sistemas);

#endif // MAIN_H_INCLUDED
