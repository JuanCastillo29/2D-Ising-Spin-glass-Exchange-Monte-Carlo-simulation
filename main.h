#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

typedef struct{
    int value;
    int vecinos[4];
}Espin;

typedef struct{
    Espin *lattice;
    int Size;
    double Beta;
}EspinLattice;

void Input( int *input, double *TMax, double *TMin);

void InicializarSistema(EspinLattice *sistema, int L);
void FreeSystemMemory(EspinLattice *sistemas);

#endif // MAIN_H_INCLUDED
