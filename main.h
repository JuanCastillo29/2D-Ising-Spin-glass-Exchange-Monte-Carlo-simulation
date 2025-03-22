#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

typedef struct{
    int value;
    int *vecinos;
}Espin;

typedef struct{
    Espin *lattice;
    int Size;
    double Beta;
}EspinLattice;

void Input( int *input, double *TMax, double *TMin);

#endif // MAIN_H_INCLUDED
