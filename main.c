#include <stdio.h>
#include <stdlib.h>

#include "main.h"

void Input( int *input, double *TMax, double *TMin){
    FILE *file = fopen("Input.txt", "r");
    if(!file){
        printf("Error opening input file.\n");
        exit(EXIT_FAILURE);
    }
    /*Input order: Lattice size, Number of temperatures which are going to be simulated,
    number of Monte Carlo steps between replica exchanges, number of Monte carlo steps between measures
    Total number of Monte Carlo steps, maximum temperature, minimum temperature.*/
    if (fscanf(file, "%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n", (input), (input +1), (input + 2), (input +3), (input + 4),
                TMax, TMin) != 7){
        printf("Error reading input parameters.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }
    fclose(file);
}


void InicializarSistema(EspinLattice *sistema, int L){
    int N = L*L;
    sistema->Size = N;
    sistema->lattice = (Espin *)malloc(N * sizeof(Espin));
    if (sistema->lattice == NULL) {
        printf("Error at memory allocation.\n");
        exit(1);
    }
    //Starts with a random configuration.
    for (int i = 0; i < N; i++) {
        sistema->lattice[i].value = (rand() % 2) ? 1 : -1;
    }

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = i * L + j;  // Convert (i, j) to 1D index

            // Assign neighbors with periodic boundary conditions
            sistema->lattice[idx].vecinos[0] = i * L + (j + 1) % L;         // Right
            sistema->lattice[idx].vecinos[1] = i * L + (j - 1 + L) % L;     // Left
            sistema->lattice[idx].vecinos[2] = ((i + 1) % L) * L + j;       // Down
            sistema->lattice[idx].vecinos[3] = ((i - 1 + L) % L) * L + j;   // Up
        }
    }
}

void FreeSystemMemory(EspinLattice *sistema) {
    free(sistema->lattice);
}

int main()
{
    int inputdata[5];
    double TMax, TMin;
    Input(inputdata, &TMax, &TMin);

    EspinLattice *sistemas = (EspinLattice *)malloc(2*inputdata[1]* sizeof(EspinLattice));
    if (sistemas == NULL) {
        printf("Error at memory allocation.\n");
        exit(1);
    }

    // Inicializar cada sistema
    for (int i = 0; i < 2*inputdata[1]; i++) {
        InicializarSistema(&sistemas[i], inputdata[0]);
    }
    //Simulación aquí:

    // Liberar memoria
    for (int i = 0; i < 2*inputdata[1]; i++) {
        FreeSystemMemory(&sistemas[i]);
    }
    free(sistemas);

    return 0;
}
