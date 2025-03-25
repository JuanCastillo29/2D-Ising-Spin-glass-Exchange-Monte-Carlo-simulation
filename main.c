#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"

//Generador de números de Parisi-Rapuano.
void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}

//Generador de los números aleatorios según distribución gaussiana por el método de Box-Muller
double DistrGauss(void){
    double d1, d2;
    //Se generan los numeros planos necesarios
    d1 = Random();
    d2 = Random();
    if(d1 == 0)
        d1 = d2;
    //Se usa la fórmula para crear dos números independientes según una gaussiana
    return -sqrt(-2*log(d1))*cos(2*Pi*d2);
}

void Input( int *input, double *TMax, double *TMin){
    FILE *file = fopen("Input.txt", "r");
    if(!file){
        printf("Error opening input file.\n");
        exit(EXIT_FAILURE);
    }
    /*Input order: Lattice size, Number of temperatures which are going to be simulated,
    number of Monte Carlo steps between replica exchanges, number of Monte carlo steps between measures
    Total number of Monte Carlo steps, maximum temperature, minimum temperature, seed of random number generator*/
    if (fscanf(file, "%d\t%d\t%d\t%d\t%d\t%lf\t%lf %d\n", (input), (input +1), (input + 2), (input +3), (input + 4),
                TMax, TMin, (input + 5)) != 8){
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
    int spinvalue;
    //Starts with a random configuration.
    for (int i = 0; i < N; i++) {
        if(Random()<0.5)
            spinvalue = -1;
        else
            spinvalue = +1;
        sistema->lattice[i].value = spinvalue;
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

void Actualizarvecinos(Espin *spin, Espin *vecino, double Beta){

}

void PasoMonteCarlo(EspinLattice *sistema){
    int n;
    Espin spin;

    for(int i=0; i< sistema->Size; i++){
        n = Random()*sistema->Size;
        spin = sistema->lattice[n];

        //With a probability given by the Metropolis algorithm, I change the spin.
        if(Random()< spin.prob){

            //The change in the spin also changes the field that its neighbours see.
            for(int j=0; j<4; j++)
                Actualizarvecinos(&spin, &sistema->lattice[spin.vecinos[j]], sistema->Beta);

            //I change the system
            sistema->Energy += spin.dE;
            spin.value = -spin.value;
            spin.dE  = -spin.dE;
        }
    }
}

void CambioTemperatura(EspinLattice *sistema1, EspinLattice *sistema2){

}

int main()
{
    int inputdata[6];
    double TMax, TMin;
    Input(inputdata, &TMax, &TMin);
    ini_ran(inputdata[5]);

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
