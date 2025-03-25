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
    //Starts with a random configuration.
    for (int i = 0; i < N; i++) {
        if(/*Random()<0.5*/ i%2==0)
            sistema->lattice[i].value  = -1;
        else
            sistema->lattice[i].value  = +1;
    }

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = i * L + j;  // Convert (i, j) to 1D index

            // Assign neighbors with periodic boundary conditions
            sistema->lattice[idx].vecinos[0] = i * L + (j + 1) % L;         // Right
            sistema->lattice[idx].vecinos[2] = i * L + (j - 1 + L) % L;     // Left
            sistema->lattice[idx].vecinos[3] = ((i + 1) % L) * L + j;       // Down
            sistema->lattice[idx].vecinos[1] = ((i - 1 + L) % L) * L + j;   // Up
        }
    }
}

double Localfield(EspinLattice *sistema, double *J,unsigned int *n, int index){
    double h = 0;
    h += J[index]* sistema->lattice[*(n)].value; //Right bond contribution.
    h += J[sistema->Size + index]* sistema->lattice[*(n+1)].value; //Up bond contribution.
    h += J[*(n+2)]* sistema->lattice[*(n+2)].value; //Left bond contribution.
    h += J[sistema->Size + *(n+3)]* sistema->lattice[*(n+3)].value; //Down bond contribution.
    return h;
}

void InicializarEnergias(EspinLattice *sistema, double *J){
    Espin *spin;
    sistema->Energy=0;
    double ESpin;
    for(int i=0; i< sistema->Size; i++){
        spin = &sistema->lattice[i];
        ESpin = spin->value*Localfield(sistema, J, spin->vecinos, i);
        sistema->Energy += ESpin/2;
        spin->dE = -2*ESpin;
        spin->prob=exp(-(spin->dE)*(sistema->Beta));

    }
}

void FreeSystemMemory(EspinLattice *sistema) {
    free(sistema->lattice);
}

void Actualizarvecinos(Espin *spin, Espin *vecino, double Beta){

}

void PasoMonteCarlo(EspinLattice *sistema){
    int n;
    Espin *spin;

    for(int i=0; i< sistema->Size; i++){
        n = Random()*sistema->Size;
        spin = &sistema->lattice[n];

        //With a probability given by the Metropolis algorithm, I change the spin.
        if(Random()< spin->prob){
            //The change in the spin also changes the field that its neighbours see.
            for(int j=0; j<4; j++)
                Actualizarvecinos(spin, &sistema->lattice[spin->vecinos[j]], sistema->Beta);

            //I change the system
            sistema->Energy += spin->dE;
            spin->value = -spin->value;
            spin->dE  = -spin->dE;
        }
    }
}

void CambioTemperatura(EspinLattice *sistema1, EspinLattice *sistema2){
    double db = sistema2->Beta - sistema1->Beta, dE = sistema2->Energy - sistema1->Energy;
    double Prob = exp(-db*dE);
    if(Random()<Prob){
        double aux = sistema2->Beta;
        sistema2->Beta = sistema1->Beta;
        sistema1->Beta = aux;
        for(int i = 0; i<sistema1->Size; i++){
            sistema1->lattice[i].prob = exp(-(sistema1->Beta)*sistema1->lattice[i].dE);
            sistema2->lattice[i].prob = exp(-(sistema2->Beta)*sistema2->lattice[i].dE);
        }
    }
}

void Simulacion(EspinLattice *sistemas, int Ttime, int NReplicas, int TMedida, int TExchange){
    for(int t=0; t<Ttime; t++){
        for(int j=0; j<NReplicas; j++)
            PasoMonteCarlo(&sistemas[j]);
       /* if(t%TMedida)
        */
        /*if(t%TExchange){
            for(int i=0; i<NReplicas -2; i++)
                CambioTemperatura(&sistemas[i], &sistemas[i+2]);
        */

    }
}

int main()
{
    int inputdata[6];
    double TMax, TMin;
    Input(inputdata, &TMax, &TMin);
    ini_ran(inputdata[5]);
    double J[200];
    for(int i=0; i<200; i++)
        J[i]=1;
    EspinLattice sistema;
    InicializarSistema(&sistema, 10);
    InicializarEnergias(&sistema, J);
    printf("%lf\n", sistema.Energy);
    FreeSystemMemory(&sistema);
    /*EspinLattice *sistemas = (EspinLattice *)malloc(2*inputdata[1]* sizeof(EspinLattice));
    if (sistemas == NULL) {
        printf("Error at memory allocation.\n");
        exit(1);
    }
    double J[inputdata[0]*inputdata[0]];
    for(int i=0; i<inputdata[0]*inputdata[0]; i++)
        J[i] = DistrGauss();

    // Inicializar cada sistema
    for (int i = 0; i < 2*inputdata[1]; i++) {
        InicializarSistema(&sistemas[i], inputdata[0]);
        InicializarEnergias(&sistemas[i], J);
    }
    //Simulación aquí:


    // Liberar memoria
    for (int i = 0; i < 2*inputdata[1]; i++) {
        FreeSystemMemory(&sistemas[i]);
    }
    free(sistemas);
*/
    return 0;
}
