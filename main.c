#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "plots.h"
#include "analisis.h"

//Generador de n�meros de Parisi-Rapuano.
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


//Generador de los n�meros aleatorios seg�n distribuci�n gaussiana por el m�todo de Box-Muller
double DistrGauss(void){
    double d1, d2;
    //Se generan los numeros planos necesarios
    d1 = Random();
    d2 = Random();
    if(d1 == 0)
        d1 = d2;
    //Se usa la f�rmula para crear dos n�meros independientes seg�n una gaussiana
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


void InicializarSistema2D(EspinLattice *sistema, int L, double T){
    sistema->Beta = 1.0/T;
    int N = L*L;
    sistema->Size = N;
    sistema->lattice = (Espin *)malloc(N * sizeof(Espin));
    if (sistema->lattice == NULL) {
        printf("Error at memory allocation.\n");
        exit(1);
    }
    //Starts with a random configuration.
    for (int i = 0; i < N; i++) {
        if(Random()<0.5)
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
        sistema->Energy += ESpin/2; //I have to divide by two to avoid double counting.
        spin->dE = -2*ESpin;
        spin->prob=exp(-(spin->dE)*(sistema->Beta));

    }
}

void FreeSystemMemory(EspinLattice *sistema, TimeEvolStorage *fichero) {
    free(sistema->lattice);
    sistema->lattice = NULL;
    fclose(fichero->f);
    fichero->f = NULL;
}



void Actualizarvecinos(Espin *spin, Espin *vecino, double Beta, double *J, int bond_index){
    // Get the correct bond strength J_{ij}
    double Jij = J[bond_index];

    // Compute the change in the local field due to the spin flip
    double dH = -2 * Jij * spin->value;

    // Update neighbor's dE (change in energy if it flips)
    vecino->dE -= 2 * vecino->value * dH;

    // Update the transition probability using the Metropolis factor
    vecino->prob = exp(-vecino->dE * Beta);
}

void PasoMonteCarlo(EspinLattice *sistema, double *J){
    int n;
    Espin *spin;
    int bond_index;
    for(int i=0; i< sistema->Size; i++){
        n = Random()*sistema->Size;
        spin = &sistema->lattice[n];

        //With a probability given by the Metropolis algorithm, I change the spin.
        if(Random()< spin->prob){
            //The change in the spin also changes the field that its neighbours see.
            for(int j=0; j<4; j++){
                // Determine the correct bond index
                if (j == 0) bond_index = n;  // Right
                else if (j == 1) bond_index = sistema->Size + n;  // Up
                else if (j == 2) bond_index = spin->vecinos[2];  // Left
                else if (j == 3) bond_index = sistema->Size + spin->vecinos[3];  // Down
                Actualizarvecinos(spin, &sistema->lattice[spin->vecinos[j]], sistema->Beta, J,  bond_index);
            }
            //I change the system
            sistema->Energy += spin->dE;
            spin->value = -spin->value;
            spin->dE  = -spin->dE;
            spin->prob = 1.0/spin->prob;
        }
    }
}

void CambioTemperatura(TimeEvolStorage *fichero1, EspinLattice *sistema, TimeEvolStorage *fichero2 ){
    EspinLattice *sistema1 = &sistema[fichero1->Systemindex], *sistema2 = &sistema[fichero2->Systemindex];
    double db = sistema2->Beta - sistema1->Beta, dE = sistema2->Energy - sistema1->Energy;
    double Prob = exp(-db*dE);
    if(Random()<Prob){
        double aux = sistema2->Beta;
        sistema2->Beta = sistema1->Beta;
        sistema1->Beta = aux;
        int aux2 = fichero2->Systemindex;
        fichero2->Systemindex = fichero1->Systemindex;
        fichero1->Systemindex = aux2;
        for(int i = 0; i<sistema1->Size; i++){
            sistema1->lattice[i].prob = exp(-(sistema1->Beta)*sistema1->lattice[i].dE);
            sistema2->lattice[i].prob = exp(-(sistema2->Beta)*sistema2->lattice[i].dE);
        }
    }
}

void InicializarFicheros(double T, TimeEvolStorage *E, int i){
    char namefile[50];
    sprintf(namefile, "Energias/Temperatura=%lf.dat", T);
    E->f = fopen(namefile, "w+");
    E->T = T;
    E->Systemindex = i;
    if(E->f == NULL)
        printf("Error at openning the file: %s.", namefile);
}

void Simulacion(EspinLattice *sistemas, int Ttime, int NReplicas, int TMedida, int TExchange,double *J, TimeEvolStorage *ficheros){
    for(int t=0; t<Ttime; t++){
        for(int j=0; j<NReplicas; j++)
            PasoMonteCarlo(&sistemas[j], J);
        if(t%TMedida == 0){
            for(int j=0; j<NReplicas; j++)
                fprintf(ficheros[j].f, "%d \t %lf\n", t, sistemas[ficheros[j].Systemindex].Energy/sistemas[ficheros[j].Systemindex].Size);
        }
        if(t%TExchange == 0){
            for(int i= (t/TExchange)%2; i<NReplicas -1; i+=2)
                CambioTemperatura(&ficheros[i], sistemas, &ficheros[i+1]);
        }
    }
}

int main()
{
    int inputdata[6];
    double TMax, TMin;
    Input(inputdata, &TMax, &TMin);
    ini_ran(inputdata[5]);
    double dT = (TMax - TMin)/(inputdata[1]-1);
    /*int flag;
    printf("Please select the type of simulation you wish to perform:\n");
    printf("Enter '1' for a standard Ising model simulation or '2' for a spin glass simulation.\n");
    scanf("%d", &flag);
    if(flag != 1 && flag != 2){
        printf("Not valid input");
        return 0;
    }*/
    EspinLattice *sistemas = (EspinLattice *)malloc(inputdata[1]* sizeof(EspinLattice));
    if (sistemas == NULL) {
        printf("Error at memory allocation.\n");
        exit(1);
    }

    TimeEvolStorage *ficheros = (TimeEvolStorage *)malloc(inputdata[1]*sizeof(TimeEvolStorage));
    if(ficheros == NULL){
        printf("Error at memory allocation.\n");
        exit(1);
    }

    double *J = (double *)malloc(2*inputdata[0]*inputdata[0]*sizeof(double));
    if (J == NULL) {
        printf("Error at memory allocation.\n");
        exit(1);
    }
    for(int i=0; i<2*inputdata[0]*inputdata[0]; i++)
        J[i] = -1;

    // Inicializar cada sistema
    for (int i = 0; i < inputdata[1]; i++) {
        InicializarSistema2D(&sistemas[i], inputdata[0], TMin + dT * i);
        InicializarFicheros(TMin + dT * i, &ficheros[i], i);
        InicializarEnergias(&sistemas[i], J);
    }
    //Simulaci�n aqu�:
    Simulacion(sistemas, inputdata[4], inputdata[1], inputdata[3], inputdata[2], J, ficheros);

    for(int i=0; i< inputdata[1]; i++)
        binning(&ficheros[i],inputdata[4]*1.0/inputdata[3], 100);


    // Liberar memoria
    for (int i = 0; i <inputdata[1]; i++) {
        FreeSystemMemory(&sistemas[i], &ficheros[i]);
    }
    free(J);
    J = NULL;
    free(ficheros);
    ficheros = NULL;
    free(sistemas);
    sistemas = NULL;
    return 0;
}
