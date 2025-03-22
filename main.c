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



int main()
{
    int inputdata[5];
    double TMax, TMin;
    Input(inputdata, &TMax, &TMin);
    printf("%lf\n", TMin);
    printf("%d\n", inputdata[2] );
    return 0;
}
