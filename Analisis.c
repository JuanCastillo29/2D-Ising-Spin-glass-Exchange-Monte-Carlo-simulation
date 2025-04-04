#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "main.h"
#include "Analisis.h"
#include "plots.h"

void FreeSystemMemoryBinning(BinAnalysis *Bins, FILE *f){
    free(Bins);
    Bins = NULL;
    fclose(f);
    f = NULL;
}

void binning(TimeEvolStorage *TimeSeries, int LData, int TermalData) {
    char filename2[100];
    sprintf(filename2, "Binning_Analysis/T=%.3lf.dat", TimeSeries->T);
    FILE *f = fopen(filename2, "w");
    if (f == NULL) {
        fprintf(stderr, "Could not open file: %s\n", filename2);
        return;
    }

    int DataLength = LData - TermalData;
    int NDivisions = (int)log2((double)DataLength / 100);
    if (NDivisions <= 0) {
        fprintf(stderr, "Error: Not enough data for binning analysis.\n");
        fclose(f);
        return;
    }

    rewind(TimeSeries->f);
    fflush(TimeSeries->f);

    // Skip thermalization lines
    char buffer[1024];
    for (int i = 0; i < TermalData; i++) {
        if (!fgets(buffer, sizeof(buffer), TimeSeries->f)) {
            fprintf(stderr, "Error: Failed to skip thermalization data.\n");
            fclose(f);
            return;
        }
    }

    BinAnalysis *Bins = (BinAnalysis *)malloc(NDivisions * sizeof(BinAnalysis));
    if (!Bins) {
        perror("Memory allocation error.");
        fclose(f);
        return;
    }

    for (int i = 0; i < NDivisions; i++) {
        Bins[i].BinSize = (int)pow(2, i);
        Bins[i].NBins = DataLength / Bins[i].BinSize;
        Bins[i].varianza = 0;
        Bins[i].data = (double *)calloc(Bins[i].NBins, sizeof(double));
    }

    // Temporary array to store raw data
    double *raw_data = (double *)malloc(DataLength * sizeof(double));
    if (!raw_data) {
        perror("Memory allocation failed for raw_data.");
        fclose(f);
        free(Bins);
        return;
    }

    // Read data into raw_data[]
    for (int i = 0; i < DataLength; i++) {
        if (fgets(buffer, sizeof(buffer), TimeSeries->f)) {
            sscanf(buffer, "%*d\t%lf", &raw_data[i]);
        } else {
            fprintf(stderr, "Error reading data at line %d\n", TermalData + i);
            fclose(f);
            free(raw_data);
            free(Bins);
            return;
        }
    }

    // Fill bins
    for (int j = 0; j < NDivisions; j++) {
        int bin_size = Bins[j].BinSize;
        for (int i = 0; i < Bins[j].NBins; i++) {
            double sum = 0.0;
            for (int k = 0; k < bin_size; k++) {
                sum += raw_data[i * bin_size + k];
            }
            Bins[j].data[i] = sum / bin_size;
        }

        // Compute mean of the binned data
        double mean = 0.0;
        for (int i = 0; i < Bins[j].NBins; i++) {
            mean += Bins[j].data[i];
        }
        mean /= Bins[j].NBins;

        // Compute variance
        double variance = 0.0;
        for (int i = 0; i < Bins[j].NBins; i++) {
            variance += pow(Bins[j].data[i] - mean, 2);
        }
        variance /= (Bins[j].NBins - 1);

        fprintf(f, "%d\t%.10lf\n", bin_size, sqrt(variance/Bins[j].NBins));

        free(Bins[j].data);
        Bins[j].data = NULL;
    }

    printf("Binning analysis completed successfully.\n");

    fclose(f);
    free(Bins);
    free(raw_data);
    plotBinningAnalysis(TimeSeries);
}
