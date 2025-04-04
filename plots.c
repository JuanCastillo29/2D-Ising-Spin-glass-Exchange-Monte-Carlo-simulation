#include <stdio.h>
#include <stdlib.h>

#include "main.h"

void PlotIsingModel(FILE *f) {
    // Flush and rewind before reading
    fflush(f);
    rewind(f); // Go to the beginning of the file
    if (f == NULL) {
        printf("Error: File pointer is NULL.\n");
        return;
    }

    // Create a temporary file to store the data
    char tempFilename[] = "temp_gnuplot_data.txt";
    FILE *tempFile = fopen(tempFilename, "w");
    if (tempFile == NULL) {
        printf("Error: Could not create temporary file.\n");
        return;
    }

    // Copy content from f to tempFile

    if (fgetc(f) == EOF) {
    printf("Error: Input file is empty.\n");
    fclose(tempFile);
    return;
    }
    rewind(f);
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), f) != NULL) {
        fputs(buffer, tempFile);
    }

    fclose(tempFile); // Close temp file after writing

    // Open GNUplot and plot the temporary file
    FILE *gnuplotPipe = _popen("gnuplot -persist", "w");
    if (gnuplotPipe == NULL) {
        printf("Error: Could not open GNUplot.\n");
        return;
    }
    fprintf(gnuplotPipe, "set xlabel 't [MCs]'\n");
    fprintf(gnuplotPipe, "set ylabel 'E/N'\n");
    fprintf(gnuplotPipe, "plot '%s' using 1:2 pt 7 title 'kT = 2.0'\n", tempFilename);
    fflush(gnuplotPipe);

    _pclose(gnuplotPipe); // Close GNUplot

    // Optional: Remove the temp file after plotting
    remove(tempFilename);
}

// === Write Gnuplot Script and Export to PDF ===

void plotBinningAnalysis(TimeEvolStorage *TimeSeries){
    FILE *gp = fopen("plot_binning.gp", "w");
    if (gp != NULL) {
        fprintf(gp,
            "set terminal pdfcairo size 4,3 enhanced font 'Arial,12'\n"
            "set output 'Binning_Analysis/binning_plot_T=%.3f.pdf'\n"
            "set title 'Binning Analysis for T = %.3f'\n"
            "set xlabel '{/:Bold Bin Size}'\n"
            "set ylabel '{/: Average Standard Deviation}'\n"
            "set grid\n"
            "set logscale x\n"
            "unset key\n"
            "plot 'Binning_Analysis/T=%.3f.dat' using 1:2 pt 7 lc rgb 'blue'\n",
            TimeSeries->T, TimeSeries->T, TimeSeries->T
        );
        fclose(gp);

        // Run gnuplot to generate PDF
        int ret = system("gnuplot plot_binning.gp");
        if (ret != 0) {
            fprintf(stderr, "Gnuplot execution failed.\n");
        }   else {
            printf("Plot exported to PDF: Binning_Analysis/binning_plot_T=%.3f.pdf\n", TimeSeries->T);
        }
    }   else {
        fprintf(stderr, "Could not create Gnuplot script.\n");
    }
}
