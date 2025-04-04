#ifndef ANALISIS_H_INCLUDED
#define ANALISIS_H_INCLUDED

typedef struct{
    int BinSize;
    int NBins;
    double *data;
    double varianza;
}BinAnalysis;

void binning(TimeEvolStorage *TimeSeries, int LData, int TermalData);

#endif // ANALISIS_H_INCLUDED
