#include <math.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdio.h>
#include <string.h>

using namespace std;

#define BROJ_UZORAKA 32
#define MAX 64
#define SEMPLING_RATE 8000
#define NUM_FILTERS 20
#define BR_COEFF 1
/*
 * Compute the center frequency (fc) of the specified filter band (l) (Eq. 4)
 * This where the mel-frequency scaling occurs. Filters are specified so that their
 * center frequencies are equally spaced on the mel scale
 * Used for internal computation only - not the be called directly
 */
float GetCenterFrequency(uint16_t filterBand)
{
    float centerFrequency = 0.0f;
    float exponent;

    if(filterBand == 0)
    {
        centerFrequency = 0;
    }
    else if(filterBand >= 1 && filterBand <= 14)
    {
        centerFrequency = (200.0f * filterBand) / 3.0f;
    }
    else
    {
        exponent = filterBand - 14.0f;
        centerFrequency = pow(1.0711703, exponent);
        centerFrequency *= 1073.4;
    }

    return centerFrequency;
}

/*
 * Compute the band-dependent magnitude factor for the given filter band (Eq. 3)
 * Used for internal computation only - not the be called directly
 */
float GetMagnitudeFactor(uint16_t filterBand)
{
    float magnitudeFactor = 0.0f;

    if(filterBand >= 1 && filterBand <= 14)
    {
        magnitudeFactor = 0.015;
    }
    else if(filterBand >= 15 && filterBand <= 48)
    {
        magnitudeFactor = 2.0f / (GetCenterFrequency(filterBand + 1) - GetCenterFrequency(filterBand -1));
    }

    return magnitudeFactor;
}
/*
 * Compute the filter parameter for the specified frequency and filter bands (Eq. 2)
 * Used for internal computation only - not the be called directly
 */
float GetFilterParameter(uint16_t frequencyBand, uint16_t filterBand)
{
    float filterParameter = 0.0f;

    float boundary = (frequencyBand * SEMPLING_RATE) / BROJ_UZORAKA;		// k * Fs / N
    float prevCenterFrequency = GetCenterFrequency(filterBand - 1);		// fc(l - 1) etc.
    float thisCenterFrequency = GetCenterFrequency(filterBand);
    float nextCenterFrequency = GetCenterFrequency(filterBand + 1);

    if(boundary >= 0 && boundary < prevCenterFrequency)
    {
        filterParameter = 0.0f;
    }
    else if(boundary >= prevCenterFrequency && boundary < thisCenterFrequency)
    {
        filterParameter = (boundary - prevCenterFrequency) / (thisCenterFrequency - prevCenterFrequency);
        filterParameter *= GetMagnitudeFactor(filterBand);
    }
    else if(boundary >= thisCenterFrequency && boundary < nextCenterFrequency)
    {
        filterParameter = (boundary - nextCenterFrequency) / (thisCenterFrequency - nextCenterFrequency);
        filterParameter *= GetMagnitudeFactor(filterBand);
    }
    else if(boundary >= nextCenterFrequency && boundary < SEMPLING_RATE)
    {
        filterParameter = 0.0f;
    }

    return filterParameter;
}
/*
 * Computes the Normalization Factor (Equation 6)
 * Used for internal computation only - not to be called directly
 */
float NormalizationFactor(int m)
{
    float normalizationFactor = 0.0f;

    if(m == 0)
    {
        normalizationFactor = sqrt(1.0f / NUM_FILTERS);
    }
    else
    {
        normalizationFactor = sqrt(2.0f / NUM_FILTERS);
    }

    return normalizationFactor;
}
/*
 * Computes the specified (mth) MFCC
 *
 * spectralData - array of floats containing the results of FFT computation. This data is already assumed to be purely real
 * samplingRate - the rate that the original time-series data was sampled at (i.e 44100)
 * NumFilters - the number of filters to use in the computation. Recommended value = 48
 * binSize - the size of the spectralData array, usually a power of 2
 * m - The mth MFCC coefficient to compute
 *
 */
float GetCoefficient(float* spectralData, uint16_t m)
{
    float result = 0.0f;
    float outerSum = 0.0f;
    float innerSum = 0.0f;
    uint16_t k, l;

    result = NormalizationFactor(m);

    for(l = 1; l <= NUM_FILTERS; l++)
    {
        // Compute inner sum
        innerSum = 0.0f;
        for(k = 0; k < BROJ_UZORAKA - 1; k++)
        {
            innerSum += fabs(spectralData[k] * GetFilterParameter(k, l));
        }

        if(innerSum > 0.0f)
        {
            innerSum = log(innerSum); // The log of 0 is undefined, so don't use it
        }

        innerSum = innerSum * cos(((m * 3.14) / NUM_FILTERS) * (l - 0.5f));

        outerSum += innerSum;
    }

    result *= outerSum;

    return result;
}

int main(void)
{
    cout << "{";
    for(int i = 1; i <= NUM_FILTERS; i++)
    {
        //cout << "**** " << i << endl;
        //cout << "{";
        for(int j = 0; j < BROJ_UZORAKA - 1; j++)
        {
            cout << GetFilterParameter(j, i) << ", ";
        }
        //cout << "}," << endl;
    }
    cout << "}";

    /*cout << "{";
    for(int i = 1; i <= BR_COEFF; i++)
    {
        //cout << "**** " << i << endl;
        //cout << "{";
        for(int j = 1; j <= NUM_FILTERS; j++)
        {
            cout << cos(((i * 3.14) / NUM_FILTERS) * (j - 0.5f)) << ", ";
        }
        //cout << "}," << endl;
    }
    cout << "}";*/

    return 0;
}





