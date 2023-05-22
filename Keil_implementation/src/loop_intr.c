#include <arm_math.h> 
#include "audio.h"
#include "LOOKUP.h"

#define M_PI 3.1415926535897932384
#define BROJ_UZORAKA 128
#define MAX 128
#define VRIJEME 5
#define BROJ_PJESAMA 6

float32_t x, y, state[N1];
arm_fir_instance_f32 S;


typedef struct {
  float real;
  float imag;
} COMPLEX;
int i, j;

int log2_moj(int N)    /*function to calculate the log2_moj(.) of int numbers*/
{
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}

int reverse(int N, int n)    //calculating revers number
{
  int j, p = 0;
  for(j = 1; j <= log2_moj(N); j++) {
    if(n & (1 << (log2_moj(N) - j)))
      p |= 1 << (j - 1);
  }
  return p;
}

void ordina(COMPLEX *f1, int N) //using the reverse order in the array
{
  COMPLEX f2[MAX];
  for(i = 0; i < N; i++)
    f2[i] = f1[reverse(N, i)];
  for(j = 0; j < N; j++)
    f1[j] = f2[j];
}

COMPLEX polar(float r, float fi){
	COMPLEX rez;
	rez.real = r * arm_cos_f32(fi);
	rez.imag = r * arm_sin_f32(fi);
	return rez;
}

COMPLEX pomnozi_complex(COMPLEX a, COMPLEX b) {
	COMPLEX rez;
	rez.real = a.real * b.real - a.imag * b.imag;
	rez.imag = a.real * b.imag + a.imag * b.real;
	return rez;
}

COMPLEX brzo_stepenuj(COMPLEX a, int16_t b){
	COMPLEX rez;
	rez.real = 1;
	rez.imag = 0;
	
	if(b == 0) return rez;
	if(b == 1) return a;
	
	rez = brzo_stepenuj(a, b >> 1);
	rez = pomnozi_complex(rez, rez);
	
	if(b % 2 == 1) rez = pomnozi_complex(rez, a);
	
	return rez;
}
COMPLEX W[BROJ_UZORAKA*2];
uint16_t n;
uint16_t a;
void transform( COMPLEX *f, int N) //
{
  ordina(f, N);    //first: reverse order
  W[1] = polar(1., -2. * M_PI / N);
  W[0].real = 1;
	W[0].imag = 0;
	n = 1;
	a = N >> 1;
  for(i = 2; i < N / 2; i++)
    W[i] = brzo_stepenuj(W[1], i);
  
  for(j = 0; j < log2_moj(N); j++) {
    for(i = 0; i < N; i++) {
      if(!(i & n)) {
        COMPLEX temp = f[i];
        COMPLEX Temp = pomnozi_complex(W[(i * a) % (n * a)], f[i + n]);
        f[i].real = temp.real + Temp.real;
				f[i].imag = temp.imag  + Temp.imag ;
        f[i + n].real = temp.real - Temp.real;
				f[i + n].imag = temp.imag - Temp.imag;
      }
    }
    n = n << 1;
    a = a >> 1;
  }
}
float izlaz[BROJ_UZORAKA];
void FFT(COMPLEX* f, int N, float d)
{
  transform(f, N);
  for(i = 0; i < N; i++) {
    f[i].real *= d; //multiplying by step
		f[i].imag *= d;
		izlaz[i] = sqrt(f[i].real * f[i].real + f[i].imag * f[i].imag);
		
		if(i <= 20 || i >= 107) {
			if(izlaz[i] < 5) izlaz[i] = 0;
			else if(izlaz[i] < 10) izlaz[i] = 1;
			else if(izlaz[i] < 15) izlaz[i] = 2;
			else if(izlaz[i] < 20) izlaz[i] = 3;
			else if(izlaz[i] < 25) izlaz[i] = 4;
			else if(izlaz[i] < 30) izlaz[i] = 5;
			else if(izlaz[i] < 35) izlaz[i] = 6;
			else izlaz[i] = 7;
		}
		if((i > 20 && i <= 32) || (i >= 95 && i < 107)) {
			if(izlaz[i] < 4) izlaz[i] = 0;
			else if(izlaz[i] < 8) izlaz[i] = 1;
			else if(izlaz[i] < 12) izlaz[i] = 2;
			else if(izlaz[i] < 16) izlaz[i] = 3;
			else if(izlaz[i] < 20) izlaz[i] = 4;
			else if(izlaz[i] < 24) izlaz[i] = 5;
			else if(izlaz[i] < 28) izlaz[i] = 6;
			else izlaz[i] = 7;
		}
		if((i > 32 && i <= 64) || (i >= 64 && i < 95)) {
			if(izlaz[i] < 3) izlaz[i] = 0;
			else if(izlaz[i] < 6) izlaz[i] = 1;
			else if(izlaz[i] < 9) izlaz[i] = 2;
			else if(izlaz[i] < 12) izlaz[i] = 3;
			else if(izlaz[i] < 15) izlaz[i] = 4;
			else if(izlaz[i] < 18) izlaz[i] = 5;
			else if(izlaz[i] < 21) izlaz[i] = 6;
			else izlaz[i] = 7;
		}
		
		
		
		//if(izlaz[i] < 10) izlaz[i] = 0;
		//else							izlaz[i] = 100;
	}
}

volatile int16_t audio_chR=0;    //16 bits audio data channel right
volatile int16_t audio_chL=0;    //16 bits audio data channel left

COMPLEX ulaz[BROJ_UZORAKA];

int br = 0; 
uint16_t br_prozora = 0; 

/*
 * Computes the Normalization Factor (Equation 6)
 * Used for internal computation only - not to be called directly
 */
float NormalizationFactor(int NumFilters, int m)
{
    float normalizationFactor = 0.0f;

    if(m == 0)
    {
        normalizationFactor = sqrt(1.0f / NumFilters);
    }
    else
    {
        normalizationFactor = sqrt(2.0f / NumFilters);
    }

    return normalizationFactor;
}

float GetCoefficient(float* spectralData, uint16_t samplingRate, uint16_t NumFilters, uint16_t binSize, uint16_t m)
{
    float result = 0.0f;
    float outerSum = 0.0f;
    float innerSum = 0.0f;
    uint16_t k, l;

    result = NormalizationFactor(NumFilters, m);

    for(l = 1; l <= NumFilters; l++)
    {
        // Compute inner sum
        innerSum = 0.0f;
         for(k = 0; k < binSize - 1; k++)
        {
					if(h1[(l-1)*(BROJ_UZORAKA - 1) + k]!=0)
            innerSum += fabs(spectralData[k] * h1[(l-1)*(BROJ_UZORAKA - 1) + k]);
						//GetFilterParameter(samplingRate, binSize, k, l));
        }

        if(innerSum > 0.0f)
        {
            innerSum = log(innerSum); // The log of 0 is undefined, so don't use it
        }

				
        //innerSum = innerSum * kosinus[(m - 1) * NumFilters + l - 1];
				innerSum = innerSum * arm_cos_f32(((m * 3.14) / NumFilters) * (l - 0.5f));

        outerSum += innerSum;
        //cout << "outerSum: " << outerSum << endl;
    }

    result *= outerSum;

    return result;
}
const int mfcc_velicina = VRIJEME * 8000 / BROJ_UZORAKA;
float mfcc1[2 + mfcc_velicina];
float mfcc2[2 + mfcc_velicina];
float mfcc3[2 + mfcc_velicina];
float mfcc4[2 + mfcc_velicina];
float mfcc5[2 + mfcc_velicina];

_Bool z = 0;
_Bool z2 = 0;

float cor1[2 * H_duzina - 1];
float max_cor1;
uint32_t pos_cor1;

float cor2[2 * H_duzina - 1];
float max_cor2;
uint32_t pos_cor2;

float cor3[2 * H_duzina - 1];
float max_cor3;
uint32_t pos_cor3;

float cor4[2 * H_duzina - 1];
float max_cor4;
uint32_t pos_cor4;

float cor5[2 * H_duzina - 1];
float max_cor5;
uint32_t pos_cor5;


float niz_max_cor1[BROJ_PJESAMA];
float niz_max_cor2[BROJ_PJESAMA];
float niz_max_cor3[BROJ_PJESAMA];
float niz_max_cor4[BROJ_PJESAMA];
float niz_max_cor5[BROJ_PJESAMA];

float niz_max[BROJ_PJESAMA];

void I2S_HANDLER(void) {   /****** I2S Interruption Handler *****/
	audio_IN = i2s_rx();	
	audio_chL = (audio_IN & 0x0000FFFF);       //Separate 16 bits channel left
	audio_chR = ((audio_IN >>16)& 0x0000FFFF); //Separate 16 bits channel right
	
	x = audio_chL;
	arm_fir_f32(&S, &x, &y, 1);
	y=x;
	ulaz[br].real = (y);
	ulaz[br].imag = 0;

	if(br >= BROJ_UZORAKA){
			z = 1;
			br = 0;
			br_prozora++;
	}
	if(br_prozora >= mfcc_velicina){
		br_prozora = 0;
		z2 = 1;
	}
	br++;
	
	audio_OUT = ((audio_chL <<16 & 0xFFFF0000)) + (audio_chL & 0x0000FFFF);	//Put the two channels toguether again
	i2s_tx(audio_OUT);
}

/*
void diskretizuj(){
	for(i = 0; i < BROJ_UZORAKA; i++) {
		
	}
	
	return;
}
*/
int main(void)
{
	arm_fir_init_f32(&S,N1,h,state,1);
	gpio_set_mode(P15, Output);
  audio_init (hz8000, line_in, intr, I2S_HANDLER);
  while(1){		
		if(z == 1)
		{
			z = 0;
			gpio_set(P15, HIGH);
			
			FFT(ulaz, BROJ_UZORAKA, 0.000125);
			if(z2 == 0) {
				mfcc1[br_prozora] = GetCoefficient(izlaz, 8000, 48, BROJ_UZORAKA, 1);
				mfcc2[br_prozora] = GetCoefficient(izlaz, 8000, 48, BROJ_UZORAKA, 2);
				mfcc3[br_prozora] = GetCoefficient(izlaz, 8000, 48, BROJ_UZORAKA, 3);
				mfcc4[br_prozora] = GetCoefficient(izlaz, 8000, 48, BROJ_UZORAKA, 4);
				mfcc5[br_prozora] = GetCoefficient(izlaz, 8000, 48, BROJ_UZORAKA, 5);
			}
			else{
				arm_correlate_f32(pirati1, H_duzina, &(mfcc1[2]), mfcc_velicina, cor1);
				arm_max_f32(cor1, 2 * H_duzina - 1, &(niz_max_cor1[0]), &pos_cor1);
				
				arm_correlate_f32(pirati2, H_duzina, &(mfcc2[2]), mfcc_velicina, cor2);
				arm_max_f32(cor2, 2 * H_duzina - 1, &(niz_max_cor2[0]), &pos_cor2);
				
				arm_correlate_f32(pirati3, H_duzina, &(mfcc3[2]), mfcc_velicina, cor3);
				arm_max_f32(cor3, 2 * H_duzina - 1, &(niz_max_cor3[0]), &pos_cor3);
				
				arm_correlate_f32(pirati4, H_duzina, &(mfcc4[2]), mfcc_velicina, cor4);
				arm_max_f32(cor4, 2 * H_duzina - 1, &(niz_max_cor4[0]), &pos_cor4);
				
				arm_correlate_f32(pirati5, H_duzina, &(mfcc5[2]), mfcc_velicina, cor5);
				arm_max_f32(cor5, 2 * H_duzina - 1, &(niz_max_cor5[0]), &pos_cor5);
				
				
				
				arm_correlate_f32(acdc1, H_duzina, &(mfcc1[2]), mfcc_velicina, cor1);
				arm_max_f32(cor1, 2 * H_duzina - 1, &(niz_max_cor1[1]), &pos_cor1);
				
				arm_correlate_f32(acdc2, H_duzina, &(mfcc2[2]), mfcc_velicina, cor2);
				arm_max_f32(cor2, 2 * H_duzina - 1, &(niz_max_cor2[1]), &pos_cor2);
				
				arm_correlate_f32(acdc3, H_duzina, &(mfcc3[2]), mfcc_velicina, cor3);
				arm_max_f32(cor3, 2 * H_duzina - 1, &(niz_max_cor3[1]), &pos_cor3);
				
				arm_correlate_f32(acdc4, H_duzina, &(mfcc4[2]), mfcc_velicina, cor4);
				arm_max_f32(cor4, 2 * H_duzina - 1, &(niz_max_cor4[1]), &pos_cor4);
				
				arm_correlate_f32(acdc5, H_duzina, &(mfcc5[2]), mfcc_velicina, cor5);
				arm_max_f32(cor5, 2 * H_duzina - 1, &(niz_max_cor5[1]), &pos_cor5);
				
				
				
				arm_correlate_f32(rasputin1, H_duzina, &(mfcc1[2]), mfcc_velicina, cor1);
				arm_max_f32(cor1, 2 * H_duzina - 1, &(niz_max_cor1[2]), &pos_cor1);
				
				arm_correlate_f32(rasputin2, H_duzina, &(mfcc2[2]), mfcc_velicina, cor2);
				arm_max_f32(cor2, 2 * H_duzina - 1, &(niz_max_cor2[2]), &pos_cor2);
				
				arm_correlate_f32(rasputin3, H_duzina, &(mfcc3[2]), mfcc_velicina, cor3);
				arm_max_f32(cor3, 2 * H_duzina - 1, &(niz_max_cor3[2]), &pos_cor3);
				
				arm_correlate_f32(rasputin4, H_duzina, &(mfcc4[2]), mfcc_velicina, cor4);
				arm_max_f32(cor4, 2 * H_duzina - 1, &(niz_max_cor4[2]), &pos_cor4);
				
				arm_correlate_f32(rasputin5, H_duzina, &(mfcc5[2]), mfcc_velicina, cor5);
				arm_max_f32(cor5, 2 * H_duzina - 1, &(niz_max_cor5[2]), &pos_cor5);
				
				
				
				arm_correlate_f32(mario1, H_duzina, &(mfcc1[2]), mfcc_velicina, cor1);
				arm_max_f32(cor1, 2 * H_duzina - 1, &(niz_max_cor1[3]), &pos_cor1);
				
				arm_correlate_f32(mario2, H_duzina, &(mfcc2[2]), mfcc_velicina, cor2);
				arm_max_f32(cor2, 2 * H_duzina - 1, &(niz_max_cor2[3]), &pos_cor2);
				
				arm_correlate_f32(mario3, H_duzina, &(mfcc3[2]), mfcc_velicina, cor3);
				arm_max_f32(cor3, 2 * H_duzina - 1, &(niz_max_cor3[3]), &pos_cor3);
				
				arm_correlate_f32(mario4, H_duzina, &(mfcc4[2]), mfcc_velicina, cor4);
				arm_max_f32(cor4, 2 * H_duzina - 1, &(niz_max_cor4[3]), &pos_cor4);
				
				arm_correlate_f32(mario5, H_duzina, &(mfcc5[2]), mfcc_velicina, cor5);
				arm_max_f32(cor5, 2 * H_duzina - 1, &(niz_max_cor5[3]), &pos_cor5);
				
				
				
				arm_correlate_f32(tetris1, H_duzina, &(mfcc1[2]), mfcc_velicina, cor1);
				arm_max_f32(cor1, 2 * H_duzina - 1, &(niz_max_cor1[4]), &pos_cor1);
				
				arm_correlate_f32(tetris2, H_duzina, &(mfcc2[2]), mfcc_velicina, cor2);
				arm_max_f32(cor2, 2 * H_duzina - 1, &(niz_max_cor2[4]), &pos_cor2);
				
				arm_correlate_f32(tetris3, H_duzina, &(mfcc3[2]), mfcc_velicina, cor3);
				arm_max_f32(cor3, 2 * H_duzina - 1, &(niz_max_cor3[4]), &pos_cor3);
				
				arm_correlate_f32(tetris4, H_duzina, &(mfcc4[2]), mfcc_velicina, cor4);
				arm_max_f32(cor4, 2 * H_duzina - 1, &(niz_max_cor4[4]), &pos_cor4);
				
				arm_correlate_f32(tetris5, H_duzina, &(mfcc5[2]), mfcc_velicina, cor5);
				arm_max_f32(cor5, 2 * H_duzina - 1, &(niz_max_cor5[4]), &pos_cor5);
				
				
				
				arm_correlate_f32(dp1, H_duzina, &(mfcc1[2]), mfcc_velicina, cor1);
				arm_max_f32(cor1, 2 * H_duzina - 1, &(niz_max_cor1[5]), &pos_cor1);
				
				arm_correlate_f32(dp2, H_duzina, &(mfcc2[2]), mfcc_velicina, cor2);
				arm_max_f32(cor2, 2 * H_duzina - 1, &(niz_max_cor2[5]), &pos_cor2);
				
				arm_correlate_f32(dp3, H_duzina, &(mfcc3[2]), mfcc_velicina, cor3);
				arm_max_f32(cor3, 2 * H_duzina - 1, &(niz_max_cor3[5]), &pos_cor3);
				
				arm_correlate_f32(dp4, H_duzina, &(mfcc4[2]), mfcc_velicina, cor4);
				arm_max_f32(cor4, 2 * H_duzina - 1, &(niz_max_cor4[5]), &pos_cor4);
				
				arm_correlate_f32(dp5, H_duzina, &(mfcc5[2]), mfcc_velicina, cor5);
				arm_max_f32(cor5, 2 * H_duzina - 1, &(niz_max_cor5[5]), &pos_cor5);
				
				for(i = 0; i < BROJ_PJESAMA; i++) {
					niz_max[i] = niz_max_cor1[i] + niz_max_cor2[i] + niz_max_cor3[i] + niz_max_cor4[i] + niz_max_cor5[i];
				}
				
 				z2 = 0;
			}
			
			gpio_set(P15, LOW);
		}
			
			//gpio_set(P15, z);
	}
}
/*
SAVE tetris1.dat 0x1FFDA240, 0x1FFDABDC
SAVE tetris2.dat 0x1FFDAC0C, 0x1FFDB5BC
SAVE tetris3.dat 0x1FFDB5D8, 0x1FFDBF88
SAVE tetris4.dat 0x1FFDBFA4, 0x1FFDC954 
SAVE tetris5.dat 0x1FFDC970, 0x1FFDD320
*/