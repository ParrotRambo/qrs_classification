#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wfdb/wfdb.h>
#include <wfdb/ecgmap.h>

typedef enum {false, true} bool;

#define K 20
#define FS 360
#define BUFFSIZE 600
#define MAX_LEADS 12

#define M 5
#define N 1
#define BETA 0.02

static WFDB_Siginfo s[MAX_LEADS];

unsigned stklen = 40000;
long nrSamps;
double *sig0;
double *sig4;

int *qrs;
int *qrsClass;
int currentQRS;
int totalQRS;

int m;
double dsum, dmin, dmax;
int openLeads;
long samples;
double normCnst = 32;

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

bool isQRS(int n)
{
  for(int i = 0; i < totalQRS; i++)
    if(qrs[i] == n)
      return true;
  return false;
}

double classifyQRS()
{
  double signal[BUFFSIZE], dcblock[BUFFSIZE], y[BUFFSIZE], yb[BUFFSIZE], zb[BUFFSIZE], f[M], x[M];
  int current;
  long unsigned int i, j, sample = 0;
  double threshold = 0;

  // Initialize vector of weighting or filter coefficients
  for(i = 0; i < M; i++)
  {
    if(i >= N && i < M - N)
      f[i] = 1. / (M - 2 * N);
    else
      f[i] = 0;
  }

  // Main loop
  do
  {
    if (sample >= BUFFSIZE)
		{
			for (i = 0; i < BUFFSIZE - 1; i++)
			{
				signal[i] = signal[i+1];
				dcblock[i] = dcblock[i+1];
				y[i] = y[i+1];
        yb[i] = yb[i+1];
				zb[i] = zb[i+1];
			}
			current = BUFFSIZE - 1;
		}
		else
		{
			current = sample;
		}

    if (sample >= samples)
			break;

		signal[current] = sig0[sample];
		sample++; // Update sample counter

    // DC Block filter
		if (current >= 1)
			dcblock[current] = signal[current] - signal[current-1] + 0.995*dcblock[current-1];
		else
			dcblock[current] = 0;

    // Highpass filter
    if(sample >= M)
    {
      y[current] = 0;

      for(i = 0; i < M; i++)
        x[i] = dcblock[current - i];
      qsort(x, M, sizeof(double), cmpfunc);

      for(i = 0; i < M; i++)
        y[current] += f[i] * x[i];

      yb[current] = dcblock[current - ((M + 1) / 2)] - y[current];
    }
    else
    {
      y[current] = 0;
      yb[current] = 0;
    }

    // Lowpass filter
    zb[current] = 0;
    for(i = 0; i < K && i < sample; i++)
      zb[current] += yb[current - i] * yb[current - i];
    //if(sample < 50)
    //  printf("%ld %f\n", i, zb[current]);

    // If the current sample is QRS activate thresholding mechanism
    if(isQRS(sample))
    {
      // If the current sample larger than the threshold, then this is a Normal QRS
      // Update the threshold
      if(zb[current] >= threshold)
      {
        qrsClass[currentQRS] = 1;
        threshold = BETA * 0.4 * zb[current] + (1 - BETA) * threshold;
      }
      // If the current sample smaller than the threshold, then this is not a Normal QRS
      else
      {
        qrsClass[currentQRS] = 0;
        //printf("%ld\n", sample);
      }
      currentQRS++;
    }

  } while(sample < samples);
}

void readFile(char *record)
{
  char *recordcpy = NULL;
  recordcpy = (char*) malloc(sizeof(char) * (strlen(record) + 2));
  strcpy(recordcpy, record);
  FILE *fp = fopen(strcat(recordcpy, ".txt"), "r");
  char ch;

  totalQRS = 0;
  currentQRS = 0;

  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      totalQRS++;
    }
  }

  qrs = (int*) malloc(sizeof(int) * totalQRS);
  qrsClass = (int*) malloc(sizeof(int) * totalQRS);

  fseek(fp, 0, SEEK_SET);
  for(int i = 0; i < totalQRS; i++)
    fscanf(fp, "%*s %d %*c %*d %*d %*d", &qrs[i]);

  fclose(fp);
}

int openRecord(char *record)
{
  return (isigopen(record, s, MAX_LEADS));
}

long ReadBuffer(long annot)
{
  long i, j;
  WFDB_Sample vec[MAX_LEADS];

  // A global variable s of type WFDB_Siginfo contains
  // all signal information and is initialized using the
  // isigopen (...) call (refer to function openRecord
  // of this frame).
  // Belos is an example of displaying the signal description,
  // gain and baseline value for each signal in the opened
  // signal file.
  for (j = 0; j < openLeads; j++)
    fprintf(stderr, "%s %lf %d\n", s[j].desc, s[j].gain, s[j].baseline);

  for (i = 0; i < nrSamps; i++)
  {
    if (getvec(vec) < 0)
      break;
    sig0[i] = 1000 * (vec[0] - s[0].baseline) / s[0].gain;
    /* Transfer of signal into uV (this equals to multiplying signal
       by 5), omit multiplication by 1000 to convert signal into mV. */
  }
  return (nrSamps);
}

void writeQRS(char *record, int chan, char *ann)
{
  /*
  WFDB_Anninfo annIFO;
  WFDB_Annotation annot;
  int i;

  annIFO.name = ann;
  annIFO.stat = WFDB_WRITE;
  if (annopen(record, & annIFO, 1) < 0)
  {
    fprintf(stderr, "Error opening QRS file\n");
    return;
  }
  annot.subtyp = annot.chan = annot.num = 0;
  annot.aux = NULL;
  for (i = 0; i < samples; i++)
  {
    if (sig4[i] != 0) {
      annot.anntyp = sig4[i];
      annot.time = i;
      if (putann(0, & annot) < 0) break;
    }
  }
  */
  char *recordcpy = NULL;
  recordcpy = (char*) malloc(sizeof(char) * (strlen(record) + 2));
  strcpy(recordcpy, record);
  FILE *fp = fopen(strcat(recordcpy, ".cls"), "w");

  for(int i = 0; i < totalQRS; i++)
  {
    char class;
    if(qrsClass[i] == 1)
      class = 'N';
    else
      class = 'V';

    fprintf(fp, "0:00:00.00 %d %c 0 0 0\n", qrs[i], class);
  }

  fclose(fp);
}

void readQRS(char *record, char *ant, int chan)
{
  WFDB_Anninfo annIFO;
  WFDB_Annotation annot;
  long i;

  annIFO.name = ant;
  annIFO.stat = WFDB_READ;

  if (annopen(record, & annIFO, 1) < 0)
  {
    fprintf(stderr, "Error opening QRS file\n");
    return;
  }

  annot.subtyp = annot.chan = annot.num = 0;
  annot.aux = NULL;

  for (i = 0; i < samples; i++)
    sig4[i] = 0;

  i = 0;
  while (getann(0, & annot) == 0)
  {
    if (annot.time > samples)
    {
      fprintf(stderr, "Error reading annotation times\n");
      return;
    }
    sig4[annot.time] = annot.anntyp;
  }
}

int main(int argc, char *argv[])
{
  long i;
  char *record = NULL;
  char *annotator = NULL;
  double thresh;

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] != '-')
    {
      fprintf(stderr, "Error parsing command line\n");
      exit(1);
    }

    switch (argv[i][1])
    {
    case 'r':
      i++;
      record = (char*) malloc(sizeof(char) * (strlen(argv[i]) + 2));
      strcpy(record, argv[i]);
      break;
    case 'a':
      i++;
      annotator = (char*) malloc(sizeof(char) * (strlen(argv[i]) + 2));
      strcpy(annotator, argv[i]);
      break;
    case 'n':
      i++;
      normCnst = atof(argv[i]);
      break;
    default:
      fprintf(stderr, "Wrong switch\n");
      exit(2);
    }
  }

  if (record == NULL)
  {
    fprintf(stderr, "No record was specified, exiting\n");
    exit(2);
  }

  if ((openLeads = openRecord(record)) < 0)
  {
    fprintf(stderr, "Error opening record, exiting\n");
    exit(3);
  }

  fprintf(stderr, "Record opened\n");
  nrSamps = s -> nsamp - 1;
  sig0 = (double*) malloc(sizeof(double) * nrSamps);
  sig4 = (double*) malloc(sizeof(double) * nrSamps);

  if ((samples = ReadBuffer(0)) <= 0)
  {
    fprintf(stderr, "Error opening record, exiting\n");
    exit(3);
  }

  fprintf(stderr, "Data read\n");
  if (samples != nrSamps)
    fprintf(stderr, "Sample count differs\n");

  readFile(record);
  classifyQRS();
  fprintf(stderr, "Data analyzed\n");
  if (annotator == NULL)
  {
    annotator = (char*) malloc(sizeof(char) * 5);
    sprintf(annotator, "qrs");
  }

  writeQRS(record, 4, annotator);
  fprintf(stderr, "Annotations written\n");
  wfdbquit();
  fprintf(stderr, "Record closed\n");
  free(sig0);
  free(sig4);
  free(annotator);
  return 0;
}
