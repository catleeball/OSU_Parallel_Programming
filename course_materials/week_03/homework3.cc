// OSU cs575 project 3:
//   - https://web.engr.oregonstate.edu//~mjb/cs575/Projects/proj03.html
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES

#ifndef DEBUG
#define DEBUG	false
#endif


// Printed stats
int	  NowYear;	  // 2021 - 2026
int  	NowMonth;	  // 0 - 11
float	NowPrecip;	// inches of rain per month
float	NowTemp;		// temperature this month
float	NowHeight;	// grain height in inches
int	  NowNumDeer;	// number of deer in the current population

// Const vals for simulation.
const float GRAIN_GROWS_PER_MONTH   = 9.0;
const float ONE_DEER_EATS_PER_MONTH = 1.0;
const float AVG_PRECIP_PER_MONTH    = 7.0;	// average
const float AMP_PRECIP_PER_MONTH    = 6.0;	// plus or minus
const float RANDOM_PRECIP           = 2.0;	// plus or minus noise
const float AVG_TEMP                = 60.0;	// average
const float AMP_TEMP                =	20.0;	// plus or minus
const float RANDOM_TEMP             = 10.0;	// plus or minus noise
const float MIDTEMP                 = 40.0;
const float MIDPRECIP               = 10.0;



unsigned int seed = 0;

float Ranf( unsigned int *seedp,  float low, float high) {
  float r = (float) rand_r( seedp ); // 0 - RAND_MAX
  return(low + r * ( high - low ) / (float)RAND_MAX);
}

int Ranf( unsigned int *seedp, int ilow, int ihigh ) {
  float low  = (float)ilow;
  float high = (float)ihigh + 0.9999f;
  return (int)(Ranf(seedp, low,high));
}

// float x = Ranf( &seed, -1.f, 1.f );



void weather() {
  float ang = (  30.*(float)NowMonth + 15.  ) * ( M_PI / 180. );
  float temp = AVG_TEMP - AMP_TEMP * cos( ang );
  NowTemp = temp + Ranf(&seed, -RANDOM_TEMP, RANDOM_TEMP);

  float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin( ang );
  NowPrecip = precip + Ranf(&seed, -RANDOM_PRECIP, RANDOM_PRECIP);
  if( NowPrecip < 0. )
	  NowPrecip = 0.;
}

void iterMonth() {
  // Advance time one month, calculate the monthly weather. 
  if () //TODO: here's where I left off
  NowMonth++;
  NowYear++;
  weather();
}

void deer() {
  while (NowYear < 2027) {
    if (DEBUG) fprintf(stderr, "NowYear is: %d", NowYear);
    // Spacing just for funzies.
    int      nextNumDeer       = NowNumDeer                         ;
    int      carryingCapacity  = (int)(NowHeight )                  ;
    if      (nextNumDeer       < carryingCapacity) { nextNumDeer++  ; }
    else if (nextNumDeer       > carryingCapacity) { nextNumDeer--  ; }
    if      (nextNumDeer       < 0               ) { nextNumDeer = 0; }

    // DoneComputing
    #pragma omp barrier

    if (DEBUG) fprintf(stderr, "nextNumDeer is: %d", nextNumDeer);
    NowNumDeer = nextNumDeer;

    // DoneAssigning
    #pragma omp barrier

    if (DEBUG) fprintf(stderr, "nextNumDeer is: %d", nextNumDeer);
    // DonePrinting
    #pragma omp barrier
  }
}


void grain() {
  while (NowYear < 2027) {
    // TODO: Figure out why factors as functions didnt work
    float tempFactor   = exp(-1 * pow((NowTemp   - MIDTEMP), 2) / 10.);
    float precipFactor = exp(-1 * pow((NowPrecip - MIDTEMP), 2) / 10.);
    float nextHeight   = NowHeight;
    nextHeight += tempFactor * precipFactor * GRAIN_GROWS_PER_MONTH;
    nextHeight -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
    if (nextHeight < 0.) nextHeight = 0.;


    if (DEBUG) fprintf(stderr, "nextHeight is: %.f", nextHeight);

    // DoneComputing
    #pragma omp barrier

    if (DEBUG) fprintf(stderr, "nextHeight is: %.f", nextHeight);

    NowHeight = nextHeight;

    // Done Assigning
    #pragma omp barrier

    // Done Printing
    #pragma omp barrier
  }
}

void watcher() {
  while (NowYear < 2027) {
    // DoneComputing
    #pragma omp barrier
    // DoneAssigning
    #pragma omp barrier
    printf("%-d\t%-2d\t%-.5f\t%-.5f\t%-.5f\t%-d\t\n",
           NowYear, NowMonth, NowTemp, NowPrecip, NowHeight, NowNumDeer);

    // DonePrinting
    #pragma omp barrier
    NowYear++;
  }
}

// void myAgent() {}

void init() {
  NowYear    = 2021;
  NowMonth   = 0;
  NowNumDeer = 1;
  NowHeight  = 1.;
  weather();
}

int main(int argc, const char** argv) {

  #ifndef _OPENMP
  	fprintf(stderr, "No OpenMP support!\n");
	  return 1;
  #endif

  init();
  // Print initial state of simulation.
  printf("%-d\t%-2d\t%-.5f\t%-.5f\t%-.5f\t%-d\t\n",
         NowYear, NowMonth, NowTemp, NowPrecip, NowHeight, NowNumDeer);

  omp_set_num_threads(3);
  #pragma omp parallel sections //private(NowHeight,NowNumDeer,NowPrecip,NowTemp,NowMonth,NowYear),shared(stderr)
  {
    // Each section/thread has three barriers: DoneComputing, DoneAssigning,
    // and DonePrinting. Details on each barrier in project 3 webpage under
    // "Use of Threads", and in project 3 slides on slide 5.
	  #pragma omp section
	  { deer(); }

	  #pragma omp section
  	{ grain(); }

	  #pragma omp section
  	{ watcher(); }

//  	#pragma omp section
//	  { myAgent(); }
  }
  // implied barrier, all funcs must return for any of them to get past here

  return 0;
}

// TODO: Imagine scenario where you add another thread worker, but don't put
// exactly 3 barriers in it, or make a mistake and one thread doesn't do all
// the work it needs to have done by a given barrier. State of each thread will
// de-sync and cause confusing behavior and a troublesome bug to catch.
//
// Consider mutex lock to synchronize? Have watcher validate thread states?
// Have threads do state validation when they hit the barrier?

