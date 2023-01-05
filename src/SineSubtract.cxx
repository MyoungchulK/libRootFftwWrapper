#ifndef DONT_HAVE_MINUIT2
#include "SineSubtract.h" 


#include "TGraph.h" 
#include "TCanvas.h" 
#include "TStyle.h" 
#include "TMutex.h" 
#include "TLinearFitter.h" 
#include "DigitalFilter.h" 
#include "TSpectrum.h" 

#ifndef __APPLE__
#include <malloc.h>
#define SINCOS sincos 
#else
#include <malloc/malloc.h>
#define SINCOS __sincos
#endif

#ifdef SINE_SUBTRACT_PROFILE
#include "TStopwatch.h"
#endif

#include "TMath.h"
#include "FFTtools.h"
#include "TF1.h" 
#include "TH2.h"

#ifdef SINE_SUBTRACT_USE_FLOATS
#define VEC_T float
#else
#define VEC_T double
#endif

#ifdef ENABLE_VECTORIZE
#include "vectormath_trig.h" 

//assume AVX is present, otherwise it'll be simulated 

#ifdef SINE_SUBTRACT_USE_FLOATS
#define VEC Vec8f 
#define VEC_N 8
#else
#define VEC Vec4d 
#define VEC_N 4 
#endif


#endif 



static TMutex fnLock; 
static int fnCount = 0; 
static int counter; 




#include "TError.h" 
//extern int gErrorIgnoreLevel; // who ordered that? 

static double guessPhase(const TGraph * g, double freq) 
{
  double phase; 
  FFTtools::dftAtFreq(g,freq,&phase,0); 
  return phase; 
}



static double normalize_angle(double phi)
{
  return  phi - 2*TMath::Pi() * FFTtools::fast_floor((phi+ TMath::Pi()) / (2.*TMath::Pi())); 
}


static __thread TLinearFitter * fitter = 0; 
static __thread int fitter_order = 0; 



static const char * order_strings[] = { "1", "1++x", "1++x++x*x","1++x++x*x++x*x*x","1++x++x*x++x*x*x++x*x*x*x" }; 


static void computeEnvelopeFit(int order, int N, const double * x, const double * y, double * out) 
{

  if (!fitter || fitter_order != order) 
  {
    fnLock.Lock(); 
    if (!fitter) 
    {
      fitter = new TLinearFitter(1,order_strings[order],""); 
    }
    else
    {
      fitter->SetFormula(order_strings[order]); 
    }
    fnLock.UnLock(); 
  }

  fitter->ClearPoints(); 
  fitter->AssignData(N,1,(double*)x,(double*)y); 
  fitter->Eval(); 

  double max = 0; 
  for (int i = 0; i < N; i++)
  {
    out[i] = fitter->GetParameter(0);
    double xx = 1; 
    for (int j = 1; j <= order; j++) 
    {
      xx *= x[i]; 
      out[i] += fitter->GetParameter(j) * xx; 
    }
    if (out[i] > max) max = out[i]; 
  }

  for (int i = 0; i < N; i++) 
  {
    out[i]/=max; 
  }
}

/* wrapper around hilbert envelope because need to interpolate to calculate it properly */
static void computeEnvelopeHilbert(const TGraph * g, double * out, double dt) 
{

  TGraph * ig = FFTtools::getInterpolatedGraph(g, dt); 
  TGraph * hilbert = FFTtools::getHilbertEnvelope(ig); 
  double xmax = ig->GetX()[ig->GetN()-1]; 

  for (int i = 0; i < g->GetN(); i++)
  {
    if (g->GetX()[i] < xmax)
    {
      out[i] = FFTtools::evalEvenGraph(hilbert, g->GetX()[i]); 
    }
    else
    {
      out[i] = ig->GetY()[ig->GetN()-1] ; 
    }
  }

  delete ig; 
  delete hilbert; 
}



FFTtools::SineFitter::SineFitter() : fDoEvalRecord(false), f(this)
{
  min.SetFunction(f);
  verbose = false; 

}


void FFTtools::SineFitter::deleteEvalRecords(){
  for(unsigned i=0; i < grEvalRecords.size(); i++){
    if(grEvalRecords[i]){
      delete grEvalRecords[i];
      grEvalRecords[i] = NULL;
    }
  }
  grEvalRecords.clear();
}

FFTtools::SineFitter::~SineFitter()
{

  deleteEvalRecords();
}

void FFTtools::SineFitter::setGuess(double fg, int ntrace, const double * phg, double ampg)
{
  freq = fg; 
  phase.clear(); 
  amp.clear(); 
  phase_err.clear(); 
  amp_err.clear(); 

  phase.insert(phase.end(), phg, phg + ntrace); 
  amp.insert(amp.end(), ntrace, ampg); 
  phase_err.insert(phase_err.end(), ntrace, 0); 
  amp_err.insert(amp_err.end(), ntrace, 0); 
}

ROOT::Math::IBaseFunctionMultiDim* FFTtools::SineFitter::SineFitFn::Clone() const
{
  SineFitFn * fn = new SineFitFn; 
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
  fn->setXY(nt,ns,xp,yp,wgt,envp); //incredibly inefficient, but quick kludge for now. We're probalby never going to clone anyway. 
#else
  fn->setXY(nt,ns,x,y,wgt,env); 
#endif
  fn->fContainer = fContainer;
  return fn; 
}

FFTtools::SineFitter::SineFitFn::SineFitFn(SineFitter* parent)
    : ROOT::Math::IGradientFunctionMultiDim(), fContainer(parent)
{
  ns = NULL; 
  nt = 0; 
  x= 0; 
  y = 0; 
  env = 0; 
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
  xp = 0; 
  yp = 0; 
  envp = 0; 
#endif

}


//here's hoping... 
__attribute__((optimize("unroll-loops")))
double FFTtools::SineFitter::SineFitFn::DoDerivative(const double * p, unsigned int coord) const
{

  double w = 2 *TMath::Pi() * p[0]; 
  double deriv = 0; 
  int type = coord == 0 ? 0 : (1 + ((coord-1) % 2)); 

  for (int ti = 0; ti < nt; ti++) 
  {

    if (type > 0 && (int(coord) - 1)/2 != ti) continue; 
    __builtin_prefetch(x[ti]); 
    __builtin_prefetch(y[ti]); 
    double ph = normalize_angle(p[1 + 2*ti]); 
    double A = p[2+2*ti]; 
    if (wgt) A*=wgt[ti]; 

#ifdef ENABLE_VECTORIZE 

    VEC vecx(0); 
    VEC vecy(0); 
    VEC vece(1); 
    VEC vec_ph(ph); 
    VEC vec_w(w); 
    VEC vec_A(A); 
    VEC vec_cos(0); 
    VEC dYdp(0); 
    int leftover = ns[ti] % VEC_N;
    int nit = ns[ti]/VEC_N + (leftover ? 1 : 0); 

#ifdef SINE_SUBTRACT_FORCE_ALIGNED
    double * xx = (double*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    double * yy = (double*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    double * envelope =   (env && env[ti]) ?  (double*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#elif defined(SINE_SUBTRACT_USE_FLOATS)
    float * xx = (float*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    float * yy = (float*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    float * envelope =   (env && env[ti]) ?  (float*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#else
    double * xx = (double*) x[ti]; 
    double * yy = (double*) y[ti]; 
    double * envelope =   (env && env[ti]) ? (double*) env[ti] : 0; 
#endif
    

    for (int i = 0; i < nit; i++)
    {
      if (i < nit-1 || !leftover)
      {
#if defined(SINE_SUBTRACT_FORCE_ALIGNED) || defined(SINE_SUBTRACT_USE_FLOATS)
        vecx.load_a(xx+VEC_N*i); 
        vecy.load(yy+VEC_N*i); 
#else
        vecx.load(xx+VEC_N*i); 
        vecy.load(yy+VEC_N*i); 
#endif
      }
      else
      {
        vecx.load_partial(leftover, xx+VEC_N*i); 
        vecy.load_partial(leftover, yy+VEC_N*i); 
      }

      VEC vec_ang = mul_add(vecx, vec_w, vec_ph); 
      VEC vec_sin = sincos(&vec_cos, vec_ang); 

      if (envelope)
      {

        if (i < nit -1 || !leftover)
        {
          vece.load(envelope + VEC_N*i); 
        }
        else
        {
          vece.load_partial(leftover, envelope + VEC_N*i); 
        }
        vec_sin *= vece; 
        vec_cos *= vece; 
      }


      VEC vec_Y = mul_sub(vec_sin, vec_A, vecy); 

      switch (type) 
      {
        case 0: 
          dYdp = vec_A *vecx*vec_cos * (2 * M_PI); 
          break; 
        case 1: 
          dYdp = vec_A *vec_cos; 
          break; 
        case 2: 
          dYdp = vec_sin; 
          break; 
      }


#ifndef SINE_SUBTRACT_DONT_HORIZONTAL_ADD
      if (i == nit-1 && leftover) //hopefully this gets unrolled? 
      {
        vec_Y.cutoff(leftover); 
      }

      deriv +=horizontal_add( vec_Y * dYdp)/ns[ti]; 
#else

      VEC ans = vec_Y * dYdp; 
      VEC_T ans_v[VEC_N] __attribute__((aligned(sizeof(VEC))); 
      ans.store_a(ans_v); 

      int vecn = (i == nit-1 && leftover) ? leftover : VEC_N;
      for (int j = 0; j < vecn; j++) 
      {
        deriv += ans_v[j]/ns[ti]; 
      }

#endif
    }

#else
    for (int i = 0; i < ns[ti]; i++)
    {
      double t = x[ti][i]; 
      double sinang = sin(w*t + ph); 
      double Y = A *sinang; 
      if (env && env[ti]) Y*= env[ti][i]; 

      double dYdp = 0; 

      switch(type) 
      {
        case 0: 
          dYdp = A*t*cos(w*t + ph) * 2 * TMath::Pi(); 
          break; 
        case 1: 
          dYdp = A*cos(w*t + ph); 
          break; 
        case 2:
          dYdp = sinang; 
          break; 
      }
      if (env && env[ti]) sinang*= env[ti][i]; 

      deriv += ((Y-y[ti][i]) * dYdp)/ns[ti]; 
    }
#endif
  }

  deriv *=2; 

//  printf("dP/d%s (f=%f, ph=%f, A=%f) = %f\n", coord == 0 ? "f" : coord == 1? "ph" : "A", p[0], ph, A, deriv);  
  return deriv/nt; 
}



// yeah, like the compiler is really going to do what I want ... may have to unroll this manually 
__attribute__((optimize("unroll-loops")))
double FFTtools::SineFitter::SineFitFn::DoEval(const double * p) const
{

  VEC_T w = 2 *TMath::Pi() * p[0]; 
  VEC_T power = 0; 


  for (int ti = 0; ti < nt; ti++)
  {
    VEC_T ph = normalize_angle(p[1+2*ti]); 
    VEC_T A = p[2+2*ti]; 
    if (wgt) A*= wgt[ti]; 

#ifdef ENABLE_VECTORIZE 
    VEC vecx(0); 
    VEC vecy(0); 
    VEC vece(1); 
    VEC vec_ph(ph); 
    VEC vec_w(w); 
    VEC vec_A(A); 

#ifdef SINE_SUBTRACT_FORCE_ALIGNED
    double * xx = (double*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    double * yy = (double*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    double * envelope =   (env && env[ti]) ?  (double*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#elif defined(SINE_SUBTRACT_USE_FLOATS)
    float * xx = (float*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    float * yy = (float*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    float * envelope =   (env && env[ti]) ?  (float*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#else
    double * xx = (double*) x[ti]; 
    double * yy = (double*) y[ti]; 
    double * envelope =   (env && env[ti]) ?  (double*) env[ti] : 0; 
#endif


    int leftover = ns[ti] % VEC_N;
    int nit = ns[ti]/VEC_N + (leftover ? 1 : 0); 

    for (int i = 0; i < nit; i++)
    {
      if (i < nit -1 || !leftover)
      {
        vecx.load(xx+VEC_N*i); 
        vecy.load(yy+VEC_N*i); 
      }
      else
      {
        vecx.load_partial(leftover, xx+VEC_N*i); 
        vecy.load_partial(leftover, yy+VEC_N*i); 
      }


      VEC vec_ang = mul_add(vecx, vec_w, vec_ph); 
      VEC vec_sin = sin(vec_ang); 

      //TODO if branch-prediction hurts, we should have two versions of this loop 
      if (envelope)
      {

        if (i < nit -1 || !leftover)
        {
          vece.load(envelope + VEC_N*i); 
        }
        else
        {
          vece.load_partial(leftover, envelope + VEC_N*i); 
        }
        vec_sin *= vece; 
      }

      VEC vec_Y = mul_sub(vec_sin, vec_A, vecy); 

      VEC vec_Y2 = square(vec_Y); 
#ifndef SINE_SUBTRACT_DONT_HORIZONTAL_ADD
      if (i == nit-1 && leftover) //hopefully this gets unrolled? 
      {
        vec_Y2.cutoff(leftover); 
      }

      power += horizontal_add(vec_Y2)/ns[ti]; 
#else 
      int vecn = (i == nit-1 && leftover) ? leftover : VEC_N;
      for (int j = 0; j< vecn; j++)
      {
        power += vec_Y2[j]/ns[ti]; 
      }
#endif
    }
#else
    for (int i = 0; i < ns[ti]; i++)
    {
      VEC_T Y = A * sin(w*x[ti][i] + ph); 
      if (env && env[ti]) Y*= env[ti][i]; 
      Y-=y[ti][i]; 
      power += Y*Y/ns[ti]; 
    }
#endif 

  }


  double evalResult = power/nt;
  if(fContainer && fContainer->GetDoEvalRecord()){
    TGraph* grEval = fContainer->getEvalRecordGraph();
    if(grEval){
      int numEvals = grEval->GetN();
      grEval->SetPoint(numEvals, numEvals, evalResult);
    }
  }

//  printf("P(f=%f, ph=%f, A=%f) = %f\n", p[0], ph, A, power);

  return evalResult;

}

FFTtools::SineFitter::SineFitFn::~SineFitFn()
{
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
setXY(0,0,0,0); 
#endif

if (ns) delete [] ns; 

}

void FFTtools::SineFitter::SineFitFn::setXY(int ntraces, const int * nsamples, const double ** xx, const double ** yy, const double * wset, const double ** envset) 
{
  wgt = wset; 
#ifdef SINE_SUBTRACT_USE_FLOATS
  // we have to convert everything to floats... 
  if (nt!=ntraces)
  {

    float ** newx = (float**) malloc(ntraces * sizeof(float*)); 
    float ** newy = (float**) malloc(ntraces * sizeof(float*)); 
    float ** newenv = 0; 
    if (envset) 
      new_env = (float**)  emalloc(ntraces * sizeof(float*)); 
    for (int i = 0; i < nt; i++) 
    {
      if (i < ntraces)
      {
        newx[i] = x[i]; 
        newy[i] = y[i]; 
        if (env) 
          newenv[i] = env[i]; 
      }
      else
      {
        free(x[i]); 
        free(y[i]); 
        if (env) 
          free(env[i]); 
      }
    }

    for (int i = nt; i < ntraces; i++) 
    {
      newx[i] = 0; 
      newy[i] = 0; 
      if (newenv) 
        newenv[i] = 0; 
    }


    free(x); 
    free(y); 
    if (env) 
      free(env); 
    x = newx; 
    y = newy; 
    env = newenv; 
  }
   
  for (int i = 0; i < ntraces; i++) 
  {
    if (i >= nt || ns[i]!=nsamples[i])
    {
       if (x[i]) free(x[i]); 
       x[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(float)); 
       if (y[i]) free(y[i]); 
       y[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(float)); 
       if (env && env[i]) free(env[i]); 
       if (env) env[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(float)); 
    }
  }


  for (int j = 0; j < ntraces; j++) 
  {
    for (int i = 0; i < nsamples[j]; i++) 
    {
       x[j][i] = xx[j][i]; 
       y[j][i] = yy[j][i]; 
       if (env && env[j])
         env[j][i] = envset[j][i];  
    }
  }

  xp = xx; 
  yp = yy; 
  envp = envset; 
  
#elif SINE_SUBTRACT_FORCE_ALIGNED
  if (nt!=ntraces)
  {

    double ** newx = (double**) malloc(ntraces * sizeof(float*)); 
    double ** newy = (double**) malloc(ntraces * sizeof(float*)); 

    double ** newenv = 0;

    if (envset) 
      newenv = (double**) malloc(ntraces * sizeof(float*)); 

    for (int i = 0; i < nt; i++) 
    {
      if (i < ntraces)
      {
        newx[i] = x[i]; 
        newy[i] = y[i]; 
        if (env) newenv[i] = env[i]; 
      }
      else
      {
        free(x[i]); 
        free(y[i]); 
        if (env) free(env[i]); 
      }
    }

    for (int i = nt; i < ntraces; i++) 
    {
      newx[i] = 0; 
      newy[i] = 0; 
      if (newenv) newenv[i] = 0; 
    }


    free(x); 
    free(y); 
    if (env) 
      free(env); 
    x = newx; 
    y = newy; 

    env = newenv; 
  }

  for (int i = 0; i < ntraces; i++) 
  {
    if (i >= nt || ns[i]!=nsamples[i])
    {
       if (x[i]) free(x[i]); 
       x[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(double)); 
       if (y[i]) free(y[i]); 
       y[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(double)); 
       if (env && env[i]) free(env[i]); 
       if (env) env[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(double)); 
    }
  }



  for (int j = 0; j < ntraces; j++) 
  {
    memcpy(x[j], xx[j], nsamples[j] * sizeof(double)); 
    memcpy(y[j], yy[j], nsamples[j] * sizeof(double)); 
    if (envset && envset[j])
    {
      memcpy(env[j], envset[j], nsamples[j] * sizeof(double)); 
    }
  }

  xp = xx; 
  yp = yy; 
  envp = env; 
 
#else
  x = xx; 
  y = yy;
  env = envset; 

#endif

  if (nsamples)
  {
    if (!ns) 
    {
      ns = new int[ntraces]; 
    }
    else if (nt != ntraces)
    {
      delete [] ns; 
      ns = new int[ntraces];  
    }
    memcpy(ns,nsamples, ntraces*sizeof(int)); 
  }
  else 
  {
    if (ns) 
    {
      delete [] ns; 
    }

    ns=0; 
  }

  nt = ntraces; 
}


//stupid function that converts an array to a string 
static char * arrString(int n, const  double * arr, int sigfigs =8 ) 
{
  int size = 2 + (n-1) + (sigfigs+5) * n + 2;
  
  char * ret = new char[size]; 
  char format[8]; 
  sprintf(format,"%%.%dg",sigfigs); 

  int ctr = 0; 
  ctr += sprintf(ret + ctr, "("); 
  for (int i = 0; i < n; i++) 
  {
    ctr+=sprintf(ret + ctr, format, arr[i]); 
    if (i < n-1)  ctr += sprintf(ret + ctr, ","); 
  }
  ctr+= sprintf(ret+ctr, ")"); 

  assert(ctr < size); 

  return ret; 
}





void FFTtools::SineFitter::doFit(int ntraces, const int * nsamples, const double ** x, const double **y, const double * w, const double ** env)
{
  f.setXY(ntraces,nsamples,x,y,w,env); 

  if (verbose) 
  {
    char * Phstr =arrString(ntraces,&phase[0]); 
    printf("Guesses: f= %f, A=%f, ph=%s", freq, amp[0], Phstr); 
    delete [] Phstr; 
    double p[1 + 2*ntraces]; 
    p[0] = freq; 

    for (int i =0; i < ntraces; i++) 
    {
      p[1+2*i] = phase[i]; 
      p[2+2*i] = amp[0]; 
    }
    printf("Guess power is %f. Gradient is (", f.DoEval(p));
    for (int i = 0; i < 2 * ntraces+1; i++) 
    {
      printf("%f", f.DoDerivative(p,i)); 
      if (i < 2*ntraces) printf(","); 
    }
    printf(")\n"); 

  }

  //Estimate Nyquist to set frequency cut 
  // Right now this just uses the first graph, although it would probably be smarter
  // to use the average or something 
  
  double dt = (x[0][nsamples[0]-1]  - x[0][0]) / (nsamples[0]-1); 
  double fnyq = 1. / (2 * dt); 
  double df = fnyq / nsamples[0]; 



  min.Clear(); 

  min.SetFunction(f); 


  if (limits.max_n_df_relative_to_guess > 0)
  {
    min.SetLimitedVariable(0, "f",freq, df  * limits.freq_start_error, freq-df * limits.max_n_df_relative_to_guess, freq+df * limits.max_n_df_relative_to_guess); 
  }
  else if (limits.max_n_df_relative_to_guess == 0)
  {
    min.SetFixedVariable(0,"f",freq); 
  }
  else 
  {
    min.SetVariable(0, "f",freq, df * limits.freq_start_error); 
  }


  double damp = limits.amp_start_error* amp[0]; 


  for (int i = 0; i < ntraces; i++) 
  {
    min.SetVariable(1+2*i, TString::Format("phi%d",i).Data(),phase[i], limits.phase_start_error); 

    if (limits.maxA_relative_to_guess > 0 && limits.minA_relative_to_guess > 0)
    {
       min.SetLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.minA_relative_to_guess, amp[0]*limits.maxA_relative_to_guess); 
    }
    else if (limits.maxA_relative_to_guess > 0)
    {
       min.SetUpperLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.maxA_relative_to_guess); 
    }
    else if (limits.minA_relative_to_guess > 0)
    {
       min.SetLowerLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.minA_relative_to_guess);
    }
    else if (limits.minA_relative_to_guess < 0 && limits.maxA_relative_to_guess > 0)
    {
       min.SetVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp); 
    }
    else
    {
      min.SetFixedVariable(2 + 2 * i, TString::Format("A%d",i).Data(), amp[0]); 
    }
  }


  min.SetPrintLevel(verbose ? 1 : 0); 

  int old_level = gErrorIgnoreLevel; 
  if (!verbose)
  {
    gErrorIgnoreLevel = 1001; 
  }
  if(fDoEvalRecord){
    grEvalRecords.push_back(new TGraph());
  }
  min.Minimize(); 
  gErrorIgnoreLevel = old_level; 
  if(verbose)  min.PrintResults(); 

  freq = min.X()[0]; 
  freq_err = min.Errors()[0]; 
  for (int i = 0; i < ntraces; i++) 
  {
    phase[i] = normalize_angle(min.X()[1+2*i]); 
    amp[i] = min.X()[2+2*i]; 
    phase_err[i] = min.Errors()[1+2*i]; 
    amp_err[i] =  min.Errors()[2+2*i]; 
  }
}

FFTtools::SineSubtract::SineSubtract(const TGraph * gmp, int maxiter, bool store)
  : abs_maxiter(0), maxiter(maxiter), max_successful_iter(0),
  min_power_reduction(gmp->GetMean(2)), store(store),
  power_estimator_params(0), peak_option_params(0), envelope_option_params(0) 

{
  g_min_power = gmp; 
  setPowerSpectrumEstimator(FFT); 
  setPeakFindingOption(NEIGHBORFACTOR); 
  setEnvelopeOption(ENV_NONE); 
  high_factor = 1; 
  verbose = false; 
  tmin = 0; 
  tmax = 0; 
  nfail_exponent =0.5; 
  id = counter++; 
#ifdef ENABLE_VECTORIZE
  no_subnormals();  // Unlikely to occur, but worth it if they do
#endif

  
}



FFTtools::SineSubtract::SineSubtract(int fft_len, double base_fft[][16], int maxiter, double min_power_reduction, bool store)
//FFTtools::SineSubtract::SineSubtract(int maxiter, double min_power_reduction, bool store)
  : abs_maxiter(0), maxiter(maxiter), max_successful_iter(0),
    min_power_reduction(min_power_reduction), store(store),
    power_estimator_params(0), peak_option_params(0), envelope_option_params(0) 

{

  //! this should be dBm/Hz scale
  baseline_fft.resize(16);
  for (int j = 0; j < 16; j++) {
    for (int i = 0; i < fft_len; i++) {
      baseline_fft[j].push_back(10 * log10(base_fft[i][j])); 
    }
  }

  setPowerSpectrumEstimator(FFT); 
  setPeakFindingOption(NEIGHBORFACTOR); 
  setEnvelopeOption(ENV_NONE); 
  high_factor = 1; 
  verbose = false; 
  tmin = 0; 
  tmax = 0; 
  g_min_power =0; 
  nfail_exponent =0.5; 
  id = counter++; 
#ifdef ENABLE_VECTORIZE
  no_subnormals();  // Unlikely to occur, but worth it if they do
#endif
}

void FFTtools::SineSubtractResult::append(const SineSubtractResult *r) 
{
  powers.insert(powers.end(), r->powers.begin(), r->powers.end()); 
  freqs.insert(freqs.end(), r->freqs.begin(), r->freqs.end()); 
  freqs_errs.insert(freqs_errs.end(), r->freqs_errs.begin(), r->freqs_errs.end()); 

  for (size_t i = 0; i < phases.size(); i++) 
  {
    phases[i].insert(phases[i].end(), r->phases[i].begin(), r->phases[i].end()); 
    amps[i].insert(amps[i].end(), r->amps[i].begin(), r->amps[i].end()); 
    phases_errs[i].insert(phases_errs[i].end(), r->phases_errs[i].begin(), r->phases_errs[i].end()); 
    amps_errs[i].insert(amps_errs[i].end(), r->amps_errs[i].begin(), r->amps_errs[i].end()); 
  }
}


void FFTtools::SineSubtractResult::clear() 
{
  powers.clear(); 
  freqs.clear(); 
  phases.clear(); 
  amps.clear(); 
  freqs_errs.clear(); 
  phases_errs.clear(); 
  amps_errs.clear(); 
}


void FFTtools::SineSubtract::reset() 
{

  r.clear(); 
  for (unsigned i = 0; i < gs.size(); i++)
  {
    for (unsigned j = 0; j < gs[i].size(); j++) 
    {
      delete gs[i][j]; 
    }
    for (unsigned j = 0; j < env_gs[i].size(); j++) 
    {
      delete env_gs[i][j]; 
    }
  }
  gs.clear(); 
  env_gs.clear(); 
  for (unsigned i = 0; i < spectra.size(); i++) delete spectra[i]; 
  spectra.clear(); 


}
/*
TGraph * FFTtools::SineSubtract::subtractCW(const TGraph * g, double dt, const SineSubtractResult* result) 
{
  TGraph * gcopy = new TGraph(g->GetN(), g->GetX(), g->GetY()); 
  gcopy->SetTitle(TString::Format("%s (subtracted)", g->GetTitle())); 
  subtractCW(1,&gcopy,dt, NULL, result);
  return gcopy; 
}
*/

double * FFTtools::SineSubtract::subtractCW(int N, const double * y, double dt, int pad_len, int * bad_idxs, int bad_len, double thres, int ant, double * yout,  const SineSubtractResult* result) 
{
  TGraph g(N); 
  for (int i = 0; i < N;i ++) 
  {
    g.GetX()[i] = i*dt; 
    g.GetY()[i] = y[i]; 
  }
  TGraph *gptr = &g;
  subtractCW(1, &gptr, dt, pad_len, bad_idxs, bad_len, thres, ant, NULL, result); 
  if (!yout) 
  {
    yout = new double[N]; 
  }
  memcpy(yout, g.GetY(), N*sizeof(double)); 
  return yout; 
}

void FFTtools::SineSubtract::subtractCW(int ntraces, TGraph ** g, double dt, int pad_len, int * bad_idxs, int bad_len, double thres, int ant, const double * w, const SineSubtractResult* result) 
{


#ifdef SINE_SUBTRACT_PROFILE
  TStopwatch sw;
  fitter.SetDoEvalRecord(true);
#endif
  fitter.deleteEvalRecords();
  reset(); 

  //! t-domain configuration
  if (dt<=0) dt = g[0]->GetX()[1] - g[0]->GetX()[0];
  int low = tmin < 0 || tmin >= g[0]->GetN() ? 0 : tmin; ///< t0 setting. in case user want to test only part of wf, pls set the tmin and tmax
  int high[ntraces]; ///< t last setting
  int Nuse[ntraces]; ///< number of wf bins
  int NuseMax = 0; ///< number of final wf bins 
  for (int ti = 0; ti < ntraces; ti++) {
    //! set t last and number of bins
    high[ti] = tmax <= 0 || tmax > g[ti]->GetN() ? g[ti]->GetN() : tmax; ///< do we have a setting for specific time range?
    Nuse[ti] = high[ti] - low;
    if (Nuse[ti] > NuseMax) NuseMax = Nuse[ti];
  }
  int fin_len = pad_len ? pad_len : NuseMax; ///< let set practical wf length we will use

  //! f-domain configuration
  int spectrum_N = pad_len ?  pad_len / 2 + 1 : NuseMax / 2 + 1; ///< rfft length
  double df = pad_len ? 1./(pad_len * dt) : 1./(NuseMax * dt);
  r.phases.insert(r.phases.end(),ntraces, std::vector<double>());
  r.phases_errs.insert(r.phases_errs.end(),ntraces, std::vector<double>());
  r.amps.insert(r.amps.end(),ntraces, std::vector<double>());
  r.amps_errs.insert(r.amps_errs.end(),ntraces, std::vector<double>());

  //! wf padding
  TGraph * gPadded[ntraces]; ///< preparing empty pad for storing subtraction results
  memset(gPadded,0,sizeof(gPadded));
  for (int ti = 0; ti < ntraces; ti++) {
    if (Nuse[ti] < NuseMax || pad_len > NuseMax) {
      gPadded[ti] = new TGraph(pad_len);
      memcpy(gPadded[ti]->GetX(), g[ti]->GetX(), Nuse[ti] *sizeof(double));
      memcpy(gPadded[ti]->GetY(), g[ti]->GetY(), Nuse[ti] *sizeof(double));
      for (int i = Nuse[ti]; i < gPadded[ti]->GetN(); i++) {
        gPadded[ti]->GetX()[i] = gPadded[ti]->GetX()[i-1] + dt;
        gPadded[ti]->GetY()[i] = 0; //should be unnecessary
      }
    }
  }

  //! choose padded or not padded wf
  TGraph * ig[ntraces];
  for (int ti = 0; ti  < ntraces; ti++) {
    ig[ti] = gPadded[ti] ? gPadded[ti] : g[ti]; ///< choose padded or not padded wf
  }

  //! counting tries and failes
  int bads = bad_len; ///< numper of bad frequencies
  std::vector<int> nfails(spectrum_N); ///< number of times when subtraction is not lowring amplitude of fft below threshold in each frequency
  std::vector<int> nskips(spectrum_N); ///< specific frequency that we tried enough. if element is bigger than 0, we are skipping (giving up) that frequency 
  int max_i = 0; ///< bad frequency from random search
  while(bads != 0 or max_i != -1) 
  {
    //! pick guess index from input bad frequency index
    int bad_idx;
    int bad_idx_idx = -1;
    for (int b = 0; b < bad_len; b++) {
      if (bad_idxs[b] > -1) {
        bad_idx = bad_idxs[b];
        bad_idx_idx = b;
        break;
      }
    }

    //! making fft and pick up guess
    double guess_ph[ntraces];
    double guess_A = 0;
    double guess_f;
    for (int ti = 0; ti  < ntraces; ti++)
    {
      FFTWComplex * the_fft = FFTtools::doFFT(fin_len, ig[ti]->GetY() + low);

      //! Do random search after testing input bad frequencies
      if (bad_idx_idx == -1) {
        TGraph * dB_spectra = new TGraph(spectrum_N);
        for (int i = 0; i < spectrum_N; i++) {
          dB_spectra->GetX()[i] = df * i;
          dB_spectra->GetY()[i] = 10 * log10(sqrt(the_fft[i].getAbsSq() / NuseMax * dt)); ///< save in dB scale
        }
        max_i = findMaxFreq(spectrum_N, dB_spectra->GetX(), dB_spectra->GetY(), &nfails[0], &nskips[0], thres, ant); ///< peak finding function
        if (max_i == -1) break; ///< there is nothing to look around
        else {
          bad_idx = max_i;
          //std::cout<<"urgh... nothing perfect..."<<std::endl;      
        }
        delete dB_spectra;
      }

      guess_f = bad_idx * df;
      guess_ph[ti] = the_fft[bad_idx].getPhase();
      guess_A += sqrt(the_fft[bad_idx].getAbsSq() / NuseMax);
      if (bad_idx > 0 && bad_idx < spectrum_N - 1) guess_A *= sqrt(2);
      delete [] the_fft;
    }

    if (max_i == -1) break; ///< there is nothing to look around. stop the while loop  
 
    //! minimization
    const double * x[ntraces]; ///< take out all value from Tgraph to array
    const double * y[ntraces];
    TGraph * xy_crop = new TGraph(NuseMax); ///< trim out zero pad for minimization
    for (int ti = 0; ti < ntraces; ti++) {
      for (int i = 0; i < NuseMax; i++) {
        xy_crop->GetX()[i] = ig[ti]->GetX()[i]+low;
        xy_crop->GetY()[i] = ig[ti]->GetY()[i]+low;
      }
      x[ti] = xy_crop->GetX();
      y[ti] = xy_crop->GetY();
      //x[ti] = ig[ti]->GetX();
      //y[ti] = ig[ti]->GetY();
    }
    fitter.setGuess(guess_f, ntraces, guess_ph, guess_A); ///< lots of things are happening inside of this function...
    fitter.doFit(ntraces, Nuse, x,y,w, 0); ///< also here too...
    delete xy_crop;

    //! subtraction
    //! bring down the fit results into earth
    double freq_temp = fitter.getFreq();
    double amp_temp[ntraces];
    double phase_temp[ntraces];
    double sub_A = 0;
    r.freqs.push_back(freq_temp); ///< store fit parameters
    r.freqs_errs.push_back(fitter.getFreqErr());
    for (int ti = 0; ti < ntraces; ti++)
    {
      amp_temp[ti] = fitter.getAmp()[ti];
      phase_temp[ti] = fitter.getPhase()[ti];
      r.phases[ti].push_back(phase_temp[ti]);
      r.phases_errs[ti].push_back(fitter.getPhaseErr()[ti]);
      r.amps[ti].push_back(amp_temp[ti]);
      r.amps_errs[ti].push_back( fitter.getAmpErr()[ti]);
      double A = w ? w[ti] * amp_temp[ti] : amp_temp[ti];

      //! subtraction the fit results as a sine wave
      for (int i = 0; i < fin_len; i++) {
        if (i < NuseMax) ig[ti]->GetY()[i] -= A * sin(2*TMath::Pi() * freq_temp *ig[ti]->GetX()[i] + phase_temp[ti]); ///< apply subtraction only in original wf
        else ig[ti]->GetY()[i] = 0;
        //ig[ti]->GetY()[i] -= A * sin(2*TMath::Pi() * freq_temp *ig[ti]->GetX()[i] + phase_temp[ti]);
      }
    
      //! fft after subtraction
      FFTWComplex * the_fft = FFTtools::doFFT(fin_len, ig[ti]->GetY() + low);
      sub_A += sqrt(the_fft[bad_idx].getAbsSq() / NuseMax);
      if (bad_idx > 0 && bad_idx < spectrum_N - 1) sub_A *= sqrt(2);
      delete [] the_fft;
    }

    //! threhold check
    double ratio = 10 * log10(sub_A / sqrt(2) * sqrt(dt)) - baseline_fft[ant][bad_idx]; ///< do it in dB scale
    /*
    std::cout<<"!!!!!!!!!!!!!!!!freq_idx: "<<bad_idx<<std::endl;
    std::cout<<"freqs: "<<guess_f<<std::endl;

    std::cout<<"!!freq_fit: "<<freq_temp<<std::endl;
    std::cout<<"!!amp_fit: "<<amp_temp[0]<<std::endl;
    std::cout<<"!!amp_guess: "<<guess_A<<std::endl;
    std::cout<<"!!amp_sub: "<<sub_A<<std::endl;
    std::cout<<"!!amp_base: "<<baseline_fft[ant][bad_idx]*sqrt(2)<<std::endl;
    std::cout<<"!!amp_ratio: "<<sub_A/guess_A<<std::endl;

    std::cout<<"guess_A: "<<10 * log10(guess_A / sqrt(2) * sqrt(dt))<<std::endl; 
    std::cout<<"sub_A: "<<10 * log10(sub_A / sqrt(2) * sqrt(dt))<<std::endl; 
    std::cout<<"base: "<<10 * log10(baseline_fft[ant][bad_idx])<<std::endl; 
    std::cout<<"guess_A - base: "<<10 * (log10(guess_A / sqrt(2) * sqrt(dt)) - log10(baseline_fft[ant][bad_idx]))<<std::endl;
    std::cout<<"sub_A - base: "<<ratio<<std::endl;
    */
    if (ratio < thres) {
      //! if it is subtracted well, exclude bad freqeuncy that we just used
      if (bad_idx_idx != -1) {
        bad_idxs[bad_idx_idx] = -1; 
        bads -= 1;
      } else {
        nskips[bad_idx]++;
      }
    } else {
      //! if it didnt work well, give bad freqeuncy to another chance. But within penalty 
      nfails[bad_idx]++;
      if (nfails[bad_idx] > 2) { ///< had enough. Lets give up
        if (bad_idx_idx != -1) {
          bad_idxs[bad_idx_idx] = -1;
          bads -= 1;
        } else {
          nskips[bad_idx]++;
          //std::cout<<"dont try me again!"<<std::endl;
        }
      }
    } 
  }///< end of while loop

  //! cleanup the mass
  for (int ti = 0; ti < ntraces; ti++) { 
    for (int i = 0; i < g[ti]->GetN(); i++) {
      g[ti]->GetY()[i] = ig[ti]->GetY()[i]; ///< time to replace otiginal wf
    }
    if (gPadded[ti]) delete gPadded[ti]; 
  }
#ifdef SINE_SUBTRACT_PROFILE
  printf("Time for SineSubtract::subtractCW(): "); 
  sw.Print("u"); 
  printf("nattempts: %d\n",nattempts); 

#endif
}


int FFTtools::SineSubtract::getNSines() const 
{
  return r.freqs.size(); 
}

void FFTtools::SineSubtract::makeSlides(const char * title, const char * fileprefix, const char * outdir , const char * format, bool standalone) const 
{

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  Int_t orig_msz = gStyle->GetMarkerSize();
  Int_t orig_mst = gStyle->GetMarkerStyle();
  Int_t orig_lt  = gStyle->GetLineWidth();

  gStyle->SetMarkerSize(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineWidth(orig_lt*4);

  int orig_stat = gStyle->GetOptStat(); 
  gStyle->SetOptStat(0); 

 
  TCanvas canvas("slidecanvas","SlideCanvas",4000,3000); 
  canvas.Divide(3,1); 

  TGraph g; 
  int niter = spectra.size(); 
  TH2I poweraxis("poweraxis","Iteration 0", 10, 0,niter, 10,0, r.powers[0]*1.1); 
  poweraxis.GetXaxis()->SetTitle("iteration"); 
  poweraxis.GetYaxis()->SetTitle("power"); 
  if (gs.size() > 1) 
  {
    canvas.cd(2)->Divide(1,gs.size()); 
  }

  FILE * texfile = fopen(TString::Format("%s/%s.tex", outdir, fileprefix),"w"); 

  if(standalone)
  {
    fprintf(texfile,"\\documentclass[hyperref={pdfpagelabels=false}]{beamer} \\mode<presentation> { \\usetheme{Boadilla} \\usecolortheme{beaver} }\n"); 
    fprintf(texfile, "\\setbeamertemplate{navigation symbols}{}\n"); 
    fprintf(texfile,"\\begin{document}\n"); 
  }

  for (int i = 0; i < niter; i++) 
  {

    canvas.cd(1); 
    poweraxis.SetTitle(i < niter -1 ? TString::Format("Iteration %d, Sub Freq = %g Ghz",i, r.freqs[i]) : TString::Format("Iteration %d (final)", i)); 

    poweraxis.Draw(); 
    g.SetPoint(i, i, r.powers[i]); 
    g.Draw("lp"); 
    canvas.cd(2); 

    int old_gs_width[gs.size()];

    for (size_t j = 0; j < gs.size(); j++) 
    {
      if (gs.size() > 1) 
      {
        canvas.cd(2)->cd(j+1); 
      }
      old_gs_width[j]= gs[j][i]->GetLineWidth(); 
      gs[j][i]->SetLineWidth(3); 
      gs[j][i]->Draw("al"); 

      if ((int) env_gs[j].size() > i) 
      {
        env_gs[j][i]->Draw("lsame"); 
      }
    }


    canvas.cd(3)->SetLogy(); 

    int old_spectrum_width = spectra[i]->GetLineWidth(); 
    spectra[i]->SetLineWidth(5); 
    spectra[i]->Draw("alp"); 


    TString canvfile = TString::Format("%s/%s_%d.%s", outdir, fileprefix, i, format); 
    canvas.SaveAs(canvfile); 

    spectra[i]->SetLineWidth(old_spectrum_width); 
    for (size_t j = 0; j < gs.size(); j++) 
    {
      gs[j][i]->SetLineWidth(old_gs_width[j]); 
    }

    fprintf(texfile, "\\begin{frame}\n"); 
    fprintf(texfile, "\t\\frametitle{%s (iteration %d)}\n", title, i); 
    fprintf(texfile, "\t\\begin{center}\n"); 
    fprintf(texfile, "\t\t\\includegraphics[width=4.2in]{%s_%d}\n",fileprefix,i); 
    fprintf(texfile, "\t\\end{center}\n"); 
    fprintf(texfile, "\\end{frame}\n\n"); 
    fflush(texfile); 
  }

  if (standalone)
  {
    fprintf(texfile,"\\end{document}\n"); 
  }
  fclose(texfile); 

  gROOT->ForceStyle(kFALSE);
  gROOT->SetBatch(kFALSE);

  gStyle->SetMarkerSize(orig_msz);
  gStyle->SetMarkerStyle(orig_mst);
  gStyle->SetLineWidth(orig_lt);
  gStyle->SetOptStat(orig_stat); 
}

void FFTtools::SineSubtract::makePlots(TCanvas * cpower, TCanvas * cw, int ncols) const 
{

  int nplots = spectra.size() * (1 + gs.size()); 
  if (!cpower) cpower = new TCanvas("sinesubtractpowerplot","SineSubtract Power Evolution",1000,600); 
  cpower->Clear(); 

  if (!cw) cw = new TCanvas("sinesubtractplots","SineSubtract Waveforms/Spectra",1000,600); 
  cw->Clear(); 
  cw->Divide(ncols, (nplots-1)/ncols +1); 


  std::vector<double> iterations(r.powers.size()); 
  for (unsigned i = 0; i < r.powers.size(); i++) iterations[i] = i;
  TGraph *gpower = new TGraph(r.powers.size(), &iterations[0], &r.powers[0]); 
  cpower->cd(); 
  gpower->SetTitle("Power vs. iteration"); 
  gpower->Draw("alp"); 


  for (size_t i = 0; i < spectra.size(); i++) 
  {
    cw->cd((1+gs.size())*i+1)->SetLogy(); 
    ((TGraph*)spectra[i])->Draw("alp"); 
    for (size_t j = 0; j < gs.size(); j++)
    {
      cw->cd((1+gs.size())*i+2+j); 
      ((TGraph*)gs[j][i])->Draw("alp"); 
      if (env_gs[j].size() > i) 
        ((TGraph*)env_gs[j][i])->Draw("lsame"); 
    }
  }

}


FFTtools::SineSubtract::~SineSubtract()
{
  delete [] power_estimator_params; 
  delete [] peak_option_params; 
  reset(); 
}




bool FFTtools::SineSubtract::allowedFreq(double freq, double df) const
{
  if (!fmin.size()) return true; 
 
  for (unsigned j = 0; j < fmin.size(); j++) 
  {
     if (freq+df >= fmin[j] && freq-df <= fmax[j]) 
     {
       return true; 
     }
  }

  return false; 
}

static __thread FFTtools::SavitzkyGolayFilter* savgol = 0; 
static __thread TSpectrum* tspect = 0; 



int FFTtools::SineSubtract::findMaxFreq(int Nfreq, const double * freq, const double * mag, const int * nfails, const int * nskips, double thres, int ant) const
{

  int max_i =-1; 
  double max = 0; 
  double df = freq[1]-freq[0]; 

  if (peak_option == NEIGHBORFACTOR) { // we can just iterate through and do this
    for (int i = 0; i < Nfreq;i++) {
      if (nskips[i]) continue;
      if (!allowedFreq(freq[i],df)) continue;

      double val = mag[i] - baseline_fft[ant][i];
      if (val < thres) continue;

      double adjval = val;
      if (nfails[i]) adjval /= pow(1 + nfails[i],nfail_exponent);       

      if (adjval > max) {
        max = adjval;
        max_i = i;
      }
    }
  }
  /*
  if (peak_option == GLOBALMAX || peak_option == NEIGHBORFACTOR)   // we can just iterate through and do this
  {
    for (int i = 0; i < Nfreq;i++)
    {

      //TODO this is a very inefficient way to do this 
      if (!allowedFreq(freq[i],df)) continue; 

      double val =mag[i]; 
      //check neighbors if we are doing that
      if (peak_option == NEIGHBORFACTOR) 
      {

        double neigh_low  = i == 0 ? DBL_MAX : mag[i-1]; 
        double neigh_high  = i == Nfreq-1 ? DBL_MAX : mag[i+1]; 
        if (val * peak_option_params[NF_NEIGHBOR_FACTOR] <= TMath::Min(neigh_low, neigh_high)) 
        {
          continue;
        }
      }

      double adjval = val;
      // adjust the value
      if (nfails[i]) adjval /= pow(1 + nfails[i],nfail_exponent); 


      if (adjval > max) 
      {
        max = adjval; 
        max_i = i; 
      }
    }
  }
  */
  else if (peak_option == TSPECTRUM)
  {
#if ROOT_VERSION_CODE <= ROOT_VERSION(6,0,0)
    fprintf(stderr,"TSPECTRUM option not supported for ROOT < 6 due to API change\n"); 
#else
    if (!tspect) tspect = new TSpectrum(32); //32 peaks should be good enough for anyone! 

    double out[Nfreq]; 

    int npeaks = tspect->SearchHighRes((double*)mag,out,Nfreq,
        peak_option_params[TS_SIGMA], peak_option_params[TS_THRESHOLD], true,
        (int) peak_option_params[TS_NDECONV_ITERATIONS], true, (int)
        peak_option_params[TS_AVERAGE_WINDOW]) ; 


    for (int p = 0; p < npeaks; p++) 
    {
       int i= tspect->GetPositionX()[p] +0.5;  // I guess round to nearest bin? 
       if (!allowedFreq(freq[i],df)) continue; 

        double val =mag[i]; 
        if (nfails[i]) val /= pow(1 + nfails[i],nfail_exponent); 

        if (val > max) 
        {
          max = val; 
          max_i = i; 
        }
    }

    // it's possible that none of our peaks are valid, in which case we will just pick the maximum from the
    // deconvolved array :( 

    if (max_i < 0) 
    {
      for (int i = 0; i < Nfreq; i++) 
      {
        if (!allowedFreq(freq[i],df)) continue; 
        double val =mag[i]; 
        if (nfails[i]) val /= pow(1 + nfails[i],nfail_exponent); 
        if (val > max) 
        {
          max = val; 
          max_i = i; 
        }
      }
    }
#endif
  }

  else if (peak_option == SAVGOLSUB)
  {

    //check if we need a new savgol 
    if (!savgol || 
        savgol->getOrder() != (int) peak_option_params[SGS_ORDER] || 
        savgol->getWidthLeft() != (int) peak_option_params[SGS_WIDTH] )
    {
      if (savgol) delete savgol; 
      savgol = new FFTtools::SavitzkyGolayFilter(peak_option_params[SGS_ORDER], peak_option_params[SGS_WIDTH]); 
    }
    
     
    double bg[Nfreq]; 
    savgol->filterOut(Nfreq, mag, bg); 
    
    for (int i = 0; i < Nfreq; i++)
    {
      //TODO this is a very inefficient way of doing this 
      if (!allowedFreq(freq[i],df)) continue; 
      double val = mag[i] - bg[i]; 

      if (nfails[i]) val /= pow(1 + nfails[i], nfail_exponent); 
      if (val > max) 
      {
        max = val; 
        max_i = i; 
      }
    }
  }


  return max_i; 
}


static const double default_FFT_params[] = {0.} ; 
static const double default_LOMBSCARGLE_params[] = {2.} ; 
static double dummy = 0; 

void FFTtools::SineSubtract::setPowerSpectrumEstimator(PowerSpectrumEstimator estimator, const double * p) 
{
  power_estimator = estimator; 
  int N = estimator == FFT? FFT_NPARAMS : 
          estimator == LOMBSCARGLE? LS_NPARAMS :
          0; 
  
  if (power_estimator_params) delete [] power_estimator_params; 
  power_estimator_params = new double[N]; 

  if (N) 
  {
    memcpy(power_estimator_params,
           p? p : 
           estimator == FFT?  default_FFT_params  : 
           estimator == LOMBSCARGLE?  default_LOMBSCARGLE_params  : 
           &dummy,//to avoid compiler warning 
           sizeof(double) *N); 
  }

}

static const double default_NEIGHBORFACTOR_params[] = {0.15}; 
static const double default_TSPECTRUM_params[] = {2,0.05,3,3}; 
static const double default_SAVGOLSUB_params[] = {3,10}; 



void FFTtools::SineSubtract::setPeakFindingOption(PeakFindingOption option, const double * p)
{
  peak_option = option; 

  int N = option == GLOBALMAX? 0 : 
          option == NEIGHBORFACTOR? NF_NPARAMS :         
          option == TSPECTRUM? TS_NPARAMS : 
          option == SAVGOLSUB ? SGS_NPARAMS : 
          0; 

  if (peak_option_params) delete[] peak_option_params; 

  peak_option_params = new double[N]; 

  if (N) 
  {
    memcpy(peak_option_params, 
           p? p : 
           option == NEIGHBORFACTOR?  default_NEIGHBORFACTOR_params : 
           option == TSPECTRUM?  default_TSPECTRUM_params : 
           option == SAVGOLSUB?  default_SAVGOLSUB_params : 
           &dummy ,  //to avoid compiler warning
           sizeof(double) *N); 
  }

}

static const double default_ENV_HILBERT_params[] = {3}; 
static const double default_ENV_RMS_params[] = {3,5}; 
static const double default_ENV_PEAK_params[] = {3,5}; 

void FFTtools::SineSubtract::setEnvelopeOption(EnvelopeOption option, const double * p) 
{
  envelope_option = option; 

  int N = option == ENV_NONE? 0 : 
          option == ENV_HILBERT? ENV_HILBERT_NPARAMS :         
          option == ENV_RMS ? ENV_RMS_NPARAMS : 
          option == ENV_PEAK ? ENV_PEAK_NPARAMS : 
          0; 

  if (envelope_option_params) delete[] envelope_option_params; 

  envelope_option_params = new double[N]; 

  if (N) 
  {
    memcpy(envelope_option_params, 
           p? p : 
           option == ENV_HILBERT?  default_ENV_HILBERT_params : 
           option == ENV_RMS?  default_ENV_RMS_params : 
           option == ENV_PEAK?  default_ENV_PEAK_params : 
           &dummy ,  //to avoid compiler warning
           sizeof(double) *N); 
  }
}




ClassImp(FFTtools::SineSubtractResult); 

#endif

