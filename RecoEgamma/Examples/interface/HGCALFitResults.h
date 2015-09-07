#ifndef HGCALFitResults_h
#define HGCALFitResults_h

//===================================================================
// Purpose: HGCAL em shower fit results interface
// Author: Claude Charlot - LLR- Ecole Polytechnique
// 02/2015
//===================================================================

class HGCALFitResults
{

  public:

  HGCALFitResults(double chi2, int ndf, double norm, double alpha, double
   invbeta): chi2_(chi2), ndf_(ndf), normalization_(norm), alpha_(alpha),
   invbeta_(invbeta) 
  {
  }
  
  const double chi2() {return chi2_;}
  const double ndf() {return ndf_;}
  const double normalization() {return normalization_;}
  const double alpha() {return alpha_;}
  const double invbeta() {return invbeta_;}

private:
  
  // longitudinal fit results
  double chi2_;
  int ndf_;
  double normalization_;
  double alpha_;
  double invbeta_;

};
#endif
