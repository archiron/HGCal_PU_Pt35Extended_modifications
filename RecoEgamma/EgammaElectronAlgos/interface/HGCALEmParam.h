#ifndef HGCALEmParam_h
#define HGCALEmParam_h


class HGCALEmParam
{

  public:

  HGCALEmParam(double radiationLength, double moliereRadius, double criticalEnergy):
  radiationLength_(radiationLength), moliereRadius_(moliereRadius), criticalEnergy_(criticalEnergy) 
  {
    meant0_=-1.396;
    meant1_=1.007;
    meanalpha0_=-0.0433;
    meanalpha1_=0.540;
    sigmalnt0_=-2.506;
    sigmalnt1_=1.245;
    sigmalnalpha0_=-0.08442;
    sigmalnalpha1_=0.7904;
    corrlnalphalnt0_=0.7858;
    corrlnalphalnt1_=-0.0232;    
  };
  
  const double meanT(double lny) {return meant0_+meant1_*lny;}
  const double meanAlpha(double lny) {return meanalpha0_+meanalpha1_*lny;}
  const double sigmaLnT(double lny) {return 1./(sigmalnt0_+sigmalnt1_*lny);}
  const double sigmaLnAlpha(double lny) {return 1./(sigmalnalpha0_+sigmalnalpha1_*lny);}
  const double correlationAlphaT(double lny) {return corrlnalphalnt0_+corrlnalphalnt1_*lny;}

  const double getRadiationLength() {return radiationLength_;}
  const double getMoliereRadius() {return moliereRadius_;}
  const double getCriticalEnergy() {return criticalEnergy_;}

private:
  
  // hgcal parameters
  double radiationLength_;
  double moliereRadius_;
  double criticalEnergy_;

  // longitudinal parametrisation
  double meant0_;
  double meant1_;
  double meanalpha0_;
  double meanalpha1_;
  double sigmalnt0_;
  double sigmalnt1_;
  double sigmalnalpha0_;
  double sigmalnalpha1_;
  double corrlnalphalnt0_;
  double corrlnalphalnt1_;

};
#endif
