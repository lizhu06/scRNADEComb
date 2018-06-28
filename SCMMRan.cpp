#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_FACTOR(group);
  DATA_VECTOR(y);
  DATA_MATRIX(Xu); // designed matrix for mu
  DATA_MATRIX(Xp); // designed matrix for p
  DATA_VECTOR(osu); // offset for mu
  DATA_VECTOR(osp); // offset for p
  PARAMETER(logtheta); // the inverse of dispersion
  PARAMETER(logsd_mu); // sd of RE
  PARAMETER(logsd_p); // sd of RE
  PARAMETER_VECTOR(b_mu); // RE
  PARAMETER_VECTOR(b_p); // RE
  PARAMETER_VECTOR(BETAu);
  PARAMETER_VECTOR(BETAp);

  Type theta = exp(logtheta);
  Type sd_mu = exp(logsd_mu);
  Type sd_p = exp(logsd_p);

  ADREPORT(theta);
  ADREPORT(sd_mu);
  ADREPORT(sd_p);

  int nobs=y.size();
  int ngroups=b_mu.size();
  
  vector<Type> etau = Xu*BETAu;
  vector<Type> etap = Xp*BETAp;

  Type nll = 0;
  for(int j=0;j<ngroups;j++){
    nll-=(dnorm(b_mu[j],Type(0),sd_mu,true));
  }

  for(int j=0;j<ngroups;j++){
    nll-=(dnorm(b_p[j],Type(0),sd_p,true));
  }


  for(int i=0;i<nobs;i++){
    int j=group[i];
    Type p1 = exp(osp[i]+etap[i]+b_p[j]);
    Type ratio = p1/(Type(1)+p1);
    Type meanNB = exp(osu[i]+etau[i]+b_mu[j]);
    //Type meanNB = exp(osu[i]+etau[i]);
    Type varNB = meanNB + meanNB*meanNB/theta;
    int ind;
    if(y[i]==0){
      ind = 1;
    }else{
      ind = 0;
    }
    //nll-=log(ratio*ind+(1-ratio)*dnbinom2(y[i],meanNB,varNB,0)); model zero part
    nll-=log((1-ratio)*ind+ratio*dnbinom2(y[i],meanNB,varNB,0)); //model NB part
  }
 return nll;
}

