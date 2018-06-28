#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_MATRIX(Xu); // designed matrix for mu
  DATA_MATRIX(Xp); // designed matrix for p
  DATA_VECTOR(osu); // offset for mu
  DATA_VECTOR(osp); // offset for p
  PARAMETER(logtheta); // the inverse of dispersion
  //PARAMETER(logsd); // sd of RE
  //PARAMETER_VECTOR(b); // RE
  PARAMETER_VECTOR(BETAu);
  PARAMETER_VECTOR(BETAp);

  Type theta = exp(logtheta);
  //Type sd = exp(logsd);
  
  ADREPORT(theta);
  //ADREPORT(sd);

  int nobs=y.size();
  //int ngroups=b.size();
  
  vector<Type> etau = Xu*BETAu;
  vector<Type> etap = Xp*BETAp;

  Type nll = 0;
  //for(int j=0;j<ngroups;j++){
  //  nll-=(dnorm(b[j],Type(0),sd,true));
  //}
  for(int i=0;i<nobs;i++){
    //int j=group[i];
    Type p1 = exp(osp[i]+etap[i]);
    Type ratio = p1/(Type(1)+p1);
    //Type meanNB = exp(osu[i]+etau[i]+b[j]);
    Type meanNB = exp(osu[i]+etau[i]);
    Type varNB = meanNB + meanNB*meanNB/theta; //theta is 1/phi
    int ind;
    if(y[i]==0){
      ind = 1;
    }else{
      ind = 0;
    }
    //nll-=log(ratio*ind+(1-ratio)*dnbinom2(y[i],meanNB,varNB,0)); model zero
    nll-=log((1-ratio)*ind+ratio*dnbinom2(y[i],meanNB,varNB,0)); //model non-zero part
  }
 return nll;
}

