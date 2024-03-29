 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * This code was autogenerated by RooClassFactory                            * 
  *****************************************************************************/ 

 // Your description goes here... 

#include "LIP/Top/interface/HeavyFlavorDiffPDF.h" 

 ClassImp(HeavyFlavorDiffPDF) 

 HeavyFlavorDiffPDF::HeavyFlavorDiffPDF(const char *name, const char *title, 
				RooAbsReal& _bmult,
				RooAbsReal& _r,
				RooAbsReal& _eb1, RooAbsReal& _eb2,
				RooAbsReal& _eq1, RooAbsReal& _eq2,
				RooAbsReal& _fcorrect,
				RooAbsReal& _fttbar,
				RooAbsReal& _fsingletop,
				RooAbsReal& _catocc):
   RooAbsPdf(name,title), 
   bmult("bmult","bmult",this,_bmult),
   r("r","r",this,_r),
   eb1("eb1","eb1",this,_eb1),
   eb2("eb2","eb2",this,_eb2),
   eq1("eq1","eq1",this,_eq1),
   eq2("eq2","eq2",this,_eq2),
   fcorrect("fcorrect","fcorrect",this,_fcorrect),
   fttbar("fttbar","fttbar",this,_fttbar),
   fsingletop("fsingletop","fsingletop",this,_fsingletop),
   catocc("catocc","catocc",this,_catocc)
 { 
   //default value
   setCategory(1);
 }


HeavyFlavorDiffPDF::HeavyFlavorDiffPDF(const HeavyFlavorDiffPDF& other, const char* name) :  
   RooAbsPdf(other,name), 
   bmult("bmult",this,other.bmult),
   r("r",this,other.r),
   eb1("eb1",this,other.eb1),
   eb2("eb2",this,other.eb2),
   eq1("eq1",this,other.eq1),
   eq2("eq2",this,other.eq2),
   fcorrect("fcorrect",this,other.fcorrect),
   fttbar("fttbar",this,other.fttbar),
   fsingletop("fsingletop",this,other.fsingletop),
   catocc("catocc",this,other.catocc),
   icat(other.icat)
 { 

 } 



 Double_t HeavyFlavorDiffPDF::evaluate() const 
 { 
   Int_t nBtags=int(bmult);
   if(nBtags<0) return 0;

   //evaluate the probability functions
   double prob=_evaluate(nBtags);
   return prob;
 } 


// 
Double_t HeavyFlavorDiffPDF::_evaluate(int nBtags) const
{
  Double_t prob(0);

  //update the values for the alphas
  int nevents = catocc;
  int npairs=2*2*nevents;
	
  if((2*fttbar+fsingletop)*nevents==0) return 0;
  double alpha=(fcorrect*npairs)/((2*fttbar+fsingletop)*nevents);
  if(alpha<0 || alpha>1) { cout<< "[HeavyFlavorDiffPDF::_evaluate] unphysical value for alpha: " << alpha << endl; return 0; }
  double alpha2=pow(alpha,2)*fttbar;
  double alpha1=2*alpha*(1-alpha)*fttbar+alpha*fsingletop;
  double alpha0=1-alpha2-alpha1;

  //compute probability
  prob =  alpha2*_evaluateKernel(0,nBtags);
  prob += alpha1*_evaluateKernel(1,nBtags);
  prob += alpha0*_evaluateKernel(2,nBtags);
  
  //the probability
  return prob;
}
    

//
Double_t HeavyFlavorDiffPDF::_evaluateKernel(int nMisassignments, int nBtags, bool doNorm) const
{
  Double_t prob(0);
  
  //check ranges
  if(nMisassignments<0 || nMisassignments>2) return prob;
  if(nBtags<0 || nBtags>2) return prob;
  
  switch(nMisassignments)
    {
      //pure signal
    case 0:
      {
	switch(nBtags)
	  {
	  case 0:
	    {
	      prob =  TMath::Power(r,2)*(1-eb1)*(1-eb2)*_normalizedAcceptance(1);
	      prob += 2*r*(1-r)*((1-eb2)*(1-eq1)+(1-eb1)*(1-eq2))*_normalizedAcceptance(0.5);
	      prob += TMath::Power((1-r),2)*(1-eq1)*(1-eq2)*_normalizedAcceptance(0);
	      break;
	    }
	  case 1:
	    {
	      prob =  TMath::Power(r,2)*((1-eb1)*eb2+(1-eb2)*eb1)*_normalizedAcceptance(1);
	      prob += r*(1-r)*((1-eb1)*eq2+eb1*(1-eq2)+(1-eb1)*eq2+eb1*(1-eq2))*_normalizedAcceptance(0.5);
	      prob += TMath::Power((1-r),2)*(eq1*(1-eq2)+(1-eq1)*eq2)*_normalizedAcceptance(0);
	      break;
	    }
	  case 2:
	    {
	      prob =  TMath::Power(r,2)*eb1*eb2*_normalizedAcceptance(1);
	      prob += r*(1-r)*(eb1*eq2+eq1*eb2)*_normalizedAcceptance(0.5);
	      prob += TMath::Power((1-r),2)*eq1*eq2*_normalizedAcceptance(0);
	      break;
	    }
	  }
	break;
      }
      
      //signal with 1 misassignment
    case 1:
      {
	switch(nBtags)
	  {
	  case 0:
	    {
	      prob =  TMath::Power(r,2)*((1-eb1)*(1-eq2)+(1-eq1)*(1-eb2))*_normalizedAcceptance(1);
	      prob += r*(1-r)*( ((1-eb1)+(1-eq1))*(1-eq2) + (1-eq1)*((1-eb2)+(1-eq2)) )  *_normalizedAcceptance(0.5);
	      prob += TMath::Power((1-r),2)*(1-eq1)*(1-eq2)*_normalizedAcceptance(0);
	      break;
	    }
	  case 1:
	    {
	      prob =  TMath::Power(r,2)*( eb1*(1-eq2) + eq1*(1-eb2) ) * _normalizedAcceptance(1);
	      prob += r*(1-r)*( (eb1+eq1)*(1-eq2) + ((1-eb1)+(1-eq1))*eq2 +
				(1-eq1)*(eb2+eq2) + eq1*((1-eb2)+(1-eq2)) ) * _normalizedAcceptance(0.5);
	      prob += TMath::Power((1-r),2)*( eq1*(1-eq2) + (1-eq1)*eq2 )*_normalizedAcceptance(0);
	      break;
	    }
	  case 2:
	    {
	      prob =  TMath::Power(r,2)*(eb1*eq2 + eq1*eq2)*_normalizedAcceptance(1);
	      prob += r*(1-r)*((eb1+eq1)*eq2 + eq1*(eb2+eq2))*_normalizedAcceptance(0.5);
	      prob += TMath::Power((1-r),2)*eq1*eq2*_normalizedAcceptance(0);
	      break;
	    }
	  }
	break;
      }
      
      //pure background + totally misreconstructed signal
    case 2:
      {
	switch(nBtags)
	  {
	  case 0:
	    {
	      prob = (1-eq1)*(1-eq2);
	      break;
	    }
	  case 1:
	    {
	      prob = eq1*(1-eq2)+(1-eq1)*eq2;
	      break;
	    }
	  case 2:
	    {
	      prob = eq1*eq2;
	      break;
	    }
	  }
	break;
      }
    }
  

  Double_t norm(1.0);
  if(doNorm)
    {
      norm = _evaluateKernel(nMisassignments,0,false)
	+ _evaluateKernel(nMisassignments,1,false)
	+ _evaluateKernel(nMisassignments,2,false);
    }	
  prob /= norm;	
  
  //the result
  return prob;
}

//
Double_t HeavyFlavorDiffPDF::_normalizedAcceptance(Double_t R)  const
{
  return 1.0;//R*(1.0-lfacceptance)+1.0;
}
