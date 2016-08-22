#include <iostream>
#include <openturns/OT.hxx>
#include "otsubsetinverse/SubsetSampling.hxx"

using namespace OT;
using namespace OTSubsetInverse;

int main(int argc, char **argv)
{
  Distribution myDistribution = Normal();
  RandomVector vect(myDistribution);
  
  NumericalMathFunction limitState("X","Y","2*X");
  RandomVector output(limitState, vect);

  Event myEvent(output, ComparisonOperator(Less()), 0.0);
  
  OTSubsetInverse::SubsetSampling a(myEvent);
  std::cout << a << std::endl; 
  return 0;
}
