#include <iostream>
#include <openturns/OT.hxx>
#include "otsubsetinverse/SubsetInverseSampling.hxx"

using namespace OT;
using namespace OTSubsetInverse;

int main(int argc, char **argv)
{
  RandomGenerator::SetSeed(0);

  Distribution myDistribution = Normal();
  RandomVector vect(myDistribution);
  
  SymbolicFunction limitState("X","2*X");
  RandomVector output(limitState, vect);

  Event myEvent(output, ComparisonOperator(Less()), 0.0);
    
  Scalar finalProbability(0.0001);

  OTSubsetInverse::SubsetInverseSampling ssi(myEvent, finalProbability);
  std::cout << ssi << std::endl;
  ssi.run();
  std::cout << ssi.getResult() << std::endl;
  return 0;
}
