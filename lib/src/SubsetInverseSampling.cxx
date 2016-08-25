//                                               -*- C++ -*-
/**
 *  @brief SubsetInverseSampling
 *
 *  Copyright 2005-2016 Airbus-EDF-IMACS-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "otsubsetinverse/SubsetInverseSampling.hxx"
#include "otsubsetinverse/SubsetInverseSamplingResult.hxx"

using namespace OT;

namespace OTSubsetInverse
{

CLASSNAMEINIT(SubsetInverseSampling);

static Factory<SubsetInverseSampling> Factory_SubsetInverseSampling;

const UnsignedInteger SubsetInverseSampling::DefaultMaximumOuterSampling = 10000;
const NumericalScalar SubsetInverseSampling::DefaultTargetProbability = 0.1;
const NumericalScalar SubsetInverseSampling::DefaultProposalRange = 2.0;
const NumericalScalar SubsetInverseSampling::DefaultBetaMin = 2.0;


/* Default constructor */
SubsetInverseSampling::SubsetInverseSampling()
: Simulation()
, proposalRange_(0.)
, targetProbability_(0.)
, iSubset_(false)
, betaMin_(0.)
, keepEventSample_(false)
, numberOfSteps_(0)
, finalTargetProbability_(0.)
{
}


/* Constructor with parameters */
SubsetInverseSampling::SubsetInverseSampling(const Event & event,
                               const NumericalScalar finalTargetProbability,
                               const NumericalScalar proposalRange,
                               const NumericalScalar targetProbability)
: Simulation(event)
, proposalRange_(proposalRange)
, targetProbability_(targetProbability)
, iSubset_(false)
, betaMin_(DefaultBetaMin)
, keepEventSample_(false)
, numberOfSteps_(0)
, finalTargetProbability_(finalTargetProbability)
{
  setMaximumOuterSampling( DefaultMaximumOuterSampling );// overide simulation default outersampling
  UnsignedInteger outputDimension = event.getFunction().getOutputDimension();
  if ( outputDimension > 1 )
    throw InvalidArgumentException(HERE) << "Output dimension for SubsetInverseSampling cannot be greater than 1, here output dimension=" << outputDimension;
}


/* Virtual constructor */
SubsetInverseSampling * SubsetInverseSampling::clone() const
{
  return new SubsetInverseSampling(*this);
}


/* Performs the actual computation. */
void SubsetInverseSampling::run()
{ 
  // First, initialize some parameters
  convergenceStrategy_.clear();
  numberOfSteps_ = 0;
  thresholdPerStep_.clear();
  gammaPerStep_.clear();
  coefficientOfVariationPerStep_.clear();
  probabilityEstimatePerStep_.clear();
  eventInputSample_.clear();
  eventOutputSample_.clear();
  Collection<NumericalSample> allLevelSample;

  dimension_ = getEvent().getAntecedent()->getDimension();

  if ( getMaximumCoefficientOfVariation() != ResourceMap::GetAsNumericalScalar( "Simulation-DefaultMaximumCoefficientOfVariation" ) )
    Log::Warn(OSS() << "The maximum coefficient of variation was set. It won't be used as termination criteria.");

  if ( targetProbability_ * getMaximumOuterSampling() < 1 )
    throw InvalidArgumentException(HERE) << "maximumOuterSampling (" << getMaximumOuterSampling() << ") should be >= " << ceil( 1.0 / targetProbability_ );
  
  if ( getMaximumOuterSampling() * getBlockSize() <= 100 )
    Log::Warn(OSS() << "The number of samples per step is very low : " << getMaximumOuterSampling()*getBlockSize() << ".");

  // perform isoprobabilistic transformation (the study is done in the standard space):
  standardEvent_ = StandardEvent(getEvent());

  NumericalScalar currentCoVsquare = 0.0;
  NumericalScalar varianceEstimate = 0.0;
  NumericalScalar coefficientOfVariationSquare = 0.0;

  // allocate input/output samples
  const UnsignedInteger maximumOuterSampling = getMaximumOuterSampling();
  const UnsignedInteger blockSize = getBlockSize();
  const UnsignedInteger N = maximumOuterSampling * blockSize;
  currentPointSample_ = NumericalSample( N, dimension_ );
  currentLevelSample_ = NumericalSample( N, getEvent().getFunction().getOutputDimension() );
  allLevelSample.add(currentLevelSample_);
  
  // Step 1: sampling
  for ( UnsignedInteger i = 0; i < maximumOuterSampling; ++ i )
  {
    NumericalSample inputSample;
    if (!iSubset_) {
      // crude MC
      inputSample = standardEvent_.getAntecedent()->getDistribution().getSample( blockSize );
    }
    else {
      // conditional sampling
      TruncatedDistribution truncatedChiSquare(ChiSquare(dimension_), betaMin_ * betaMin_, TruncatedDistribution::LOWER);
      Normal normal(dimension_);
      inputSample = NumericalSample(0, dimension_);
      for ( UnsignedInteger j = 0; j < blockSize; ++ j)
      {
        NumericalPoint direction = normal.getRealization();
        NumericalScalar norm = direction.norm();
        NumericalScalar radius = sqrt(truncatedChiSquare.getRealization()[0]);
        if (fabs(norm) > 1e-12)
        {
          radius *= 1.0 / norm;
        }
        inputSample.add(direction * radius);
      }
    }
    NumericalSample blockSample( standardEvent_.getFunction()( inputSample ) );
    for ( UnsignedInteger j = 0 ; j < blockSize; ++ j )
    {
      currentPointSample_[ i*blockSize+j ] = inputSample[j];
      currentLevelSample_[ i*blockSize+j ] = blockSample[j];
    }
  }
  ++ numberOfSteps_;

  // Stop if the wanted probability if greater than the target probability per step 
  Bool stop = finalTargetProbability_ >= targetProbability_;

  if  (stop)
  {
    setTargetProbability(finalTargetProbability_);
  }

  // computation of the first intermediate threshold with the sample create with a normal distribution */
  NumericalScalar currentThreshold = computeThreshold(); 
  
  // compute monte carlo probability estimate
  NumericalScalar probabilityEstimate = computeProbability(1.0, currentThreshold);
  
  if (iSubset_)
  {
    NumericalScalar correction = 1.0 - ChiSquare(standardEvent_.getImplementation()->getAntecedent()->getDistribution().getDimension()).computeCDF(betaMin_ * betaMin_);
    probabilityEstimate *= correction;
  }

  // if there is no subset step ...
  if ( stop )
  {
    if (probabilityEstimate > 0.0)
    {
      // ... compute coefficient of variation
      coefficientOfVariationSquare = (1.0 - probabilityEstimate) / (probabilityEstimate * currentLevelSample_.getSize() * 1.0);
      // ... compute variance estimate
      varianceEstimate = coefficientOfVariationSquare * probabilityEstimate * probabilityEstimate;
    }
  }
  
  thresholdPerStep_.add( currentThreshold );
  gammaPerStep_.add(0.);
  probabilityEstimatePerStep_.add(probabilityEstimate);
  coefficientOfVariationPerStep_.add(coefficientOfVariationSquare);
  
  // as long as the conditional failure domain do not overlap the global one
  while ( !stop )
  {    
    // prepare new seeds
    initializeSeed( currentThreshold );

    // create new points using MCMC
    generatePoints( currentThreshold );

    // save the level sample
    allLevelSample.add(currentLevelSample_);

    // compute new threshold
    currentThreshold = computeThreshold();

    // compute probability estimate on the current sample and group seeds at the beginning of the work sample
    NumericalScalar currentProbabilityEstimate = computeProbability( probabilityEstimate, currentThreshold );
    
    // update probability estimate
    probabilityEstimate *= currentProbabilityEstimate;

    // update stopping criterion
    stop = finalTargetProbability_ >= probabilityEstimate;

    if (stop)
    {
      // change the target probability of the final step
      setTargetProbability(finalTargetProbability_ / probabilityEstimatePerStep_[numberOfSteps_-1]);
      // compute the final threshold
      currentThreshold = computeThreshold();
      // compute the current probability estimate 
      NumericalScalar currentProbabilityEstimate = computeProbability( probabilityEstimatePerStep_[numberOfSteps_-1], currentThreshold );
      // update probability estimate
      probabilityEstimate = probabilityEstimatePerStep_[numberOfSteps_-1] * currentProbabilityEstimate;
    }

    // update coefficient of variation 
    NumericalScalar gamma = computeVarianceGamma( currentProbabilityEstimate, currentThreshold );
    currentCoVsquare = (1.0 - currentProbabilityEstimate) / (currentProbabilityEstimate * currentLevelSample_.getSize() * 1.0);
    coefficientOfVariationSquare += (1.0 + gamma) * currentCoVsquare;

    thresholdPerStep_.add( currentThreshold );
    gammaPerStep_.add(gamma);
    probabilityEstimatePerStep_.add(probabilityEstimate);
    coefficientOfVariationPerStep_.add(sqrt(coefficientOfVariationSquare));
    
    // stop if the number of subset steps is too high, else results are not numerically defined anymore
    if ( fabs( pow( probabilityEstimate, 2.) ) < SpecFunc::MinNumericalScalar )
      throw NotDefinedException(HERE) << "Probability estimate too low: " << probabilityEstimate;

    // compute variance estimate
    varianceEstimate = coefficientOfVariationSquare * pow( probabilityEstimate, 2. );

    ++ numberOfSteps_;
  }

  // compute the threshold distribution
  // sampling of the asymptotic pf distribution
  // Truncated distribution to avoid negative probability realizations when imprecise simulation
  Distribution probabilityDistribution = TruncatedDistribution(Normal(probabilityEstimate, sqrt(varianceEstimate)), 0, TruncatedDistribution::LOWER);
  NumericalScalar sizeSample = 10000;
  NumericalSample sampleProbDistribution = probabilityDistribution.getSample(sizeSample);
  sampleThreshold_ = NumericalSample(sizeSample, 1);
  // compute the corresponding threshold for each probability
  for ( UnsignedInteger i = 0; i < sizeSample; ++ i )
  { 
    currentLevelSample_ = allLevelSample[numberOfSteps_ - 1];
    NumericalScalar newTargetProbability(sampleProbDistribution[i][0] / probabilityEstimatePerStep_[numberOfSteps_-2]);
    // change the step when the probability is greater than the previous probability estimate step
    // use the previous step data
    NumericalScalar stepBackward(1);
    while (newTargetProbability >= 1)
    { 
      newTargetProbability = sampleProbDistribution[i][0] / probabilityEstimatePerStep_[numberOfSteps_ - 2 - stepBackward];
      currentLevelSample_ = allLevelSample[numberOfSteps_ - 1 - stepBackward];
      stepBackward ++;
    }
    setTargetProbability(newTargetProbability);
    NumericalPoint threshold(1, computeThreshold());
    sampleThreshold_[i] = threshold;
  }

  //update the event with the final threshold
  Event modified_event = Event(RandomVector(getEvent().getFunction(), getEvent().getAntecedent()), getEvent().getOperator(), currentThreshold);

  setResult( SubsetInverseSamplingResult(modified_event, probabilityEstimate, varianceEstimate, numberOfSteps_ * getMaximumOuterSampling(), getBlockSize(), sqrt( coefficientOfVariationSquare ), currentThreshold) );
  
  // keep the event sample if requested
  if (keepEventSample_)
  {
    eventInputSample_ = NumericalSample(0, dimension_);  
    eventOutputSample_ = NumericalSample (0, getEvent().getFunction().getOutputDimension());
    for ( UnsignedInteger i = 0; i < currentPointSample_.getSize(); ++ i )
    {
      if ( getEvent().getOperator()( currentLevelSample_[i][0], currentThreshold ) )
      {
        eventInputSample_.add( standardEvent_.getAntecedent()->getDistribution().getInverseIsoProbabilisticTransformation()(currentPointSample_[i]) );
        eventOutputSample_.add( currentLevelSample_[i] );
      }
    }
  }
  
  // free work samples
  currentLevelSample_.clear();
  currentPointSample_.clear();
}


/* Compute the block sample */
NumericalSample SubsetInverseSampling::computeBlockSample()
{
  return NumericalSample();
}


/* Compute the new threshold corresponding to the target failure probability */
NumericalScalar SubsetInverseSampling::computeThreshold()
{
  // compute the quantile according to the event operator
  NumericalScalar ratio = getEvent().getOperator()(1.0, 2.0) ?  targetProbability_ : 1.0 - targetProbability_;
  
  NumericalScalar currentThreshold = currentLevelSample_.computeQuantile( ratio )[0];
  
  return currentThreshold;
}


NumericalScalar SubsetInverseSampling::computeProbability(NumericalScalar probabilityEstimateFactor, NumericalScalar threshold)
{
  const UnsignedInteger maximumOuterSampling = getMaximumOuterSampling();
  const UnsignedInteger blockSize = getBlockSize();
  NumericalScalar probabilityEstimate = 0.0;
  NumericalScalar varianceEstimate = 0.0;
  
  for ( UnsignedInteger i = 0; i < maximumOuterSampling; ++ i )
  {
    const NumericalScalar size = i + 1.0;
    NumericalScalar meanBlock = 0.0;
    NumericalScalar varianceBlock = 0.0;
    for ( UnsignedInteger j = 0 ; j < blockSize; ++ j )
    {
      if ( getEvent().getOperator()( currentLevelSample_[ i*blockSize+j ][0], threshold ) )
      {
        // update local mean and variance
        meanBlock += 1.0 / blockSize;
        varianceBlock += 1.0 * 1.0 / blockSize;
      }
    }
    varianceBlock -= pow( meanBlock, 2.0 );   
    
    // update global mean and variance
    varianceEstimate = (varianceBlock + (size - 1.0) * varianceEstimate) / size + (1.0 - 1.0 / size) * (probabilityEstimate - meanBlock) * (probabilityEstimate - meanBlock) / size;
    probabilityEstimate = (meanBlock + (size - 1.0) * probabilityEstimate) / size;
    
    // store convergence at each block
    NumericalPoint convergencePoint(2);
    convergencePoint[0] = probabilityEstimate * probabilityEstimateFactor;
    convergencePoint[1] = varianceEstimate * probabilityEstimateFactor * probabilityEstimateFactor / size;
    convergenceStrategy_.store(convergencePoint);
  }
  
  // cannot determine next subset domain if no variance
  const NumericalScalar epsilon = ResourceMap::GetAsNumericalScalar( "SpecFunc-Precision" );
  if ( fabs( varianceEstimate ) < epsilon )
    throw NotDefinedException(HERE) << "Null output variance";
                                          
  return probabilityEstimate;
}


/* Sort new seeds */
void SubsetInverseSampling::initializeSeed(NumericalScalar threshold)
{
  UnsignedInteger seedIndex = 0;
  const UnsignedInteger maximumOuterSampling = getMaximumOuterSampling();
  const UnsignedInteger blockSize = getBlockSize();
  for ( UnsignedInteger i = 0; i < maximumOuterSampling; ++ i )
  {
    for ( UnsignedInteger j = 0 ; j < blockSize; ++ j )
    {
      if ( getEvent().getOperator()( currentLevelSample_[ i*blockSize+j ][0], threshold ) )
      {
        // initialize seeds : they're grouped at the beginning of the sample
        currentPointSample_[ seedIndex ] = currentPointSample_[ i*blockSize+j ];
        currentLevelSample_[ seedIndex ] = currentLevelSample_[ i*blockSize+j ];
        ++ seedIndex;
      }
    }
  }
}


/* Compute the correlation on markov chains at the current state of the algorithm */
NumericalScalar SubsetInverseSampling::computeVarianceGamma(NumericalScalar currentFailureProbability, NumericalScalar threshold)
{
  const UnsignedInteger N = currentPointSample_.getSize();
  const UnsignedInteger Nc = std::max<UnsignedInteger>(1, targetProbability_ * N);
  Matrix IndicatorMatrice( Nc, N / Nc );
  NumericalPoint correlationSequence( N / Nc - 1 );
  NumericalScalar currentFailureProbability2 = pow( currentFailureProbability, 2. );
  for ( UnsignedInteger i = 0; i < N / Nc; ++ i )
  {
    for ( UnsignedInteger j = 0; j < Nc; ++ j )
    {
      IndicatorMatrice(j, i) = getEvent().getOperator()(currentLevelSample_[ i*Nc+j ][0], threshold);
    }
  }
  for ( UnsignedInteger k = 0; k < N / Nc - 1; ++ k )
  {
    for ( UnsignedInteger j = 0; j < Nc; ++ j )
    {
      for ( UnsignedInteger l = 0; l < N / Nc - k - 1; ++ l )
      {
        correlationSequence[k] += 1.0 * IndicatorMatrice(j, l) * IndicatorMatrice(j, l + (k + 1));
      }
    }
    correlationSequence[k] /= 1.0 * N - 1.0 * (k + 1) * Nc;
    correlationSequence[k] -= currentFailureProbability2;
  }
  const NumericalScalar R0 = currentFailureProbability * ( 1.0 - currentFailureProbability );
  NumericalPoint rho = ((1.0 / R0) * correlationSequence);
  NumericalScalar gamma = 0.0;
  for ( UnsignedInteger k = 0; k < N / Nc - 1; ++ k )
  {
    gamma += 2.0 * (1.0 - (k + 1) * 1.0 * Nc / N) * rho[k];
  }
  return gamma;
}


/* Iterate one step of the algorithm */
void SubsetInverseSampling::generatePoints(NumericalScalar threshold)
{  
  UnsignedInteger maximumOuterSampling = getMaximumOuterSampling();
  UnsignedInteger blockSize = getBlockSize();
  Distribution randomWalk(ComposedDistribution(ComposedDistribution::DistributionCollection(dimension_, Uniform(-0.5*proposalRange_, 0.5*proposalRange_))));
  UnsignedInteger N = currentPointSample_.getSize(); // total sample size
  UnsignedInteger Nc = targetProbability_ * N; //number of seeds (also = maximumOuterSampling*blockSize)
  
  for ( UnsignedInteger i = 0; i < maximumOuterSampling; ++ i )
  {    
    NumericalSample inputSample( blockSize, dimension_ );
    for ( UnsignedInteger j = 0; j < blockSize; ++ j )
    {
      // assign the new point to the seed, seed points being regrouped at the beginning of the sample
      if ( i*blockSize+j >= Nc )
      {
        currentPointSample_[ i*blockSize+j ] = currentPointSample_[ i*blockSize+j-Nc ];
        currentLevelSample_[ i*blockSize+j ] = currentLevelSample_[ i*blockSize+j-Nc ];     
      }
      
      // generate a new point
      NumericalPoint oldPoint( currentPointSample_[ i*blockSize+j ] );
      NumericalPoint newPoint( oldPoint + randomWalk.getRealization() );
      
      // 1. accept / reject new components
      NumericalPoint uniform( RandomGenerator::Generate(dimension_) );
      for (UnsignedInteger k = 0; k < dimension_; ++ k)
      {
        // compute ratio
        NumericalScalar ratio = std::min(1.0, exp(0.5 * (oldPoint[k] * oldPoint[k] - newPoint[k] * newPoint[k])));

        // accept new point with probability ratio
        if (ratio < uniform[k])
        {
          newPoint[k] = oldPoint[k];
        }
      }
      
      inputSample[j] = newPoint;
    }
   
    NumericalSample blockSample( standardEvent_.getFunction()( inputSample ) );

    for ( UnsignedInteger j = 0; j < getBlockSize(); ++ j )
    {
      // 2. accept the new point if in the failure domain
      if ( getEvent().getOperator()( blockSample[j][0], threshold ) )
      {
        currentPointSample_[ i*blockSize+j ] = inputSample[j];
        currentLevelSample_[ i*blockSize+j ] = blockSample[j];
      }
    }
  }
}

NumericalSample SubsetInverseSampling::getThresholdSample() const
{
  return sampleThreshold_;
}

/* Confidence Length of the threshold */
NumericalScalar SubsetInverseSampling::getThresholdConfidenceLength(const NumericalScalar level) const
{
  NumericalScalar thresholdInf = sampleThreshold_.computeQuantile((1 - level)/2)[0];
  NumericalScalar thresholdSup = sampleThreshold_.computeQuantile(level/2)[0];
  return NumericalScalar (std::max(thresholdSup, thresholdInf) - std::min(thresholdSup, thresholdInf));
}


/* Markov parameter accessor */
void SubsetInverseSampling::setProposalRange(NumericalScalar proposalRange)
{
  proposalRange_ = proposalRange;
}


NumericalScalar SubsetInverseSampling::getProposalRange() const
{
  return proposalRange_;
}


/* Ratio accessor */
void SubsetInverseSampling::setTargetProbability(NumericalScalar targetProbability)
{
  if ( (targetProbability <= 0.) || (targetProbability >= 1.) ) throw InvalidArgumentException(HERE) << "Probability should be in (0, 1)";
  targetProbability_ = targetProbability;
}

NumericalScalar SubsetInverseSampling::getTargetProbability() const
{
  return targetProbability_;
}

/* final target probability accessor */
void SubsetInverseSampling::setFinalTargetProbability(NumericalScalar finalTargetProbability)
{
  if ( (finalTargetProbability <= 0.) || (finalTargetProbability >= 1.) ) throw InvalidArgumentException(HERE) << "Probability should be in (0, 1)";
  finalTargetProbability_ = finalTargetProbability;
}

NumericalScalar SubsetInverseSampling::getFinalTargetProbability() const
{
  return finalTargetProbability_;
}

UnsignedInteger SubsetInverseSampling::getNumberOfSteps()
{
  return numberOfSteps_;
}


OT::NumericalPoint SubsetInverseSampling::getGammaPerStep() const
{
  return gammaPerStep_;
}


OT::NumericalPoint SubsetInverseSampling::getCoefficientOfVariationPerStep() const
{
  return coefficientOfVariationPerStep_;
}


OT::NumericalPoint SubsetInverseSampling::getProbabilityEstimatePerStep() const
{
  return probabilityEstimatePerStep_;
}


OT::NumericalPoint SubsetInverseSampling::getThresholdPerStep() const
{
  return thresholdPerStep_;
}


void SubsetInverseSampling::setKeepEventSample(bool keepEventSample)
{
  keepEventSample_ = keepEventSample;
}


NumericalSample SubsetInverseSampling::getEventInputSample() const
{
  return eventInputSample_;
}


NumericalSample SubsetInverseSampling::getEventOutputSample() const
{
  return eventOutputSample_;
}


void SubsetInverseSampling::setISubset(OT::Bool iSubset)
{
  iSubset_ = iSubset;
}

void SubsetInverseSampling::setBetaMin(NumericalScalar betaMin)
{
  if (betaMin <= 0.) throw InvalidArgumentException(HERE) << "Beta min should be positive";
  betaMin_ = betaMin;
}




/* String converter */
String SubsetInverseSampling::__repr__() const
{
  OSS oss;
  oss << "class=" << getClassName()
      << " derived from " << Simulation::__repr__()
      << " finalTargetProbability=" << finalTargetProbability_
      << " proposalRange=" << proposalRange_
      << " targetProbability=" << targetProbability_
      << " keepEventSample_=" << keepEventSample_;
  return oss;
}


/* Method save() stores the object through the StorageManager */
void SubsetInverseSampling::save(Advocate & adv) const
{
  Simulation::save(adv);
  adv.saveAttribute("finalTargetProbability", finalTargetProbability_);
  adv.saveAttribute("proposalRange_", proposalRange_);
  adv.saveAttribute("targetProbability_", targetProbability_);
  adv.saveAttribute("iSubset_", iSubset_);
  adv.saveAttribute("betaMin_", betaMin_);
  adv.saveAttribute("keepEventSample_", keepEventSample_);  
  
  adv.saveAttribute("numberOfSteps_", numberOfSteps_);
  adv.saveAttribute("thresholdPerStep_", thresholdPerStep_);
  adv.saveAttribute("gammaPerStep_", gammaPerStep_);
  adv.saveAttribute("coefficientOfVariationPerStep_", coefficientOfVariationPerStep_);
  adv.saveAttribute("probabilityEstimatePerStep_", probabilityEstimatePerStep_);
}


/* Method load() reloads the object from the StorageManager */
void SubsetInverseSampling::load(Advocate & adv)
{
  Simulation::load(adv);
  adv.loadAttribute("finalTargetProbability", finalTargetProbability_);
  adv.loadAttribute("proposalRange_", proposalRange_);
  adv.loadAttribute("targetProbability_", targetProbability_);
  adv.loadAttribute("keepEventSample_", keepEventSample_);
  adv.loadAttribute("iSubset_", iSubset_);
  adv.loadAttribute("betaMin_", betaMin_);
  
  adv.loadAttribute("numberOfSteps_", numberOfSteps_);
  adv.loadAttribute("thresholdPerStep_", thresholdPerStep_);
  adv.loadAttribute("gammaPerStep_", gammaPerStep_);
  adv.loadAttribute("coefficientOfVariationPerStep_", coefficientOfVariationPerStep_);
  adv.loadAttribute("probabilityEstimatePerStep_", probabilityEstimatePerStep_);
}




} /* namespace OTSubsetInverse */
