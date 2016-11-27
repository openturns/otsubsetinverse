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
#ifndef OTSUBSETINVERSE_SUBSETINVERSESAMPLING_HXX
#define OTSUBSETINVERSE_SUBSETINVERSESAMPLING_HXX

// #include <OT.hxx>
#include <openturns/Simulation.hxx>
#include <openturns/StandardEvent.hxx>
#include "otsubsetinverse/OTSubsetInverseprivate.hxx"
#include "openturns/Collection.hxx"
#include "openturns/PersistentCollection.hxx"

namespace OTSubsetInverse
{


class OTSUBSETINVERSE_API SubsetInverseSampling
: public OT::Simulation
{
CLASSNAME
public:

  typedef OT::Collection<OT::NumericalSample> NumericalSampleCollection;
  typedef OT::PersistentCollection<OT::NumericalSample> NumericalSamplePersistentCollection;

  /** Default Parameters */
  static const OT::UnsignedInteger DefaultMaximumOuterSampling;
  static const OT::NumericalScalar DefaultProposalRange;
  static const OT::NumericalScalar DefaultConditionalProbability;
  static const OT::NumericalScalar DefaultBetaMin;

  /** Default Constructor */
  SubsetInverseSampling();


  /** Constructor with parameters */
  SubsetInverseSampling(const OT::Event & event,
                 const OT::NumericalScalar targetProbability,
                 const OT::NumericalScalar proposalRange = DefaultProposalRange,
                 const OT::NumericalScalar conditionalProbability = DefaultConditionalProbability);

  /** Virtual constructor */
  virtual SubsetInverseSampling * clone() const;

  /** The range of the uniform proposal pdf */
  void setProposalRange(OT::NumericalScalar proposalRange);
  OT::NumericalScalar getProposalRange() const;

  /** Ratio parameter */
  void setConditionalProbability(OT::NumericalScalar conditionalProbability);
  OT::NumericalScalar getConditionalProbability() const;

  /** final target probability */
  void setTargetProbability(OT::NumericalScalar targetProbability);
  OT::NumericalScalar getTargetProbability() const;
  
  /** Accessor to the achieved number of steps */
  OT::UnsignedInteger getNumberOfSteps();

  OT::NumericalScalar getThresholdConfidenceLength(const OT::NumericalScalar level = OT::ResourceMap::GetAsNumericalScalar( "SimulationResultImplementation-DefaultConfidenceLevel" )) const;
    
  /** Stepwise result accessors */
  OT::NumericalPoint getThresholdPerStep() const;
  OT::NumericalPoint getGammaPerStep() const;
  OT::NumericalPoint getCoefficientOfVariationPerStep() const;
  OT::NumericalPoint getProbabilityEstimatePerStep() const;
  OT::NumericalPoint getThresholdCoefficientOfVariationPerStep() const;
  
  /** Keep event sample */
  void setKeepEventSample(bool keepEventSample);
  
  /** Event input/output sample accessor */
  OT::NumericalSample getEventInputSample() const;
  OT::NumericalSample getEventOutputSample() const;

  /** All level sample accessor*/
  NumericalSampleCollection getOutputSample() const;
  NumericalSampleCollection getInputSample() const;

  /** i-subset */
  void setISubset(OT::Bool iSubset);
  void setBetaMin(OT::NumericalScalar betaMin);
  
  /** Performs the actual computation. */
  void run();

  /** String converter */
  OT::String __repr__() const;

  /** Method save() stores the object through the StorageManager */
  virtual void save(OT::Advocate & adv) const;

  /** Method load() reloads the object from the StorageManager */
  virtual void load(OT::Advocate & adv);

private:
  /** Compute the block sample */
  OT::NumericalSample computeBlockSample();

  /** Compute the new threshold corresponding to the conditional failure probability */
  OT::NumericalScalar computeThreshold();

  /** compute probability estimate on the current sample */
  OT::NumericalScalar computeProbability(OT::NumericalScalar probabilityEstimate, OT::NumericalScalar threshold);

  /** Sort new seeds */
  void initializeSeed(OT::NumericalScalar threshold);
  
  /** Compute the correlation on markov chains at the current state of the algorithm */
  OT::NumericalScalar computeVarianceGamma(OT::NumericalScalar currentFailureProbability, OT::NumericalScalar threshold);

  /** Generate new points in the conditional failure domain */
  void generatePoints(OT::NumericalScalar threshold);

  // some parameters
  OT::NumericalScalar proposalRange_;// width of the proposal pdf
  OT::NumericalScalar conditionalProbability_;// conditional probability at each subset
  OT::Bool iSubset_;// conditional pre-sampling
  OT::NumericalScalar betaMin_;// pre-sampling hypersphere exclusion radius
  OT::Bool keepEventSample_;// do we keep the event sample ?
  OT::NumericalScalar targetProbability_;// final target probability

  // some results
  OT::UnsignedInteger numberOfSteps_;// number of subset steps
  OT::NumericalPoint thresholdPerStep_;// intermediate thresholds
  OT::NumericalPoint gammaPerStep_;// intermediate gammas
  OT::NumericalPoint coefficientOfVariationPerStep_;// intermediate COVS
  OT::NumericalPoint probabilityEstimatePerStep_;// intermediate PFs
  OT::NumericalPoint thresholdCoefficientOfVariationPerStep_;// intermediate threshold COVs
  OT::NumericalSample eventInputSample_;// event input sample
  OT::NumericalSample eventOutputSample_;// event output sample
  NumericalSamplePersistentCollection allLevelSample_; // all event output sample
  NumericalSamplePersistentCollection allPointSample_; // all event output sample
  OT::Distribution thresholdDistribution_;//distribution of the final threshold value

  // attributes used for conveniency, not to be saved/loaded
  OT::StandardEvent standardEvent_;// the algorithm happens in U
  OT::UnsignedInteger dimension_;// input dimension
  OT::NumericalSample currentPointSample_;// X
  OT::NumericalSample currentLevelSample_;//f(X)

}; /* class SubsetInverseSampling */

} /* namespace OTSubsetInverse */

#endif /* OTSUBSETINVERSE_SUBSETINVERSESAMPLING_HXX */
