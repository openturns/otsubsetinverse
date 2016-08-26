//                                               -*- C++ -*-
/**
 *  @brief SubsetInverseSamplingResult
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
#ifndef OTSUBSETINVERSE_SUBSETINVERSESAMPLINGRESULT_HXX
#define OTSUBSETINVERSE_SUBSETINVERSESAMPLINGRESULT_HXX

#include <OT.hxx>
#include "otsubsetinverse/OTSubsetInverseprivate.hxx"

namespace OTSubsetInverse
{

class OTSUBSETINVERSE_API SubsetInverseSamplingResult
: public OT::SimulationResultImplementation
{
  CLASSNAME;
public:
  
  /** Constructor with parameters */
  SubsetInverseSamplingResult();

  /** Constructor with parameters */
  SubsetInverseSamplingResult(const OT::Event & event,
                       const OT::NumericalScalar probabilityEstimate,
                       const OT::NumericalScalar varianceEstimate,
                       const OT::UnsignedInteger outerSampling,
                       const OT::UnsignedInteger blockSize,
                       const OT::NumericalScalar coefficientOfVariation = 0.0,
                       const OT::NumericalScalar threshold = 0.0);

  /** Virtual constructor */
  virtual SubsetInverseSamplingResult * clone() const;

  /** Coefficient of variation estimate accessor */
  virtual OT::NumericalScalar getCoefficientOfVariation() const;

  /** String converter */
  virtual OT::String __repr__() const;

  /** Method save() stores the object through the StorageManager */
  virtual void save(OT::Advocate & adv) const;

  /** Method load() reloads the object from the StorageManager */
  virtual void load(OT::Advocate & adv);

protected:
  OT::NumericalScalar coefficientOfVariation_;
  OT::NumericalScalar threshold_;

private:

}; /* class SubsetInverseSamplingResult */

} /* namespace OTSubsetInverse */

#endif /* OTSUBSETINVERSE_SUBSETINVERSESAMPLINGRESULT_HXX */
