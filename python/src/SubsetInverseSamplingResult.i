// SWIG file

%{
#include "otsubsetinverse/SubsetInverseSamplingResult.hxx"
%}

%include otsubsetinverse/SubsetInverseSamplingResult.hxx
namespace OTSubsetInverse { %extend SubsetInverseSamplingResult { SubsetInverseSamplingResult(const SubsetInverseSamplingResult & other) { return new OTSubsetInverse::SubsetInverseSamplingResult(other); } } }
