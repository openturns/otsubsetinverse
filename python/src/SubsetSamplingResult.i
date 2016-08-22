// SWIG file

%{
#include "otsubsetinverse/SubsetSamplingResult.hxx"
%}

%include otsubsetinverse/SubsetSamplingResult.hxx
namespace OTSubsetInverse { %extend SubsetSamplingResult { SubsetSamplingResult(const SubsetSamplingResult & other) { return new OTSubsetInverse::SubsetSamplingResult(other); } } }
