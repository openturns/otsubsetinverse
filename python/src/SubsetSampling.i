// SWIG file

%{
#include "otsubsetinverse/SubsetSampling.hxx"
%}

%include otsubsetinverse/SubsetSampling.hxx
namespace OTSubsetInverse { %extend SubsetSampling { SubsetSampling(const SubsetSampling & other) { return new OTSubsetInverse::SubsetSampling(other); } } }
