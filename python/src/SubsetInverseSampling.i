// SWIG file

%{
#include "otsubsetinverse/SubsetInverseSampling.hxx"
%}

%include SubsetInverseSampling_doc.i

%include otsubsetinverse/SubsetInverseSampling.hxx
namespace OTSubsetInverse { %extend SubsetInverseSampling { SubsetInverseSampling(const SubsetInverseSampling & other) { return new OTSubsetInverse::SubsetInverseSampling(other); } } }
