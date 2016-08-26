// SWIG file otsubsetinverse_module.i

%module(docstring="otsubsetinverse module") otsubsetinverse

%{
#include <OT.hxx>
#include <PythonWrappingFunctions.hxx>
%}

// Prerequisites needed
%include typemaps.i
%include exception.i
%ignore *::load(OT::Advocate & adv);
%ignore *::save(OT::Advocate & adv) const;

%import base_module.i
%import uncertainty_module.i

// The new classes
%include otsubsetinverse/OTSubsetInverseprivate.hxx
%include SubsetInverseSamplingResult.i
%include SubsetInverseSampling.i


