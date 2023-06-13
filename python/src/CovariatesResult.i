// SWIG file CovariatesResult.i

%{
#include "openturns/CovariatesResult.hxx"
%}

%include CovariatesResult_doc.i

%include openturns/CovariatesResult.hxx
namespace OT { %extend CovariatesResult { CovariatesResult(const CovariatesResult & other) { return new OT::CovariatesResult(other); } } }
