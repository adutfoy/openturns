//                                               -*- C++ -*-
/**
 * @brief Class for identity evaluation
 *
 *  Copyright 2005-2018 Airbus-EDF-IMACS-Phimeca
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
 *  You should have received a copy of the GNU Lesser General Public
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef OPENTURNS_MEMOIZEEVALUATION_HXX
#define OPENTURNS_MEMOIZEEVALUATION_HXX

#include "openturns/EvaluationProxy.hxx"
#include "openturns/Evaluation.hxx"
#include "openturns/HistoryStrategy.hxx"
#include "openturns/Full.hxx"
#include "openturns/Point.hxx"
#include "openturns/Sample.hxx"
#include "openturns/Indices.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class MemoizeEvaluation
 */
class OT_API MemoizeEvaluation
  : public EvaluationProxy
{
  CLASSNAME
public:

  /** Default constructor */
  MemoizeEvaluation();

  /** Parameter constructor */
  explicit MemoizeEvaluation(const Evaluation & evaluation, const HistoryStrategy & historyStrategy = Full());

  /** Virtual constructor */
  virtual MemoizeEvaluation * clone() const;

  /** Function implementation accessors */
  void setEvaluation(const Evaluation & evaluation);
  Evaluation getEvaluation() const;

  /** Comparison operator */
  Bool operator ==(const MemoizeEvaluation & other) const;

  /** String converter */
  virtual String __repr__() const;
  virtual String __str__(const String & offset = "") const;

  /* Here is the interface that all derived class must implement */

  /** Operator () */
  virtual Point operator() (const Point & inPoint) const;

  /** Operator () */
  virtual Sample operator() (const Sample & inSample) const;

  /** Get the evaluation corresponding to indices components */
  virtual Evaluation getMarginal(const Indices & indices) const;
  using EvaluationImplementation::getMarginal;

  /** Enable or disable the input/output history */
  void enableHistory() const;
  void disableHistory() const;

  /** Test the history mechanism activity */
  Bool isHistoryEnabled() const;

  /** Clear history of the input and output values */
  void clearHistory() const;

  /** Retrieve the history of the input values */
  Sample getInputHistory() const;

  /** Retrieve the history of the output values */
  Sample getOutputHistory() const;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv);

private:

  /** An history mechanism that allows to remember all the input/output associations including duplicate calls */
  mutable HistoryStrategy inputStrategy_;
  mutable HistoryStrategy outputStrategy_;

  /** Flag to activate or deactivate the history mechanism */
  mutable Bool isHistoryEnabled_;

}; /* class MemoizeEvaluation */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_MEMOIZEEVALUATION_HXX */