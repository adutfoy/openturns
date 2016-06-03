import openturns as ot

# Create the objective and constraint functions
objFunc = ot.NumericalMathFunction(['x1', 'x2'], ['2*(x1-2)^2+3*(x2-2)^2'])
constFunc = ot.NumericalMathFunction(['x1', 'x2'], ['4-(x1^2+x2^2)'])

# Create the optimization problem
myOptimPb = ot.OptimizationProblem()
myOptimPb.setObjective(objFunc)

# Add the constraint
myOptimPb.setInequalityConstraint(constFunc)

# Add the bounds
myBounds = ot.Interval([-5.0]*2, [5.0]*2)
myOptimPb.setBounds(myBounds)dr=g2.getDrawable(0)

# Solve it
mySolver = ot.Cobyla(myOptimPb)
mySolver.setStartingPoint([0.0]*2)
mySolver.run()

# Get the whole result!
result = mySolver.getResult()

# Get the optimal point
print result.getOptimalPoint()

# Get the value of the objective function and the constraint function
# at the optimal point
print result.getOptimalValue()



