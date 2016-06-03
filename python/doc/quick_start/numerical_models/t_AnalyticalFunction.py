import openturns as ot

# Create the fucntion from its analytical expression
f = ot.NumericalMathFunction(['x_1', 'x_2'], ['y'], ['cos(sqrt(x_1^2+x_2^2))'])

# Compute f on a point
print f([0.0, 1.0])

# Compute f on a sample of saize 2
mySample = [[0.0, 1.0], [2.0, 3.0]] 
print f(mySample)

# Get its gradients
print f.getGradient()

# Draw its isovalues
graph = f.draw([-3.0]*2, [3.0]*2)
