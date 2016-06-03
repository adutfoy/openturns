import openturns as ot

# Create the Normal distribution
mu = 0.0
sigma = 1.0
myNormalDist = ot.Normal(mu, sigma)

# Get a realization
myReal = myNormalDist.getRealization()

# Sample it!
N = 5
mySample = myNormalDist.getSample(N)
print(mySample)
