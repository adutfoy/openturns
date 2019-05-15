import openturns as ot
from matplotlib import pyplot as plt
from openturns.viewer import View

R = ot.CorrelationMatrix(2)
R[0, 1] = 0.8
d1 =  ot.Normal([0.0, 0.0], R)

R2 = ot.CorrelationMatrix(2)
R2[0, 1] = -0.8
d2 = ot.Normal([0.0, 0.0], R2)
collection = [d1, d2]
d2 =  ot.Normal([0.0, 0.0], R2)

dist = ot.BlockIndependentDistribution([d1, d2])
sample = dist.getSample(10000)

g12 = ot.Graph()
g12.setAxes(True)
cloud = ot.Cloud(sample.getMarginal([0,1]).rank())
cloud.setPointStyle('dot')
g12.add(cloud)
g12.setLegends([r'$(u_0, u_1)$'])
g12.setLegendPosition('topright')

g34 = ot.Graph()
g34.setAxes(True)
cloud = ot.Cloud(sample.getMarginal([2,3]).rank())
cloud.setPointStyle('dot')
g34.add(cloud)
g34.setLegends([r'$(u_2, u_3)$'])
g34.setLegendPosition('topright')

g13 = ot.Graph()
g13.setAxes(True)
cloud = ot.Cloud(sample.getMarginal([0,2]).rank())
cloud.setPointStyle('dot')
g13.add(cloud)
g13.setLegends([r'$((u_0, u_2)$'])
g13.setLegendPosition('topright')

g24 = ot.Graph()
g24.setAxes(True)
cloud = ot.Cloud(sample.getMarginal([1,3]).rank())
cloud.setPointStyle('dot')
g24.add(cloud)
g24.setLegends([r'$((u_1, u_3)$'])
g24.setLegendPosition('topright')

fig = plt.figure(figsize=(10, 10))
g12_axis = fig.add_subplot(221)
g34_axis = fig.add_subplot(222)
g13_axis = fig.add_subplot(223)
g24_axis = fig.add_subplot(224)
g12_axis.set_xlim(auto=True)
g34_axis.set_xlim(auto=True)
g13_axis.set_xlim(auto=True)
g24_axis.set_xlim(auto=True)

View(g12, figure=fig, axes=[g12_axis], add_legend=True)
View(g34, figure=fig, axes=[g34_axis], add_legend=True)
View(g13, figure=fig, axes=[g13_axis], add_legend=True)
View(g24, figure=fig, axes=[g24_axis], add_legend=True)
fig.suptitle("IndependentBolckDistrbution(2d normal, 2d normal) in the rank space: 2d samples")
