import openturns as ot
from matplotlib import pyplot as plt
from openturns.viewer import View

locations = ot.RegularGrid(0.0, 0.25, 11).getVertices().asPoint()
f = ot.SymbolicFunction(['x'], ['sin(pi_*x)'])
values = f(ot.Sample.BuildFromPoint(locations))
derivatives = [f.gradient([x])[0,0] for x in locations]
func_H1 = ot.PiecewiseHermiteEvaluation(locations, values.asPoint(), derivatives)

graph = ot.Graph(r'Piecewise Hermite interpolation of $x  \mapsto \sin(\pi x)$', r'$x$', r'$y$', True) 
cloud = ot.Cloud(ot.Sample.BuildFromPoint(locations), values)
cloud.setPointStyle("square")
graph.add(cloud)
draw_f = f.draw(-0.5, 3., 1024).getDrawable(0)
draw_f.setLineStyle('dashed')
graph.add(draw_f)
graph.add(func_H1.draw(-0.5, 3., 1024))
graph.setColors(['lightblue', "blue", "red"])

fig = plt.figure(figsize=(10, 4))
func_axis = fig.add_subplot(111)
view = View(graph, figure=fig, axes=[func_axis], add_legend=False)
