from plotly.graph_objects import Figure, Contour

def plotcontour(x, y, u, size = (700, 700)):
    a, b = size
    fig = Figure(data = Contour(z = u, x = x, y = y))
    fig.update_layout({
    "title": {"text": r"$\text{Maximum power } \frac{P_{max}}{N\cdot\omega_z^2}$", "x": 0.5, "y": 0.9, "font": {"size": 14} },
    "showlegend": True,
    "xaxis": {"title": "$g$", "showticklabels": True, "dtick": 0.2}, 
    "yaxis": {"title": "$N$", "showticklabels": True, "dtick": 1}, 
    "autosize": False, 
    "width":a, 
    "height":b})
    fig.show()