import streamlit as st
import plotly.graph_objs as go
import networkx as nx
import numpy as np
from gudhi import RipsComplex

# Distance matrix (sG3) and labels from the H7N9 triple reassortment validation
sG3 = np.array([
    [    0.,  639., 1303., 1220., 2027., 2092.],
    [  639.,    0., 1452., 1634., 2150., 2026.],
    [ 1303., 1452.,    0., 2080., 1673., 1982.],
    [ 1220., 1634., 2080.,    0., 1092., 1690.],
    [ 2027., 2150., 1673., 1092.,    0., 1526.],
    [ 2092., 2026., 1982., 1690., 1526.,    0.]
])

labels = [
    'A/brambling/Beijing/16/2012(H9N2)',
    'A/quail/Wuxi/7/2010(H9N2)',
    'A/quail/Lebanon/273/2010(H9N2)',
    'A/Shanghai/02/2013(H7N9)',
    'A/wild bird/Korea/A14/2011(H7N9)',
    'A/duck/Zhejiang/10/2011(H7N3)',
]

# Create Rips complex and compute persistence
max_edge_length = 2300
max_dimension = 5
rips_complex = RipsComplex(distance_matrix=sG3, max_edge_length=max_edge_length)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dimension)
simplex_tree.persistence()

# Predefined octahedral positions for visualization
octahedron_pos = {
    0: (0, 0, 1),    # A/brambling/Beijing (top)
    1: (1, 1, 0),    # A/quail/Wuxi
    2: (1, -1, 0),   # A/quail/Lebanon
    3: (-1, 1, 0),   # A/Shanghai
    4: (-1, -1, 0),  # A/wild bird/Korea/A14
    5: (0, 0, -1)    # A/duck/Zhejiang (bottom)
}

# Function to create 3D visualization
def create_3d_visualization(simplex_tree, vertex_names, filtration_value):
    G = nx.Graph()
    triangles = []
    tetrahedra = []
    hole_2d_simplices = []

    # Identify 2D holes for highlighting
    persistence_2 = simplex_tree.persistence_intervals_in_dimension(2)
    hole_birth = persistence_2[0][0] if len(persistence_2) > 0 else None

    # Process simplices based on filtration value
    for simplex, filt in simplex_tree.get_filtration():
        if filt <= filtration_value:
            if len(simplex) == 2:
                G.add_edge(simplex[0], simplex[1])
            elif len(simplex) == 1:
                G.add_node(simplex[0])
            elif len(simplex) == 3:
                triangles.append(simplex)
                if hole_birth is not None and abs(filt - hole_birth) < 1e-6:
                    hole_2d_simplices.append(simplex)
            elif len(simplex) == 4:
                tetrahedra.append(simplex)

    # Use predefined positions
    pos = octahedron_pos

    # Node traces
    node_x, node_y, node_z = zip(*[pos[node] for node in G.nodes()])
    node_labels = [vertex_names[node] for node in G.nodes()]
    node_colors = ['blue' if 'H9N2' in label else 'orange' if 'H7N9' in label else 'green' for label in node_labels]
    node_trace = go.Scatter3d(
        x=node_x, y=node_y, z=node_z, mode='markers+text',
        marker=dict(size=5, color=node_colors),
        text=node_labels, textposition="top center",
        hoverinfo='text'
    )

    # Edge traces
    edge_x, edge_y, edge_z = [], [], []
    edge_colors = []
    for edge in G.edges():
        x0, y0, z0 = pos[edge[0]]
        x1, y1, z1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])
        is_hole_edge = any(set(edge).issubset(triangle) for triangle in hole_2d_simplices)
        edge_colors.extend(['yellow' if is_hole_edge else 'black'] * 3)
    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode='lines',
        line=dict(width=2, color=edge_colors), hoverinfo='none'
    )

    # Triangle traces
    triangle_traces = []
    for triangle in triangles:
        x = [pos[v][0] for v in triangle]
        y = [pos[v][1] for v in triangle]
        z = [pos[v][2] for v in triangle]
        color = 'yellow' if triangle in hole_2d_simplices else 'red'
        opacity = 0.7 if triangle in hole_2d_simplices else 0.3
        triangle_trace = go.Mesh3d(
            x=x, y=y, z=z,
            color=color, opacity=opacity,
            i=[0], j=[1], k=[2]
        )
        triangle_traces.append(triangle_trace)

    # Tetrahedron traces
    tetrahedra_traces = []
    for tetrahedron in tetrahedra:
        x = [pos[v][0] for v in tetrahedron]
        y = [pos[v][1] for v in tetrahedron]
        z = [pos[v][2] for v in tetrahedron]
        tetrahedron_trace = go.Mesh3d(
            x=x, y=y, z=z,
            opacity=0.1, color='green',
            i=[0, 1, 0, 1], j=[1, 2, 2, 3], k=[2, 3, 3, 0]
        )
        tetrahedra_traces.append(tetrahedron_trace)

    # Create figure
    fig = go.Figure(data=[node_trace, edge_trace] + triangle_traces + tetrahedra_traces)
    fig.update_layout(
        title=f"Simplicial Complex Visualization (Filtration: {filtration_value})",
        scene=dict(
            xaxis_title="X", yaxis_title="Y", zaxis_title="Z",
            xaxis=dict(showgrid=False), yaxis=dict(showgrid=False), zaxis=dict(showgrid=False)
        ),
        width=800, height=600, showlegend=False
    )
    return fig

# Streamlit interface
st.title("Persistent Homology in Viral Genetic Sequences")

# English description
st.write("""
### Visualization of Topological Structures in Viral Evolution
This tool implements a 3D visualization of simplicial complexes derived from genetic distance matrices of viral sequences, as developed in the thesis "Identification of Recombination in Viruses Using Topological Data Analysis" by Jaime Salvador López Viveros. Based on persistent homology, it transforms genomic sequences into topological structures to detect horizontal gene transfer events, such as reassortment and recombination. Here, we replicate the validation experiment for the H7N9 triple reassortment (Section 5.1.1), where a distance matrix of six concatenated sequences reveals a 1-hole (H1, birth: 1673, death: 2026) and a 2-hole (H2, birth: 2027, death: 2080), indicative of reassortment among H9N2, H7N9, and H7N3 strains.

- **Nodes**: Represent viral sequences (H9N2: blue, H7N9: orange, H7N3: green).
- **Edges**: Genetic distances forming 1-simplices (yellow highlights edges of 2-holes).
- **Triangles**: 2-simplices (red, yellow for persistent 2-holes).
- **Tetrahedra**: 3-simplices (green, rare in this dataset).

Adjust the filtration value using the slider or text box to explore how topological features evolve, reflecting the formation and collapse of cycles (H1) and cavities (H2) that signal evolutionary events.
""")

# Spanish description
st.write("""
### Visualización de Estructuras Topológicas en la Evolución Viral
Esta herramienta implementa una visualización 3D de complejos simpliciales derivados de matrices de distancia genética de secuencias virales, desarrollada en la tesis "Identificación de Recombinación en Virus Utilizando Análisis Topológico de Datos" por Jaime Salvador López Viveros. Basada en homología persistente, transforma secuencias genómicas en estructuras topológicas para detectar eventos de transferencia genética horizontal, como reordenamiento y recombinación. Aquí replicamos el experimento de validación del reordenamiento triple en H7N9 (Sección 5.1.1), donde una matriz de distancia de seis secuencias concatenadas revela un 1-hoyo (H1, nacimiento: 1673, muerte: 2026) y un 2-hoyo (H2, nacimiento: 2027, muerte: 2080), indicativos de reordenamiento entre cepas H9N2, H7N9 y H7N3.

- **Nodos**: Representan secuencias virales (H9N2: azul, H7N9: naranja, H7N3: verde).
- **Aristas**: Distancias genéticas que forman 1-símplices (amarillo resalta aristas de 2-hoyos).
- **Triángulos**: 2-símplices (rojo, amarillo para 2-hoyos persistentes).
- **Tetraedros**: 3-símplices (verde, raros en este conjunto).

Ajusta el valor de filtración con el deslizador o el cuadro de texto para explorar cómo evolucionan las características topológicas, reflejando la formación y colapso de ciclos (H1) y cavidades (H2) que señalan eventos evolutivos.
""")

# Synchronized slider and text input
if 'filtration_value' not in st.session_state:
    st.session_state.filtration_value = 0

col1, col2 = st.columns([3, 1])
with col1:
    filtration_value_slider = st.slider(
        "Filtration Value (Slider) / Valor de Filtración (Deslizador)",
        min_value=0,
        max_value=2080,
        value=st.session_state.filtration_value,
        step=1,
        key="slider",
        on_change=lambda: st.session_state.update(filtration_value=st.session_state.slider)
    )

with col2:
    filtration_value_input = st.number_input(
        "Filtration Value (Manual) / Valor de Filtración (Manual)",
        min_value=0,
        max_value=2080,
        value=st.session_state.filtration_value,
        step=1,
        key="input",
        on_change=lambda: st.session_state.update(filtration_value=st.session_state.input)
    )

filtration_value = st.session_state.filtration_value

# Generate and display visualization
fig = create_3d_visualization(simplex_tree, labels, filtration_value)
st.plotly_chart(fig, use_container_width=True)

# Persistence intervals
st.write("### Persistence Intervals / Intervalos de Persistencia")
st.write("**1-Dimensional Holes (Cycles) / Agujeros 1-Dimensionales (Ciclos):**")
st.write("Birth / Nacimiento: 1673.0, Death / Muerte: 2026.0")
st.write("**2-Dimensional Holes (Cavities) / Agujeros 2-Dimensionales (Cavidades):**")
st.write("Birth / Nacimiento: 2027.0, Death / Muerte: 2080.0")