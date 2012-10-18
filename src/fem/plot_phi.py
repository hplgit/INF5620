import scitools.std as plt
from fe_approx1D import mesh

def plot_fe_mesh(nodes, elements, element_marker=[0, 0.1]):
    plt.hold('on')
    all_x_L = [nodes[elements[e][0]] for e in range(len(elements))]
    element_boundaries = all_x_L + [nodes[-1]]
    for x in element_boundaries:
        plt.plot([x, x], element_marker, 'm--')  # m gives dotted eps/pdf lines
    plt.plot(nodes, [0]*len(nodes), 'ro2')

def fe_basis_function_figure(d, target_elm=[1], n_e=3,
                             filename='tmp.pdf'):
    """
    Draw all basis functions, of degree d, associated with
    element target_elm (may be list of elements).  Add a mesh with n_e
    elements.
    """
    nodes, elements = mesh(n_e, d)
    """
    x = 1.1
    print locate_element_vectorized(x, elements, nodes)
    print locate_element_scalar(x, elements, nodes)
    x = 0.1, 0.4, 0.8
    print locate_element_vectorized(x, elements, nodes)
    """
    if isinstance(target_elm, int):
        target_elm = [target_elm]  # wrap in list

    # Draw the basis functions for element 1
    phi_drawn = []  # list of already drawn phi functions
    ymin = ymax = 0
    for e in target_elm:
        for i in elements[e]:
            if not i in phi_drawn:
                x, y = phi_glob(i, elements, nodes)
                ymax = max(ymax, max(y))
                ymin = min(ymin, min(y))
                plt.plot(x, y, '-')
                plt.hold('on')
                #legend(r'\phi_%d' % i)
                phi_drawn.append(i)

    plt.axis([nodes[0], nodes[-1], ymin-0.1, ymax+0.1])
    plot_fe_mesh(nodes, elements, element_marker=[ymin-0.1, ymax+0.1])
    plt.hold('off')
    plt.savefig(filename)

def draw_basis_functions():
    for ext in 'pdf', 'png':
        fe_basis_function_figure(d=1, target_elm=1, n_e=4,
                                 filename='fe_basis_p1_4e.%s' % ext)
        figure()
        fe_basis_function_figure(d=2, target_elm=1, n_e=4,
                                 filename='fe_basis_p2_4e.%s' % ext)
        figure()
        fe_basis_function_figure(d=1, target_elm=1, n_e=5,
                                 filename='fe_basis_p1_5e.%s' % ext)
        figure()
        fe_basis_function_figure(d=3, target_elm=1, n_e=4,
                                 filename='fe_basis_p3_4e.%s' % ext)
        figure()
        fe_basis_function_figure(d=3, target_elm=2, n_e=5,
                                 filename='fe_basis_p3_5e.%s' % ext)

def draw_sparsity_pattern(elements, num_nodes):
    import matplotlib.pyplot as plt
    sparsity_pattern = {}
    for e in elements:
        for i in range(len(e)):
            for j in range(i, len(e)):
                sparsity_pattern[(i,j)] = 1
    x = [i for i, j in sparsity_pattern]
    y = [j for i, j in sparsity_pattern]
    plt.pyplot(x, y, 'bo')
    plt.axis([0, num_nodes-1, 0, num_nodes])
    plt.savefig('tmp.pdf')
    plt.savefig('tmp.png')

if __name__ == '__main__':
    if sys.argv[1] == 'phi':
        draw_basis_functions()
    elif sys.argv[1] == 'pattern':
        num_elements = int(sys.argv[2])
        d = int(sys.argv[3])
        uniform = sys.argv[4]
        nodes, elements = mesh(num_elements, d, [0, 1])
        num_nodes = len(nodes)

        # Shuffle nodes in random order if necessary
        if uniform == 'random':
            global_node_numbering = range(0, num_nodes)
            import random
            random.shuffl(global_node_numbering)
            for e in range(len(elements)):
                for r in range(len(e)):
                    elements[e][r] = global_node_numbering[elements[e][r]]
        draw_sparsity_pattern(elements, num_elements)
