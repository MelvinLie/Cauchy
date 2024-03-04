from .cauchy_c import assemble_K
from .cauchy_c import assemble_D
from .cauchy_c import assemble_M
from .cauchy_c import assemble_S
from .cauchy_c import grad_representation_formula

import numpy as np

def neumann_2_dirichlet(mesh, g, g_std):
    '''Compute the dirichlet data for given
    (measured) Neumann data.

    :param mesh:
        The boundary mesh.

    :param g:
        The Neumann data.

    :param g_std:
        The standard deviation in the measured
        Neumann data.

    :return:
        The dirichlet data, and the posterior
        covariance matrix.
        Moreover, return also the solution matrix A.
    '''

    # compute all integral operators
    K = assemble_K(mesh)
    D = assemble_D(mesh)
    M = assemble_M(mesh)
    S = assemble_S(mesh)

    # compute the right hand side
    b = (0.5 * M.T - K.T) @ g

    # compute the Neumann data
    u = np.linalg.solve(D + S, b)

    # compute the Neumann 2 Dirichlet map
    A = np.linalg.solve( D + S, 0.5 * M.T - K.T )

    # compute the covariance matrix for the Dirichlet data
    Sigma_u = A * np.diag(g_std) * A.T

    return u, Sigma_u, A

def compute_H_matrix_neumann(mesh, points, A):
    '''Compute the obserbation matrix for the
    field evaluation (grad) at some points,
    for the Neumann 2 Dirichlet map.

    :param mesh:
        The boundary mesh.

    :param points:
        The points.

    :param A:
        The matrix for the Neumann 2 Dirichlet map.

    :return:
        The observation matrices for the two
        spatial derivatives.
    '''

    # number of nodes
    num_nodes = mesh.shape[0]

    # number of points
    num_points = points.shape[0]

    # allocate the matrices
    Hx = np.zeros((num_points, num_nodes))
    Hy = np.zeros((num_points, num_nodes))

    # make a neumann data vector
    g = np.zeros((num_nodes, ))

    # fill it
    for i in range(num_nodes):

        # set this basis function 1
        g[i] = 1.0

        # compute the field
        grad = grad_representation_formula(mesh, A[:, i], g, points)

        # set this basis function 1
        g[i] = 0.0

        # fill the matrix
        Hx[:, i] = grad[:, 0]
        Hy[:, i] = grad[:, 1]

    return Hx, Hy

def compute_variance(H_x, H_y, Sigma):
    '''Compute the standard deviation of the field evaluation for
    given observation operators and input data.

    :param Hx:
        The x gradient evaluation matrix.

    :param Hy:
        The y gradient evaluation matrix.

    :param Sigma:
        The covariance matrix.

    :return:
        The variances, i.e. the diagonals of the
        posterior covariance matrices for grad_x
        and grad_y.
    '''
    # the number of observations
    num_eval = H_x.shape[0]

    # the right hand side matrices
    R_x = Sigma @ H_x.T
    R_y = Sigma @ H_y.T

    # make the return vector
    var_x = np.zeros((num_eval, ))
    var_y = np.zeros((num_eval, ))

    # loop over all evaluations
    for i in range(num_eval):
        var_x[i] = np.sum(H_x[i, :] * R_x[:, i])
        var_y[i] = np.sum(H_y[i, :] * R_y[:, i])

    return var_x, var_y


def neumann_BEM_uq(mesh, A, g, g_std, points):
    '''Launch the standard Neumann BEM analysis and
    uncertainty quantification.
    
    :param mesh:
        The boundary mesh.

    :param A:
        The matrix for the neumann to dirichlet map.

    :param g:
        The (measured) Neumann data.

    :param g_std:
        The standard deviation in the Neumann data.

    :param points:
        The evaluation points.

    :return:
        The gradient field as well as the standard
        deviations in each evaluation point.
        Moreover, we return the solution matrix A.
    '''

    # compute the observation matrix
    H_x, H_y = compute_H_matrix_neumann(mesh, points, A)

    # compute the means
    grad_x = H_x @ g
    grad_y = H_y @ g

    # compute the variances
    var_x, var_y = compute_variance(H_x, H_y, np.diag(g_std**2))

    # assemble (M x 2) arrays
    grad = np.zeros((grad_x.shape[0], 2))
    grad[:, 0] = grad_x
    grad[:, 1] = grad_y

    var = np.zeros((var_x.shape[0], 2))
    var[:, 0] = var_x
    var[:, 1] = var_y

    return grad, var 