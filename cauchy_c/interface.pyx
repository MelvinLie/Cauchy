# distutils: sources = [cauchy_c/c-algorithms/assembly_tools.c, cauchy_c/c-algorithms/gauss_legendre.c]
# distutils: include_dirs = cauchy_c/c-algorithms/

import numpy as np
cimport numpy as np
cimport assembly_tools
from cython cimport view
from libc.stdlib cimport malloc, free

def assemble_V(mesh):

    print("Assemble Discrete Single Layer Operator.")

    N = mesh.shape[0]

    # close the mesh (append first node to the end) this is needed for some c
    # functions
    mesh = np.append(mesh, mesh[0:1, :])

    mesh_f = mesh.flatten()

    cdef np.ndarray[double, ndim=1, mode = 'c'] np_buff = np.ascontiguousarray(mesh_f, dtype = np.double)
    cdef double* mesh_p = <double*> np_buff.data

    cdef assembly_tools.disc* gamma = assembly_tools.make_parameterization(mesh_p, N, 1.0)

    cdef double* V_c = <double *> malloc( N * N * sizeof(double))
    V_c = assembly_tools.assemble_V_c(gamma)

    # convert to numpy array
    cdef view.array V_array = <double[:N*N]> V_c

    V = np.asarray(V_array)
    V.shape = (N, N)

    return V

def assemble_K(mesh):

    print("Assemble Discrete Double Layer Operator.")

    N = mesh.shape[0]

    # close the mesh (append first node to the end) this is needed for some c
    # functions
    mesh = np.append(mesh, mesh[0:1, :])

    mesh_f = mesh.flatten()


    cdef np.ndarray[double, ndim=1, mode = 'c'] np_buff = np.ascontiguousarray(mesh_f, dtype = np.double)
    cdef double* mesh_p = <double*> np_buff.data

    cdef assembly_tools.disc* gamma = assembly_tools.make_parameterization(mesh_p, N, 1.0)

    cdef double* K_c = <double *> malloc( N * N * sizeof(double))
    K_c = assembly_tools.assemble_K_c(gamma)

    # convert to numpy array
    cdef view.array K_array = <double[:N*N]> K_c

    K = np.asarray(K_array)
    K.shape = (N, N)

    return K

def assemble_D(mesh):

    print("Assemble Discrete Hypersingular Operator.")

    N = mesh.shape[0]

    # close the mesh (append first node to the end) this is needed for some c
    # functions
    mesh = np.append(mesh, mesh[0:1, :])

    mesh_f = mesh.flatten()


    cdef np.ndarray[double, ndim=1, mode = 'c'] np_buff = np.ascontiguousarray(mesh_f, dtype = np.double)
    cdef double* mesh_p = <double*> np_buff.data

    cdef assembly_tools.disc* gamma = assembly_tools.make_parameterization(mesh_p, N, 1.0)

    cdef double* D_c = <double *> malloc( N * N * sizeof(double))
    D_c = assembly_tools.assemble_D_c(gamma, NULL)

    # convert to numpy array
    cdef view.array D_array = <double[:N*N]> D_c

    D = np.asarray(D_array)
    D.shape = (N, N)

    return D


def assemble_M(mesh):

    print("Assemble Discrete Mass Matrix.")

    N = mesh.shape[0]

    # close the mesh (append first node to the end) this is needed for some c
    # functions
    mesh = np.append(mesh, mesh[0:1, :])

    mesh_f = mesh.flatten()


    cdef np.ndarray[double, ndim=1, mode = 'c'] np_buff = np.ascontiguousarray(mesh_f, dtype = np.double)
    cdef double* mesh_p = <double*> np_buff.data

    cdef assembly_tools.disc* gamma = assembly_tools.make_parameterization(mesh_p, N, 1.0)

    cdef unsigned long long *rc  = <unsigned long long *> malloc( 4*N * sizeof(unsigned long long))
    cdef double *M_c = <double *> malloc( 2 * N * sizeof(double))

    assembly_tools.assemble_M_c(gamma, rc, M_c)

    # convert to numpy array
    M = np.zeros((N, N))

    for i in range(2*N):
      M[rc[2*i], rc[2*i+1]] = M_c[i]

    return M


def assemble_S(mesh):

  # the number of nodes
  N = mesh.shape[0]

  # compute the segment lengths
  l = np.linalg.norm(mesh, axis=1)

  # make the a vector
  a = 0.5*l
  a[1:] += l[:-1]
  a[0] += l[-1]

  return a @ a.T

def grad_representation_formula(mesh, u, g, points):
    '''Evaluate the gradient of the representation formula
    for given field points.

    :param mesh:
        The boundary mesh.

    :param u:
        The Dirichlet data.

    :param g:
        The Neumann data.
    
    :param points:
        The evaluation points.

    :return:
        The gradient field.
    '''

    # print("Evaluate the gradient of the representation formula.")

    # number of nodes
    N = mesh.shape[0]

    # the number of observation points
    M = points.shape[0]

    # close the mesh (append first node to the end) this is needed for some c
    # functions
    mesh = np.append(mesh, mesh[0:1, :])

    # ===================
    # Convert python -> C
    # ===================
    cdef np.ndarray[double, ndim=1, mode = 'c'] mesh_buff = np.ascontiguousarray(mesh.flatten(), dtype = np.double)
    cdef double* mesh_c = <double*> mesh_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] points_buff = np.ascontiguousarray(points.flatten(), dtype = np.double)
    cdef double* points_c = <double*> points_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] u_buff = np.ascontiguousarray(u.flatten(), dtype = np.double)
    cdef double* u_c = <double*> u_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] g_buff = np.ascontiguousarray(g.flatten(), dtype = np.double)
    cdef double* g_c = <double*> g_buff.data

    cdef assembly_tools.disc* gamma = assembly_tools.make_parameterization(mesh_c, N, 1.0)


    # ===================
    # Launch C code
    # ===================
    cdef double *grad_c = assembly_tools.evaluate_grad_u(points_c, M, gamma, u_c, g_c)

    # ===================
    # Convert C -> python
    # ===================

    # convert to numpy array
    cdef view.array grad_array = <double[:2*M]> grad_c

    grad = np.asarray(grad_array)
    grad.shape = (M, 2)

    return grad