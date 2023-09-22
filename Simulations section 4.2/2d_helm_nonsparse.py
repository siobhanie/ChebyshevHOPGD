
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 10:14:16 2022

@author: siobhanie
"""
from __future__ import print_function
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import matplotlib as p
import scipy.io
from scipy.io import savemat
import scipy.io
from numpy import * 
from mshr import * 


parameters ['linear_algebra_backend'] = 'Eigen'
parameters ['reorder_dofs_serial'] = False


def helmholtz(f,f1,mesh,mu1,mu2): 
    
    #Set up 
    V = FunctionSpace(mesh, 'Lagrange', 1)
    u, v = TrialFunction(V), TestFunction(V)
    
    #Boundart condition (Dirichlet)
    u_D=Constant(0)
    #u_D = Expression('sin(pi*x[0])*sin(pi*x[1])',degree=2,domain=mesh) #rhs
    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, u_D, boundary)
    
    #Variational form
    L = f*v*dx
    
    a0 = (-(dot(grad(u), grad(v))))*dx
    A0, b = assemble_system(a0,L,bc);
    rows,cols,vals = as_backend_type(A0).data() 
    A0 = as_backend_type(A0).sparray()

    a2 = (u*v)*dx
    A2, b = assemble_system(a2,L,bc);
    rows,cols,vals = as_backend_type(A2).data() 
    A2 = as_backend_type(A2).sparray()
    
    a3 = f1*(u*v)*dx
    A3, b = assemble_system(a3,L,bc);
    rows,cols,vals = as_backend_type(A3).data() 
    A3 = as_backend_type(A3).sparray()
    

    a = a0 + (cos(mu1)+mu1**3)*a2 + (sin(mu2) + mu2**2)*a3
    u = Function(V)
    solve(a == L,u,bc)
    
    bvec = zeros(len(b),)
    for i in range(len(b)):
        bvec[i] = (b[i])
    b = bvec
    
    uvec = zeros(len(b),)
    for i in range(len(b)):
        uvec[i] = (u.vector()[i])
    u1 = uvec    
    print(len(b))

    normi = p.colors.Normalize(vmin=-.025, vmax=0)
    plt.jet()
    plt.margins(0,0)
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    c = plot(u)
    c.set_norm(normi)
    plt.colorbar(c)
    plt.show()


    return 

 
if __name__ == "__main__":
    
    N1 = 5
    N2 = 5
    mu1=2
    mu2=2
    #n = 270 #grid pts 
    n=270
    
    D1 = Rectangle(Point(0,0),Point(1,1))
    D2 = Circle(Point(0.2,0.2),.04)
    D3 = Circle(Point(0.8,0.8),.04)
    D4 = Circle(Point(0.5,0.5),.05) 
    domain = D1-D4
    mesh = generate_mesh(domain,n)
    
    f = Expression('sin(pi*x[0])*sin(pi*x[1])',degree=2,domain=mesh) #rhs
    f1 = Expression('x[1]',degree=2,domain=mesh)

    helmholtz(f,f1,mesh,mu1,mu2)

