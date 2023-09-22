#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 15:16:23 2022

@author: siobhanie
"""
from __future__ import print_function
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import scipy.io
from scipy.io import savemat
from numpy import * 
from mshr import * 


parameters ['linear_algebra_backend'] = 'Eigen'
parameters ['reorder_dofs_serial'] = False


def helmholtz(f,f1,mesh,mu): 
    
    #Set up 
    V = FunctionSpace(mesh, 'Lagrange', 1)
    u, v = TrialFunction(V), TestFunction(V)
    
    #Boundart condition (Dirichlet)
    u_D = Constant(0)
    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, u_D, boundary)
    
    #Variational form
    L = f*v*dx
    
    a0 = (-(dot(grad(u), grad(v))))*dx
    A0, b = assemble_system(a0,L,bc);
    rows,cols,vals = as_backend_type(A0).data() 
    A0 = as_backend_type(A0).sparray()

    a1 = (u*v)*dx
    A1, b = assemble_system(a1,L,bc);
    rows,cols,vals = as_backend_type(A1).data() 
    A1 = as_backend_type(A1).sparray()
    
    a2 = f1*(u*v)*dx
    A2, b = assemble_system(a2,L,bc);
    rows,cols,vals = as_backend_type(A2).data() 
    A2 = as_backend_type(A2).sparray()
    
    a = a0 + (sin(mu[0]))**2*a1 + (cos(mu[1]))**2*a2
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
    #print(u.vector()[1])

    savemat('param_est.mat', {'A0':A0,'A1':A1,'A2':A2,'b':b})

    # Plot solution
    plt.margins(0,0)
        
    c=plot(u)
    plt.jet()
    plt.colorbar(c)
    plt.margins(0,0)
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title('$\mu=3$',y=-.25)
    
    return 

 
if __name__ == "__main__":
    

    mu= array([.9,.1])
    n = 80 #grid pts 1
    
    D2 = Circle(Point(0.2,0.2),.04)
    D3 = Circle(Point(0.8,0.8),.04)
    D4 = Circle(Point(0.5,0.5),.08) 
    D1 = Rectangle(Point(0,0),Point(1.0,1.0))
    domain=D1 #-D4
    mesh = generate_mesh(domain,n)
    
    f = Expression('exp(-x[0]*x[1])',degree=2,domain=mesh) #rhs
    f1 = Expression('x[0]',degree=2,domain=mesh)

    helmholtz(f,f1,mesh,mu)

