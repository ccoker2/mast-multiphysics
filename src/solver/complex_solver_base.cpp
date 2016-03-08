/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


// MAST includes
#include "solver/complex_solver_base.h"
#include "base/complex_assembly_base.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/dof_map.h"

// PETSc includes
#include <petscmat.h>


MAST::ComplexSolverBase::ComplexSolverBase():
_assembly(NULL),
tol(1.0e-3),
max_iters(20) {
    
}



MAST::ComplexSolverBase::~ComplexSolverBase() {
    
}




void
MAST::ComplexSolverBase::set_assembly(MAST::ComplexAssemblyBase& assembly) {
    
    _assembly = &assembly;
}



void
MAST::ComplexSolverBase::clear_assembly() {
    
    _assembly = NULL;
}




libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::real_solution() {
    
    libMesh::System& sys = _assembly->system();
    
    std::string nm;
    nm = sys.name();
    nm += "real_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);

    return sys.get_vector(nm);
}



const libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::real_solution() const {
    
    libMesh::System& sys = _assembly->system();

    std::string nm;
    nm = sys.name();
    nm += "real_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);
    
    return sys.get_vector(nm);
}


libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::imag_solution() {
    
    libMesh::System& sys = _assembly->system();

    std::string nm;
    nm = sys.name();
    nm += "imag_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);
    
    return sys.get_vector(nm);
}


const libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::imag_solution() const {
    
    libMesh::System& sys = _assembly->system();

    std::string nm;
    nm = sys.name();
    nm += "imag_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);
    
    return sys.get_vector(nm);
}



void
MAST::ComplexSolverBase::solve() {
    
    
    //  The complex system of equations
    //     (J_R + i J_I) (x_R + i x_I) + (r_R + i r_I) = 0
    //  is rewritten as
    //     [ J_R   -J_I] {x_R}  +  {r_R}  = {0}
    //     [ J_I    J_R] {x_I}  +  {r_I}  = {0}
    //
    
    
    // continue iterations till the L2 residual of both real and imaginary
    // parts is satisfied
    
    bool
    if_cont = true,
    if_re   = true;  // keeps track of whether real or imag part is being solved
    
    Real
    res_l2  = 0.;
    
    unsigned int
    iters   = 0;
    
    libMesh::System& sys = _assembly->system();

    
    while (if_cont) {
        
        // tell the assembly object to now assemble
        if (if_re) {

            // swap the solution with the real part
            *(sys.solution) =  this->real_solution();
            _assembly->set_assemble_real_part();
            libMesh::out << "Solving Real Part: " << std::endl;
        }
        else {
            
            *(sys.solution) =  this->imag_solution();
            _assembly->set_assemble_imag_part();
            libMesh::out << "Solving Imaginary Part: " << std::endl;
        }
        sys.solution->close();
        
        // tell implicit system to solve
        sys.solve();
        
        // now copy the solution vector back
        if (if_re) {
            
            this->real_solution() = (*sys.solution);
            this->real_solution().close();
            if_re = false;
        }
        else {
            
            this->imag_solution() = (*sys.solution);
            this->imag_solution().close();
            if_re = true;
        }
        
        // check the residual
        res_l2  =  _assembly->residual_l2_norm();
        libMesh::out
        << "Iter: " << iters
        << "   Complex Residual L2-norm: " << res_l2 << std::endl;
        
        
        // increment the iteration counter
        iters++;
        
        
        if (res_l2 <= tol ||
            iters  >= max_iters) {
            
            libMesh::out
            << "Terminating complex solver iterations!" << std::endl;
            if_cont = false;
        }
    }
    
}




void
MAST::ComplexSolverBase::solve_pc_fieldsplit() {
    
    START_LOG("complex_solve()", "PetscFieldSplitSolver");

    // get reference to the system
    MAST::NonlinearSystem& sys =
    dynamic_cast<MAST::NonlinearSystem&>(_assembly->system());
    
    // create a petsc nested matrix of size 2x2
    PetscErrorCode   ierr;
    Mat              mat;
    Vec              res, sol, res_vec_R, res_vec_I, sol_vec_R, sol_vec_I;
    std::vector<IS>  is(2);
    std::vector<Mat> sub_mats(4);
    sub_mats[0] = dynamic_cast<libMesh::PetscMatrix<Real>*>(sys.matrix)->mat();   // real
    sub_mats[1] = dynamic_cast<libMesh::PetscMatrix<Real>*>(sys.matrix_A)->mat(); // -imag
    sub_mats[2] = dynamic_cast<libMesh::PetscMatrix<Real>*>(sys.matrix_B)->mat(); // imag
    sub_mats[3] = dynamic_cast<libMesh::PetscMatrix<Real>*>(sys.matrix)->mat();   // real
    
    ierr = MatCreateNest(_assembly->system().comm().get(),
                         2, NULL,
                         2, NULL,
                         &sub_mats[0],
                         &mat);
    CHKERRABORT(sys.comm().get(), ierr);


    // get the IS belonging to each block
    ierr  =    MatNestGetISs(mat, &is[0], NULL); CHKERRABORT(sys.comm().get(), ierr);
    
    
    
    // setup vector for residual
    ierr = VecCreate(sys.comm().get(), &res);                      CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecSetSizes(res, PETSC_DECIDE, 2*sys.solution->size()); CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecSetType(res, VECMPI);                                CHKERRABORT(sys.comm().get(), ierr);
    
    
    // setup vector for solution
    ierr  =   VecDuplicate(res, &sol);                             CHKERRABORT(sys.comm().get(), ierr);

    
    
    // get the residual subvectors
    ierr = VecGetSubVector(res, is[0], &res_vec_R);                CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecGetSubVector(res, is[1], &res_vec_I);                CHKERRABORT(sys.comm().get(), ierr);

    // get the solution subvectors
    ierr = VecGetSubVector(sol, is[0], &sol_vec_R);                CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecGetSubVector(sol, is[1], &sol_vec_I);                CHKERRABORT(sys.comm().get(), ierr);
    
    
    // clone vectors for use as system RHS
    std::auto_ptr<libMesh::NumericVector<Real> >
    res_R(new libMesh::PetscVector<Real>(res_vec_R, sys.comm())),
    res_I(new libMesh::PetscVector<Real>(res_vec_I, sys.comm())),
    sol_R(new libMesh::PetscVector<Real>(sol_vec_R, sys.comm())),
    sol_I(new libMesh::PetscVector<Real>(sol_vec_I, sys.comm()));
    
    sol_R->zero();
    sol_R->close();
    sol_I->zero();
    sol_I->close();
    
    // assemble the matrices
    _assembly->residual_and_jacobian_field_split(*sol_R,
                                                 *sol_I,
                                                 *res_R,
                                                 *res_I,
                                                 *sys.matrix,    // J_R
                                                 *sys.matrix_B,  // J_I
                                                 sys);

    // restore the subvectors
    //ierr = VecRestoreSubVector(res, is[0], &res_vec_R); CHKERRABORT(sys.comm().get(), ierr);
    //ierr = VecRestoreSubVector(res, is[1], &res_vec_I); CHKERRABORT(sys.comm().get(), ierr);
    
    // copy the imag Jacobian for the first row of the sys of eq.
    sys.matrix_A->zero();
    sys.matrix_A->close();
    sys.matrix_A->add(-1., *sys.matrix_B);
    sys.matrix_A->close();
    
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);  CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);    CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecAssemblyBegin(res);                      CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecAssemblyEnd(res);                        CHKERRABORT(sys.comm().get(), ierr);
    
    
    
    // now initialize the KSP and ask for solution.
    KSP        ksp;
    PC         pc;
    
    // setup the KSP
    ierr = KSPCreate(sys.comm().get(), &ksp); CHKERRABORT(sys.comm().get(), ierr);
    ierr = KSPSetOperators(ksp, mat, mat);    CHKERRABORT(sys.comm().get(), ierr);
    ierr = KSPSetFromOptions(ksp);            CHKERRABORT(sys.comm().get(), ierr);

    // setup the PC
    ierr = KSPGetPC(ksp, &pc);                CHKERRABORT(sys.comm().get(), ierr);
    ierr = PCFieldSplitSetIS(pc, NULL, is[0]);CHKERRABORT(sys.comm().get(), ierr);
    ierr = PCFieldSplitSetIS(pc, NULL, is[1]);CHKERRABORT(sys.comm().get(), ierr);
    ierr = PCSetFromOptions(pc);              CHKERRABORT(sys.comm().get(), ierr);
    
    
    
    // now solve
    ierr = KSPSolve(ksp, res, sol);
    
    
    // assemble the matrices
    _assembly->residual_and_jacobian_field_split(*sol_R,
                                                 *sol_I,
                                                 *res_R,
                                                 *res_I,
                                                 *sys.matrix,    // J_R
                                                 *sys.matrix_B,  // J_I
                                                 sys);

    // copy solutions for output
    /*this->real_solution() = *sol_R;
    this->imag_solution() = *sol_I;
    this->real_solution().close();
    this->imag_solution().close();*/
    
    
    // restore the subvectors
    ierr = VecRestoreSubVector(res, is[0], &res_vec_R); CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecRestoreSubVector(res, is[1], &res_vec_I); CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecRestoreSubVector(sol, is[0], &sol_vec_R); CHKERRABORT(sys.comm().get(), ierr);
    ierr = VecRestoreSubVector(sol, is[1], &sol_vec_I); CHKERRABORT(sys.comm().get(), ierr);

    
    // destroy the objects
    ierr = KSPDestroy(&ksp);
    ierr = MatDestroy(&mat);
    ierr = VecDestroy(&res);
    ierr = VecDestroy(&sol);
    
    STOP_LOG("complex_solve()", "PetscFieldSplitSolver");
}



void
MAST::ComplexSolverBase::solve_block_matrix()  {
    
    START_LOG("complex_solve()", "PetscComplexSolve");

    // get reference to the system
    MAST::NonlinearSystem& sys =
    dynamic_cast<MAST::NonlinearSystem&>(_assembly->system());

    libMesh::DofMap& dof_map = sys.get_dof_map();
    
    const PetscInt
    my_m = dof_map.n_dofs(),
    my_n = my_m,
    n_l  = dof_map.n_dofs_on_processor(sys.processor_id()),
    m_l  = n_l;

    const std::vector<libMesh::dof_id_type>
    & n_nz = dof_map.get_n_nz(),
    & n_oz = dof_map.get_n_oz();

    
    // create the matrix
    PetscErrorCode   ierr;
    Mat              mat;

    ierr = MatCreate(sys.comm().get(), &mat);                   CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatSetSizes(mat, 2*m_l, 2*n_l, 2*my_m, 2*my_n);      CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatSetType(mat, MATBAIJ);                            CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatSetBlockSize(mat, 2);                             CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatSeqBAIJSetPreallocation (mat, 2,
                                       0, (PetscInt*)&n_nz[0]); CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatMPIBAIJSetPreallocation (mat, 2,
                                       0, (PetscInt*)&n_nz[0],
                                       0, (PetscInt*)&n_oz[0]); CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatSetOption(mat,
                        MAT_NEW_NONZERO_ALLOCATION_ERR,
                        PETSC_TRUE);                            CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatSetFromOptions(mat);                              CHKERRABORT(sys.comm().get(), ierr);

    
    // now create the vectors
    Vec              res_vec, sol_vec;
    
    ierr = MatCreateVecs(mat, &res_vec, PETSC_NULL);            CHKERRABORT(sys.comm().get(), ierr);
    ierr = MatCreateVecs(mat, &sol_vec, PETSC_NULL);            CHKERRABORT(sys.comm().get(), ierr);

    
    std::auto_ptr<libMesh::SparseMatrix<Real> >
    jac_mat(new libMesh::PetscMatrix<Real>(mat, sys.comm()));

    std::auto_ptr<libMesh::NumericVector<Real> >
    res(new libMesh::PetscVector<Real>(res_vec, sys.comm())),
    sol(new libMesh::PetscVector<Real>(sol_vec, sys.comm()));
    
    // assemble the matrix
    _assembly->residual_and_jacobian_blocked(*sol,
                                             *res,
                                             *jac_mat,
                                             sys);
    

    // now initialize the KSP and ask for solution.
    KSP        ksp;
    PC         pc;
    
    // setup the KSP
    ierr = KSPCreate(sys.comm().get(), &ksp); CHKERRABORT(sys.comm().get(), ierr);
    ierr = KSPSetOperators(ksp, mat, mat);    CHKERRABORT(sys.comm().get(), ierr);
    ierr = KSPSetFromOptions(ksp);            CHKERRABORT(sys.comm().get(), ierr);
    
    // setup the PC
    ierr = KSPGetPC(ksp, &pc);                CHKERRABORT(sys.comm().get(), ierr);
    ierr = PCSetFromOptions(pc);              CHKERRABORT(sys.comm().get(), ierr);
    
    
    
    // now solve
    ierr = KSPSolve(ksp, res_vec, sol_vec);

    // evaluate the residual again
    _assembly->residual_and_jacobian_blocked(*sol,
                                             *res,
                                             *jac_mat,
                                             sys);
    
    
    // copy the solution to separate real and imaginary vectors
    libMesh::NumericVector<Real>
    &sol_R = this->real_solution(),
    &sol_I = this->imag_solution();
    
    unsigned int
    first = sol_R.first_local_index(),
    last  = sol_I.last_local_index();
    
    for (unsigned int i=first; i<last; i++) {
        sol_R.set(i, (*sol)(  2*i));
        sol_I.set(i, (*sol)(2*i+1));
    }
    
    sol_R.close();
    sol_I.close();
    
    STOP_LOG("complex_solve()", "PetscComplexSolve");
}




