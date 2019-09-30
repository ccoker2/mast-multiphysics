/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "base/transient_assembly_elem_operations.h"
#include "base/elem_base.h"


namespace MAST {
    bool
    transient_is_numerical_zero(const Real v, const Real eps) {
        return fabs(v) <= eps;
    }


    bool
    transient_compare(const Real v1, const Real v2, const Real tol) {

        const Real
                eps      = 1.0e-7;

        bool rval = false;

        // check to see if the values are both small enough
        // to be zero
        if (MAST::transient_is_numerical_zero(v1, eps) &&
            MAST::transient_is_numerical_zero(v2, eps))
            rval = true;
            // check to see if the absolute difference is small enough
        else if (MAST::transient_is_numerical_zero(v1-v2, eps))
            rval = true;
            // check to see if the relative difference is small enough
        else if (fabs(v1) > 0)
            rval = fabs((v1-v2)/v1) <= tol;

        return rval;
    }

    bool
    transient_compare_matrix(const RealMatrixX& m0, const RealMatrixX& m, const Real tol) {

        unsigned int
                m0_rows = (unsigned int) m0.rows(),
                m0_cols = (unsigned int) m0.cols(),
                n_passes,
                k=0;
        Real
                max_err,
                min_err,
                mean_err;
        RealVectorX
                errs   = RealVectorX::Zero(m.rows()*m.cols()),
                passes = RealVectorX::Zero(m.rows()*m.cols());


        libmesh_assert_equal_to(m0_rows,  m.rows());
        libmesh_assert_equal_to(m0_cols,  m.cols());


        bool pass = true;
        for (unsigned int i=0; i<m0_rows; i++) {
            for (unsigned int j=0; j<m0_cols; j++) {
                if (!MAST::transient_compare(m0(i, j), m(i, j), tol)) {
                    libMesh::out << "Failed comparison at (i,j) = ("
                         << i << ", " << j << ") : "
                         << "expected: " << m0(i, j) << "  , "
                         << "computed: " << m(i, j) << " : "
                         << "diff: " << m0(i, j) - m(i, j) << " , "
                         << "tol: " << tol << std::endl;
                    pass = false;
                }
                errs(k) = m0(i, j) - m(i, j);
                passes(k) = pass;
                k++;
            }
        }
        max_err = errs.cwiseAbs().maxCoeff();
        min_err = errs.cwiseAbs().minCoeff();
        mean_err = errs.mean();
        n_passes = (passes.array() != 0).count();

        libMesh::out << "Max error: " << max_err << ", "
        << " Min error: " << min_err
        << " Mean error: " << mean_err
        << " Number of elements within tolerance " << n_passes << std::endl;

        return pass;
    }
}



MAST::TransientAssemblyElemOperations::TransientAssemblyElemOperations():
MAST::AssemblyElemOperations() {
    
}



MAST::TransientAssemblyElemOperations::~TransientAssemblyElemOperations() {
    
}


void
MAST::TransientAssemblyElemOperations
::check_element_numerical_jacobian()  {

    RealVectorX
            sol,
            vel;

    sol = this->_physics_elem->MAST::ElementBase::sol(0);
    vel = this->_physics_elem->MAST::ElementBase::vel(0);
    unsigned int ndofs = (unsigned int)sol.size();

    // initial
    RealVectorX
            f_m = RealVectorX::Zero(ndofs),
            f_x = RealVectorX::Zero(ndofs);

    // perturbed
    RealVectorX
            df_m = RealVectorX::Zero(ndofs),
            df_x = RealVectorX::Zero(ndofs),
            dsol = RealVectorX::Zero(ndofs),
            dvel = RealVectorX::Zero(ndofs);

    // analytical
    RealMatrixX
            f_m_jac_xdot0 = RealMatrixX::Zero(ndofs,ndofs),
            f_m_jac0 = RealMatrixX::Zero(ndofs,ndofs),
            f_x_jac0 = RealMatrixX::Zero(ndofs,ndofs);

    // perturbed
    RealMatrixX
            f_m_jac_xdot = RealMatrixX::Zero(ndofs,ndofs),
            f_m_jac = RealMatrixX::Zero(ndofs,ndofs),
            f_x_jac = RealMatrixX::Zero(ndofs,ndofs);

    // dummy vars
    RealMatrixX
            dummy_mat1 = RealMatrixX::Zero(ndofs,ndofs),
            dummy_mat2 = RealMatrixX::Zero(ndofs,ndofs),
            dummy_mat3 = RealMatrixX::Zero(ndofs,ndofs);

    RealVectorX
            dummy_vec1 = RealVectorX::Zero(ndofs);

    // numerical
    RealMatrixX
            jac = RealMatrixX::Zero(ndofs,ndofs),
            mas = RealMatrixX::Zero(ndofs,ndofs);



    // initial
    this->set_elem_solution(sol);
    this->elem_calculations(true, f_m, f_x, f_m_jac_xdot0, f_m_jac0, f_x_jac0);
    Real delta = 1.0e-8;

    // numerical jacobian matrix
    for (unsigned int i=0; i<sol.size(); i++) {
        dsol = sol;
        dsol(i) += delta;
        this->set_elem_solution(dsol);

        this->elem_calculations(false, dummy_vec1, df_x, dummy_mat1, dummy_mat2, dummy_mat3);
        jac.col(i) = (df_x-f_x)/delta;
    }

    this->set_elem_solution(sol);
    // calculate numerical mass matrix
    for (unsigned int i=0; i<vel.size(); i++) {
        dvel = vel;
        dvel(i) += delta;
        this->set_elem_velocity(dvel);

        this->elem_calculations(false, df_m, dummy_vec1, dummy_mat1, dummy_mat2, dummy_mat3);
        mas.col(i) = (df_m-f_m)/delta;
    }


    // write the numerical and analytical jacobians
//    libMesh::out
//            << "Analytical Jacobian: " << std::endl
//            << f_x_jac0
//            << std::endl << std::endl
//            << "Numerical Jacobian: " << std::endl
//            << jac
//            << std::endl << std::endl;

    MAST::transient_compare_matrix(jac, f_x_jac0, 1.0e1);
    // set the original solution vector for the element

//    // write the numerical and analytical mass matrices
//    libMesh::out
//            << "Analytical Mass Matrix: " << std::endl
//            << f_m_jac_xdot0
//            << std::endl << std::endl
//            << "Numerical Mass Matrix: " << std::endl
//            << mas
//            << std::endl << std::endl;

    MAST::transient_compare_matrix(mas, f_m_jac_xdot0, 1.0e-1);
    // set the original solution vector for the element
    this->set_elem_solution(sol);
}
