//
// Created by Christian Coker on 10/22/19.
//
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

// C/C++ includes.
#include <iostream>

// MAST includes.
#include "examples/fluid/meshing/cylinder.h"
#include "examples/fluid/meshing/naca0012.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "examples/fluid/meshing/panel_mesh_3D.h"
#include "examples/fluid/meshing/naca0012_wing.h"
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "fluid/flight_condition.h"
#include "fluid/integrated_force_output.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "solver/stabilized_first_order_transient_sensitivity_solver.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_elem_type.h"    // ElemType
#include "libmesh/fe_type.h"           // FEFamily, Order
#include "libmesh/parallel_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/periodic_boundary.h"


// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/libmesh.h"


// Test includes
#include "tests/fluid/base/fluid_elem_initialization.h"
#include "tests/base/test_comparisons.h"


libMesh::LibMeshInit     *_libmesh_init         = nullptr;
const Real                _frac                 = 1.e-4;
const Real                _delta                = 1.e-4;
const Real                _tol                  = 1.e-5;


//int main(int argc, const char** argv)
//{
//    // create the libMeshInit function
//    _libmesh_init =
//            new libMesh::LibMeshInit(argc,
//                                     argv);
//
//    BuildFluidElem elem;
//    bool if_viscous = true;
//    elem.init(if_viscous);
//
//    RealVectorX sol;
//    elem.init_solution_for_elem(sol);
//
//    MAST::FluidElemBase::check_element_diffusion_flux_jacobian()
//
//    delete _libmesh_init;
//}


int main(int argc, const char** argv)
{

    // create the libMeshInit function
    _libmesh_init =
    new libMesh::LibMeshInit(argc,
                         argv);


    BuildFluidElem elem;
    bool if_viscous = true;

    elem.init(if_viscous);

    elem.check_element_diffusion_flux_jacobian();

    delete _libmesh_init;
}
