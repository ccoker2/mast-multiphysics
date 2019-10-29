//
// Created by coker on 10/23/19.
//
// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/libmesh.h"

// Test includes
#include "fluid/base/fluid_elem_initialization.h"
#include "base/test_comparisons.h"


libMesh::LibMeshInit     *_libmesh_init         = nullptr;
const Real                _frac                 = 1.e-4;
const Real                _delta                = 1.e-4;
const Real                _tol                  = 1e-8;


int main(int argc, const char** argv) {

    // create the libMeshInit function
    _libmesh_init = new libMesh::LibMeshInit(argc,
                                             argv);

    struct Check {
            bool            jac_xdot;
            Real            frac;
            Real            delta;
            Real            tol;
            BuildFluidElem e;
            void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
                e._fluid_elem->internal_residual(jac, f, j);
            }
    };


    Check val;
    val.jac_xdot = false;
    val.frac     = _frac;
    val.delta    = _delta;
    // a smaller tolerance is required for the internal resisudal since
    // an exact Jacobian is not computed for the stabilization terms
//     val.tol      = 1.e-2;
    val.tol      = _tol;


//    BuildFluidElem e;
//    e.init(if_viscous);

    bool if_viscous = true;
    val.e.init(if_viscous);
    val.e._delta = 0.1;
    val.e.check_jacobian(val);

//    this->_delta = 1.e-4;
//    this->init(true);

//    check_jacobian(val);

    delete _libmesh_init;
}