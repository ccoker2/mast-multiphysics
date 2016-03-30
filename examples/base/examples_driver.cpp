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
#include "examples/structural/bar_extension/bar_extension.h"
#include "examples/structural/beam_modal_analysis/beam_modal_analysis.h"
#include "examples/structural/beam_buckling_prestress/beam_column_buckling.h"
#include "examples/structural/beam_bending/beam_bending.h"
#include "examples/structural/beam_oscillating_load/beam_oscillating_load.h"
#include "examples/structural/beam_bending_with_offset/beam_bending_with_offset.h"
#include "examples/structural/beam_bending_thermal_stress_with_offset/beam_bending_thermal_stress.h"
#include "examples/structural/beam_optimization/beam_optimization.h"
#include "examples/structural/beam_optimization_single_stress_functional/beam_optimization.h"
#include "examples/structural/beam_optimization_section_offset/beam_optimization_section_offset.h"
#include "examples/structural/beam_optimization_thermal_stress/beam_optimization_thermal_stress.h"
#include "examples/structural/beam_piston_theory_flutter/beam_piston_theory_flutter.h"
#include "examples/structural/beam_piston_theory_time_accurate/beam_piston_theory_time_accurate.h"
#include "examples/structural/membrane_extension_uniaxial_stress/membrane_extension_uniaxial.h"
#include "examples/structural/membrane_extension_biaxial_stress/membrane_extension_biaxial.h"
#include "examples/structural/plate_bending/plate_bending.h"
#include "examples/structural/plate_oscillating_load/plate_oscillating_load.h"
#include "examples/structural/plate_bending_section_offset/plate_bending_section_offset.h"
#include "examples/structural/plate_bending_thermal_stress/plate_bending_thermal_stress.h"
#include "examples/structural/plate_modal_analysis/plate_modal_analysis.h"
#include "examples/structural/plate_buckling_prestress/plate_buckling_prestress.h"
#include "examples/structural/plate_piston_theory_flutter/plate_piston_theory_flutter.h"
#include "examples/structural/plate_thermally_stressed_piston_theory_flutter/plate_thermally_stressed_piston_theory_flutter.h"
#include "examples/structural/plate_optimization/plate_optimization.h"
#include "examples/structural/plate_optimization_single_stress_functional/plate_optimization_single_functional.h"
#include "examples/structural/plate_optimization_section_offset/plate_section_offset_optimization.h"
#include "examples/structural/plate_optimization_thermal_stress/plate_thermal_stress_optimization.h"
#include "examples/structural/stiffened_plate_optimization_thermal_stress/stiffened_plate_thermal_stress_optimization.h"
#include "examples/structural/stiffened_plate_optimization_piston_theory_flutter/stiffened_plate_piston_theory_flutter_optimization.h"
#include "examples/structural/topology_optim_2D/topology_optim_2D.h"
#include "optimization/npsol_optimization_interface.h"
#include "optimization/dot_optimization_interface.h"
#include "examples/fluid/panel_inviscid_analysis_2D/panel_inviscid_analysis_2d.h"
#include "examples/fluid/panel_inviscid_analysis_3D_half_domain/panel_inviscid_analysis_3D_half_domain.h"
#include "examples/fluid/panel_small_disturbance_frequency_domain_analysis_2D/panel_small_disturbance_frequency_domain_analysis_2d.h"
#include "examples/fluid/panel_small_disturbance_frequency_domain_3D/panel_small_disturbance_frequency_domain_inviscid_analysis_3D.h"
#include "examples/fluid/panel_small_disturbance_frequency_domain_3D_half_domain/panel_small_disturbance_frequency_domain_inviscid_analysis_3D_half_domain.h"
#include "examples/fsi/beam_flutter_solution/beam_euler_fsi_flutter_solution.h"
#include "examples/fsi/plate_flutter_solution/plate_euler_fsi_flutter_solution.h"
#include "examples/fsi/plate_flutter_solution_half_domain/plate_euler_fsi_half_domain_flutter_solution.h"
#include "examples/thermal/bar_transient/bar_transient.h"
#include "examples/thermal/bar_steady_state/bar_steady_state.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"


libMesh::LibMeshInit     *__init         = NULL;
MAST::FunctionEvaluation *__my_func_eval = NULL;




template <typename ValType>
void analysis(const std::string& case_name,
              bool with_sens,
              std::string& par_name)  {
    
    ValType run_case;
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, true);
        }
    }

}



template <typename ValType>
void eigenvalue_analysis(const std::string& case_name,
                         bool with_sens,
                         std::string& par_name)  {
    
    ValType run_case;
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            std::vector<Real> eig;
            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, eig);
        }
    }
    
}




template <typename ValType>
void flutter_analysis(const std::string& case_name,
                         bool with_sens,
                         std::string& par_name)  {
    
    ValType run_case;
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            std::vector<Real> eig;
            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p);
        }
    }
    
}



template <typename ValType>
void plate_analysis(const std::string& case_name,
                    const bool nonlinear,
                    bool with_sens,
                    std::string& par_name)  {
    
    ValType run_case;
    run_case.init(libMesh::QUAD4, nonlinear);
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, true);
        }
    }
    
}


template <typename ValType>
void plate_eigenvalue_analysis(const std::string& case_name,
                    const bool nonlinear,
                    bool with_sens,
                    std::string& par_name)  {
    
    ValType run_case;
    run_case.init(libMesh::QUAD4, nonlinear);
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            std::vector<Real> eig;

            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, eig);
        }
    }
    
}



template <typename ValType>
void plate_flutter_analysis(const std::string& case_name,
                            const bool nonlinear,
                            bool with_sens,
                            std::string& par_name)  {
    
    ValType run_case;
    run_case.init(libMesh::QUAD4, nonlinear);
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p);
        }
    }
    
}



template <typename ValType>
void fluid_analysis(const std::string& case_name)  {
    
    ValType run_case;
    
    libMesh::out << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    /*if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            libMesh::out
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, true);
        }
    }*/
}




template <typename ValType>
void optimization(const std::string& case_name, bool verify_grads)  {

    
    libMesh::out
    << case_name << std::endl
    << "  input.in should be provided in the working directory with"
    << " desired parameter values."
    << "  In absence of a parameter value, its default value will be used."
    << std::endl
    << "  Output per iteration is written to optimization_output.txt."
    << std::endl;
    
    GetPot infile("input.in");
    std::ofstream output;
    output.open("optimization_output.txt", std::ofstream::out);
    
    MAST::GCMMAOptimizationInterface optimizer;
    
    // create and attach sizing optimization object
    ValType func_eval(infile, output);
    __my_func_eval = &func_eval;

    if (verify_grads) {
        std::vector<Real> dvals(func_eval.n_vars());
        std::fill(dvals.begin(), dvals.end(), 0.05);
        libMesh::out << "******* Begin: Verifying gradients ***********" << std::endl;
        func_eval.verify_gradients(dvals);
        libMesh::out << "******* End: Verifying gradients ***********" << std::endl;
    }

    // attach and optimize
    optimizer.attach_function_evaluation_object(func_eval);
    optimizer.optimize();
    
    output.close();
}




template <typename ValType>
void plate_optimization(const std::string& case_name,
                        bool verify_grads,
                        bool nonlinear)  {
    
    
    libMesh::out
    << case_name << std::endl
    << "  input.in should be provided in the working directory with"
    << " desired parameter values."
    << "  In absence of a parameter value, its default value will be used."
    << std::endl
    << "  Output per iteration is written to optimization_output.txt."
    << std::endl;
    
    GetPot infile("input.in");
    std::ofstream output;
    output.open("optimization_output.txt", std::ofstream::out);
    
    MAST::GCMMAOptimizationInterface optimizer;
    
    // create and attach sizing optimization object
    ValType func_eval(output);
    __my_func_eval = &func_eval;
    func_eval.init(infile, libMesh::QUAD4, nonlinear);
    
    if (verify_grads) {
        std::vector<Real> dvals(func_eval.n_vars());
        std::fill(dvals.begin(), dvals.end(), 0.05);
        libMesh::out << "******* Begin: Verifying gradients ***********" << std::endl;
        func_eval.verify_gradients(dvals);
        libMesh::out << "******* End: Verifying gradients ***********" << std::endl;
    }
    
    // attach and optimize
    optimizer.attach_function_evaluation_object(func_eval);
    optimizer.optimize();
    
    output.close();
}




int main(int argc, char* const argv[]) {

    libMesh::LibMeshInit init(argc, argv);
    __init  = &init;
    
    // use to get arguments from the command line
    GetPot command_line(argc, argv);
    
    // look for the name that the user has requested to run
    std::string
    case_name    = command_line("run_case", ""),
    par_name     = command_line(   "param", "");
    bool
    with_sens    = command_line("with_sensitivity",    false),
    if_nonlin    = command_line("nonlinear",           false),
    verify_grads = command_line("verify_grads",        false);
    

    
    if (case_name == "bar_extension")
        analysis<MAST::BarExtension>(case_name, with_sens, par_name);
    else if (case_name == "beam_modal_analysis")
        eigenvalue_analysis<MAST::BeamModalAnalysis>(case_name, with_sens, par_name);
    else if (case_name == "beam_prestress_buckling_analysis")
        eigenvalue_analysis<MAST::BeamColumnBucklingAnalysis>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending")
        analysis<MAST::BeamBending>(case_name, with_sens, par_name);
    else if (case_name == "beam_oscillating_load")
        analysis<MAST::BeamOscillatingLoad>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending_with_offset")
        analysis<MAST::BeamBendingWithOffset>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending_thermal_stress")
        analysis<MAST::BeamBendingThermalStress>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending_optimization")
        optimization<MAST::BeamBendingSizingOptimization>(case_name, verify_grads);
    else if (case_name == "beam_bending_single_functional_optimization")
        optimization<MAST::BeamBendingSingleFunctionalSizingOptimization>
        (case_name, verify_grads);
    else if (case_name == "beam_bending_section_offset_optimization")
        optimization<MAST::BeamBendingSectionOffsetSizingOptimization>
        (case_name, verify_grads);
    else if (case_name == "beam_bending_thermal_stress_optimization")
        optimization<MAST::BeamBendingThermalStressSizingOptimization>
        (case_name, verify_grads);
    else if (case_name == "beam_piston_theory_flutter_analysis")
        flutter_analysis<MAST::BeamPistonTheoryFlutterAnalysis>(case_name, with_sens, par_name);
    else if (case_name == "beam_piston_theory_time_accurate_analysis")
        flutter_analysis<MAST::BeamPistonTheoryTimeAccurateAnalysis>(case_name, with_sens, par_name);
    else if (case_name == "membrane_extension_uniaxial")
        analysis<MAST::MembraneExtensionUniaxial>(case_name, with_sens, par_name);
    else if (case_name == "membrane_extension_biaxial")
        analysis<MAST::MembraneExtensionBiaxial>(case_name, with_sens, par_name);
    else if (case_name == "plate_modal_analysis")
        plate_eigenvalue_analysis<MAST::PlateModalAnalysis>(case_name, true, with_sens, par_name);
    else if (case_name == "plate_prestress_buckling_analysis")
        plate_eigenvalue_analysis<MAST::PlateBucklingPrestress>(case_name, true, with_sens, par_name);
    else if (case_name == "plate_bending")
        plate_analysis<MAST::PlateBending>
        (case_name, if_nonlin, with_sens, par_name);
    else if (case_name == "plate_oscillating_load")
        plate_analysis<MAST::PlateOscillatingLoad>
        (case_name, if_nonlin, with_sens, par_name);
    else if (case_name == "plate_bending_section_offset")
        plate_analysis<MAST::PlateBendingWithOffset>
        (case_name, if_nonlin, with_sens, par_name);
    else if (case_name == "plate_bending_thermal_stress")
        plate_analysis<MAST::PlateBendingThermalStress>
        (case_name, if_nonlin, with_sens, par_name);
    else if (case_name == "plate_piston_theory_flutter_analysis")
        plate_flutter_analysis<MAST::PlatePistonTheoryFlutterAnalysis>
        (case_name, false, with_sens, par_name);
    else if (case_name == "plate_thermally_stressed_piston_theory_flutter_analysis")
        plate_flutter_analysis<MAST::PlateThermallyStressedPistonTheoryFlutterAnalysis>
        (case_name, if_nonlin, with_sens, par_name);
    else if (case_name == "plate_bending_sizing_optimization")
        plate_optimization<MAST::PlateBendingSizingOptimization>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "plate_bending_single_functional_sizing_optimization")
        plate_optimization<MAST::PlateBendingSingleStressFunctionalSizingOptimization>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "plate_bending_section_offset_optimization")
        plate_optimization<MAST::PlateBendingSectionOffsetSizingOptimization>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "plate_bending_thermal_stress_optimization")
        plate_optimization<MAST::PlateBendingThermalStressSizingOptimization>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "stiffened_plate_bending_thermal_stress_optimization")
        plate_optimization<MAST::StiffenedPlateBendingThermalStressSizingOptimization>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "stiffened_plate_piston_theory_optimization")
        plate_optimization<MAST::StiffenedPlatePistonTheorySizingOptimization>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "topology_optimization_2D")
        plate_optimization<MAST::TopologyOptimization2D>
        (case_name, verify_grads, if_nonlin);
    else if (case_name == "panel_inviscid_analysis_2d")
        fluid_analysis<MAST::PanelInviscidAnalysis2D>(case_name);
    else if (case_name == "panel_inviscid_analysis_3d_half_domain")
        fluid_analysis<MAST::PanelInviscidAnalysis3DHalfDomain>(case_name);
    else if (case_name == "panel_inviscid_small_disturbance_frequency_domain_analysis_2d")
        fluid_analysis<MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis>(case_name);
    else if (case_name == "panel_inviscid_small_disturbance_frequency_domain_analysis_3d")
        fluid_analysis<MAST::PanelSmallDisturbanceFrequencyDomainInviscidAnalysis3D>(case_name);
    else if (case_name == "panel_inviscid_small_disturbance_frequency_domain_analysis_3d_half_domain")
        fluid_analysis<MAST::PanelSmallDisturbanceFrequencyDomainInviscidAnalysis3DHalfDomain>(case_name);
    else if (case_name == "beam_fsi_flutter_analysis")
        fluid_analysis<MAST::BeamEulerFSIFlutterAnalysis>(case_name);
    else if (case_name == "plate_fsi_flutter_analysis")
        fluid_analysis<MAST::PlateEulerFSIFlutterAnalysis>(case_name);
    else if (case_name == "plate_fsi_half_domain_flutter_analysis")
        fluid_analysis<MAST::PlateEulerFSIHalfDomainFlutterAnalysis>(case_name);
    else if (case_name == "bar_steady_state_conduction")
        analysis<MAST::BarSteadyState>(case_name, with_sens, par_name);
    else if (case_name == "bar_transient_conduction")
        analysis<MAST::BarTransient>(case_name, with_sens, par_name);
    else {
        libMesh::out
        << "Please run the driver with the name of example specified as: \n"
        << "   run_case=<name>"
        << "   nonlinear=<true/false>"
        << "   verify_grads=<true/false>"
        << "   with_sensitivity=<true/false>"
        << "   param=<name>\n\n"
        << "Possible values are:\n\n\n"
        << "**********************************\n"
        << "*********   STRUCTURAL   *********\n"
        << "**********************************\n"
        << "  bar_extension \n"
        << "  beam_bending \n"
        << "  beam_oscillating_load \n"
        << "  beam_modal_analysis\n"
        << "  beam_prestress_buckling_analysis\n"
        << "  beam_bending_with_offset \n"
        << "  beam_bending_thermal_stress \n"
        << "  beam_bending_optimization \n"
        << "  beam_bending_single_functional_optimization \n"
        << "  beam_bending_section_offset_optimization \n"
        << "  beam_bending_thermal_stress_optimization \n"
        << "  beam_piston_theory_flutter_analysis\n"
        << "  beam_piston_theory_time_accurate_analysis\n"
        << "  membrane_extension_uniaxial \n"
        << "  membrane_extension_biaxial \n"
        << "  plate_bending \n"
        << "  plate_oscillating_load \n"
        << "  plate_bending_section_offset \n"
        << "  plate_bending_thermal_stress \n"
        << "  plate_bending_sizing_optimization \n"
        << "  plate_bending_single_functional_sizing_optimization \n"
        << "  plate_bending_section_offset_optimization \n"
        << "  plate_bending_thermal_stress_optimization \n"
        << "  plate_modal_analysis\n"
        << "  plate_piston_theory_flutter_analysis\n"
        << "  plate_thermally_stressed_piston_theory_flutter_analysis\n"
        << "  plate_prestress_buckling_analysis\n"
        << "  stiffened_plate_bending_thermal_stress_optimization \n"
        << "  stiffened_plate_piston_theory_optimization \n"
        << "  topology_optimization_2D \n"
        << "*  The default for with_sensitivity is: false.\n"
        << "*  param is used to specify the parameter name for which sensitivity is desired.\n"
        << "*  nonlinear is used to turn on/off nonlinear stiffening in the problem.\n"
        << "*  verify_grads=true will verify the gradients of the optimization problem before calling the optimizer.\n"
        << "\n\n\n"
        << "**********************************\n"
        << "***********   FLUID   ************\n"
        << "**********************************\n"
        << "  panel_inviscid_analysis_2d \n"
        << "  panel_inviscid_analysis_3d_half_domain \n"
        << "  panel_inviscid_small_disturbance_frequency_domain_analysis_2d\n"
        << "  panel_inviscid_small_disturbance_frequency_domain_analysis_3d\n"
        << "  panel_inviscid_small_disturbance_frequency_domain_analysis_3d_half_domain\n"
        << "\n\n\n"
        << "**********************************\n"
        << "***********   CONDUCTION   *******\n"
        << "**********************************\n"
        << "  bar_steady_state_conduction \n"
        << "  bar_transient_conduction \n"
        << "\n\n\n"
        << "**********************************\n"
        << "***********   FSI     ************\n"
        << "**********************************\n"
        << "  beam_fsi_flutter_analysis \n"
        << "  plate_fsi_flutter_analysis \n"
        << "  plate_fsi_half_domain_flutter_analysis \n"
        << "\n\n\n"
        << std::endl;
    }
    
    return 0;
}
