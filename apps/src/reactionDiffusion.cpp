#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>

#include "reactionVessel.hpp"
#include "reactionInputProcessing.hpp"
#include "solverClass.hpp"

#include "UblasIncludes.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include "SimpleNonlinearEllipticSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "SimpleNewtonNonlinearSolver.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "HeatEquationWithSourceTerm.hpp"
#include "SimplePoissonEquation.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"

#include "RandomNumberGenerator.hpp"

#include "PetscTools.hpp"
#include "Exception.hpp"


//#include "PetscSetupAndFinalize.hpp"




#include "ExecutableSupport.hpp"

//#include "PetscException.hpp"


//#include <cxxtest/TestSuite.h>
//#include <cxxtest/GlobalFixture.h>

/*

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "UblasIncludes.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "HeatEquationWithSourceTerm.hpp"
#include "SimplePoissonEquation.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
//#include "FakePetscSetup.hpp"

#include "reactionVessel.hpp"
#include "reactionInputProcessing.hpp"

*/

int main(int argc, char *argv[])
{
    ExecutableSupport::StandardStartup(&argc, &argv);
    int exit_code = ExecutableSupport::EXIT_OK;
    try
    {
      /*

      std::string mixtureFilename = readMixtureComposition();
      std::cout<<"Mixture filename: "<<mixtureFilename<<std::endl;

      //std::cout<<"Reaction filename: "<<reactionFilename<<std::endl;

      //reactionVessel cell(reactionFilename);
      reactionVessel cell;
      std::cout<<"Print the stoichiometric matrix for the system: "<<std::endl;
      printVec(cell.getReactingSpeciesNames());
      printVecVec(cell.getReactingStoichiometry());
      std::cout<<"Print the mixture composition for the system: "<<std::endl;
      printVec(cell.getMixtureSpeciesNames());
      printVec(cell.getMixtureInitialConcentrations());
      //std::string newReactionString = "4C + 2Y <-> ATP | kf = 3.4 kr = 2.4";
      std::string newReactionString = "A <-> 2Z + 4C | kf = 1.4 kr = 2.4";
      std::cout<<"Add new reaction string: "+newReactionString<<" "<<std::endl;

      cell.addReactionFromString(newReactionString);
      std::cout<<"Print the stoichiometric matrix for the system: "<<std::endl;
      printVec(cell.getReactingSpeciesNames());
      printVecVec(cell.getReactingStoichiometry());
      std::cout<<"Print the mixture composition for the system: "<<std::endl;
      printVec(cell.getMixtureSpeciesNames());
      printVec(cell.getMixtureInitialConcentrations());

      std::string newReactionFile = "/home/chaste/projects/reactionDiffusion/apps/src/Data/reversible.txt";
      std::cout<<"Add reaction from file: reversible.txt"<<" "<<std::endl;
      cell.addReactionFromFile(newReactionFile);
      std::cout<<"Print the stoichiometric matrix for the system: "<<std::endl;
      printVec(cell.getReactingSpeciesNames());
      printVecVec(cell.getReactingStoichiometry());
      std::cout<<"Print the mixture composition for the system: "<<std::endl;
      printVec(cell.getMixtureSpeciesNames());
      printVec(cell.getMixtureInitialConcentrations());

      //printVecVec(cell.getKineticLabels());
      //printVecVec(cell.getReactionKinetics());
      std::cout<<std::endl;
      std::cout<<"Concentration change: "<<std::endl;
      printVec(cell.calculateConcentrationChange());

      */

      TetrahedralMesh<3,3> mesh;
      mesh.ConstructRegularSlabMesh(0.1,1.0,1.0,1.0);
      HeatEquationWithSourceTerm<3> thermalPde;
      BoundaryConditionsContainer<3,3,1> thermalBcc;
      thermalBcc.DefineConstantDirichletOnMeshBoundary(&mesh, 1.0);
      SimpleLinearParabolicSolver<3,3> solver(&mesh, &thermalPde, &thermalBcc);

      Vec initial_condition = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),1.0);
      solver.SetInitialCondition(initial_condition);

      // solver parameters
      double t_start = 0;
      double t_end = 1;
      double dt = 0.01;
      solver.SetTimes(t_start, t_end);
      solver.SetTimeStep(dt);

      solver.SetOutputDirectoryAndPrefix("whereOutput", "results");
      solver.SetOutputToTxt(true);
      solver.SetPrintingTimestepMultiple(10);

      Vec solution = solver.Solve();
      ReplicatableVector solution_repl(solution);

      SimplePoissonEquation<3,3> static_pde;
      SimpleLinearEllipticSolver<3,3> static_solver(&mesh, &static_pde, &thermalBcc);

      Vec static_solution = static_solver.Solve();
      ReplicatableVector static_solution_repl(static_solution);
      std::cout<<"Solution: "<<solution<<"\n";
      //for(unsigned int i=0; i<solution.size(); i++){
      //  std::cout<<"Solution "<<i<<": "<<solution[i]<<"\n";
      //}

    //  for(unsigned int i=0; i<solution.rGetTimes().size(); i++){
    //    std::cout<<solution.rGetTimes()[i] << " "<<solution.rGetSolutions()[i][0] <<" "
    //    <<solution.rGetSolutions()[i][1]<<"\n";

    //  }


      for(unsigned int i=0; i<static_solution_repl.GetSize(); i++){
  //      TS_ASSERT_DELTA(solution_repl[i], static_solution_repl[i], 1e-3);
      }

      PetscTools::Destroy(initial_condition);
      PetscTools::Destroy(solution);
      PetscTools::Destroy(static_solution);




    std::cout<<"MyTwoVariable: "<<std::endl;

    // run the 2/3 variable pde solvers

    TetrahedralMesh<2,2> mesh2;
    mesh2.ConstructRegularSlabMesh(0.01, 1.0, 1.0);
    BoundaryConditionsContainer<2,2,2> bcc2;
    bcc2.DefineZeroDirichletOnMeshBoundary(&mesh2,0); // zero dirichlet for u
    bcc2.DefineZeroDirichletOnMeshBoundary(&mesh2,1); // zero dirichlet for v

    MyTwoVariablePdeSolver solver2(&mesh2, &bcc2);

    Vec result2 = solver2.Solve();
    ReplicatableVector result_repl2(result2);

    for(unsigned int i=0; i<mesh2.GetNumNodes(); i++){
      double x = mesh2.GetNode(i)->GetPoint()[0];
      double y = mesh2.GetNode(i)->GetPoint()[1];
      double u = result_repl2[2*i];
      double v = result_repl2[2*i+1];
      double u_exact = sin(M_PI*x)*sin(M_PI*y);
      double v_exact = sin(2*M_PI*x)*sin(2*M_PI*y);
    }

    PetscTools::Destroy(result2);

    std::cout<<"ThreeParabolicPdesSolver: "<<std::endl;

    TetrahedralMesh<2,2> mesh3;
    mesh3.ConstructRegularSlabMesh(0.05, 1.0,1.0);

    BoundaryConditionsContainer<2,2,3> bcc3;

    bcc3.DefineZeroDirichletOnMeshBoundary(&mesh3, 1);
    bcc3.DefineZeroDirichletOnMeshBoundary(&mesh3, 2);

    ConstBoundaryCondition<2>* p_neumann_bc = new ConstBoundaryCondition<2>(1.0);

    TetrahedralMesh<2,2>::BoundaryElementIterator iter3 = mesh3.GetBoundaryElementIteratorBegin();

    while(iter3 < mesh3.GetBoundaryElementIteratorEnd()){
      if(fabs((*iter3)->CalculateCentroid()[0])<1e-6){
        bcc3.AddNeumannBoundaryCondition(*iter3, p_neumann_bc, 0);
      }
      iter3++;
    }

    ThreeParabolicPdesSolver solver3(&mesh3, &bcc3);

    solver3.SetTimeStep(0.01);
    solver3.SetTimes(0.0, 2.0);

    Vec initial_condition3 = PetscTools::CreateAndSetVec(3*mesh3.GetNumNodes(), 0.0);
    solver3.SetInitialCondition(initial_condition3);

    solver3.SetOutputDirectoryAndPrefix("threeVar","results3");

    solver3.SetOutputToTxt(true);
    solver3.SetPrintingTimestepMultiple(10);
    Vec result3 = solver3.Solve();
    ReplicatableVector result_repl3(result3);

    TrianglesMeshWriter<2,2> mesh_writer("threeVarProblem", "mesh3", false);
    mesh_writer.WriteFilesUsingMesh(mesh3);

    PetscTools::Destroy(initial_condition3);
    PetscTools::Destroy(result3);



    // solving elliptic pde
    std::cout<<"EllipticcPdesSolver: "<<std::endl;
    //TrianglesMeshReader<2,2> mesh_reader4("lib/mesh/data/square_128_elements");
    //TrianglesMeshReader<2,2> mesh_reader4("mesh/test/data/square_128_elements");
    TetrahedralMesh<2,2> mesh4;
    mesh4.ConstructRegularSlabMesh(0.1,1.0,1.0);

    MyPde pde4;
    BoundaryConditionsContainer<2,2,1> bcc4;
    ConstBoundaryCondition<2>* p_zero_boundary_condition = new ConstBoundaryCondition<2>(0.0);
    TetrahedralMesh<2,2>::BoundaryNodeIterator iter4 = mesh4.GetBoundaryNodeIteratorBegin();
    while(iter4 < mesh4.GetBoundaryNodeIteratorEnd()){

      double x = (*iter4) -> GetPoint()[0];
      double y = (*iter4) -> GetPoint()[1];
      if((x==0) || (y==0)){
        bcc4.AddDirichletBoundaryCondition(*iter4, p_zero_boundary_condition);
      }
      iter4++;
    }

    TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter4 = mesh4.GetBoundaryElementIteratorBegin();
    while(surf_iter4 != mesh4.GetBoundaryElementIteratorEnd()){
      unsigned int node_index4 = (*surf_iter4) -> GetNodeGlobalIndex(0);
      double x = mesh4.GetNode(node_index4) -> GetPoint()[0];
      double y = mesh4.GetNode(node_index4) -> GetPoint()[1];

      if((fabs(x-1.0) < 1e-6) || (fabs(y-1.0) < 1e-6)){
        bcc4.AddNeumannBoundaryCondition(*surf_iter4, p_zero_boundary_condition);
      }

      surf_iter4++;
    }

    SimpleLinearEllipticSolver<2,2> solver4(&mesh4, &pde4, &bcc4);

    Vec result4 = solver4.Solve();
    ReplicatableVector result_repl4(result4);
    OutputFileHandler output_file_handler("solvingLinearPde");
    out_stream p_file = output_file_handler.OpenOutputFile("linear_solution.txt");

    for(unsigned int i=0; i<result_repl4.GetSize(); i++){
      double x = mesh4.GetNode(i) -> rGetLocation()[0];
      double y = mesh4.GetNode(i) -> rGetLocation()[1];

      double u = result_repl4[i];

      (*p_file) << x << " " << y << " " << u << "\n";
    }

    PetscTools::Destroy(result4);

    std::cout<<"HeatParabolicPdesSolver: "<<std::endl;
    TetrahedralMesh<2,2> mesh5;
    mesh5.ConstructRegularSlabMesh(0.1,1.0,1.0);
    HeatEquationWithSourceTerm<2> pde5;

    BoundaryConditionsContainer<2,2,1> bcc5;
    bcc5.DefineConstantDirichletOnMeshBoundary(&mesh5, 1.0);
    SimpleLinearParabolicSolver<2,2> solver5(&mesh5, &pde5, &bcc5);

    Vec initial_condition5 = PetscTools::CreateAndSetVec(mesh5.GetNumNodes(),1.0);
    solver5.SetInitialCondition(initial_condition5);

    double t_start5 = 0;
    double t_end5 = 1;
    double dt5 = 0.01;
    solver5.SetTimes(t_start5, t_end5);
    solver5.SetTimeStep(dt5);

    solver5.SetOutputDirectoryAndPrefix("parablolicHeatSolver", "results");
    solver5.SetOutputToTxt(true);
    solver5.SetPrintingTimestepMultiple(1);

    Vec solution5 = solver5.Solve();
    ReplicatableVector solution_repl5(solution5);

    TrianglesMeshWriter<2,2> mesh_writer5("heatTxt", "mesh5", false);
    mesh_writer5.WriteFilesUsingMesh(mesh5);

    PetscTools::Destroy(initial_condition5);
    PetscTools::Destroy(solution5);

    std::cout<<"NonlinearWllipticPdesSolver: "<<std::endl;
    // solving nonlinear elliptic pde
    //TrianglesMeshReader<2, 2> mesh_reader6("mesh/test/data/square_128_elements");

    //TrianglesMeshReader<2, 2> mesh_reader6("/lib/mesh/test/datasquare_128_elements");
    TetrahedralMesh<2,2> mesh6;
    mesh6.ConstructRegularSlabMesh(0.1,1.0,1.0);

    MyNonlinearPde pde6;

    BoundaryConditionsContainer<2,2,1> bcc6;
    ConstBoundaryCondition<2>* p_zero_bc6 = new ConstBoundaryCondition<2>(0.0);
    for(TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter6 = mesh6.GetBoundaryNodeIteratorBegin(); node_iter6 != mesh6.GetBoundaryNodeIteratorEnd(); node_iter6++){
      if(fabs((*node_iter6) -> GetPoint()[1]) < 1e-12){
        bcc6.AddDirichletBoundaryCondition(*node_iter6, p_zero_bc6);
      }
    }

    FunctionalBoundaryCondition<2>* p_functional_bc = new FunctionalBoundaryCondition<2>(&MyNeummanFunction);
    for(TetrahedralMesh<2,2>::BoundaryElementIterator elt_iter6  =mesh6.GetBoundaryElementIteratorBegin(); elt_iter6 != mesh6.GetBoundaryElementIteratorEnd(); elt_iter6++){
      double y = (*elt_iter6)-> GetNodeLocation(0,1);
      if(fabs(y-1.0)<1e-12){
        bcc6.AddNeumannBoundaryCondition(*elt_iter6, p_functional_bc);
      }else{
        bcc6.AddNeumannBoundaryCondition(*elt_iter6, p_zero_bc6);
      }
    }

    SimpleNonlinearEllipticSolver<2,2> solver6(&mesh6, &pde6, &bcc6);

    Vec initial_guess6 = PetscTools::CreateAndSetVec(mesh6.GetNumNodes(), 0.25);

    SimpleNewtonNonlinearSolver newton_solver;
    solver6.SetNonlinearSolver(&newton_solver);
    newton_solver.SetTolerance(1e-10);
    newton_solver.SetWriteStats();

    Vec answer6 = solver6.Solve(initial_guess6);


    PetscTools::Destroy(initial_guess6);
    PetscTools::Destroy(answer6);

    std::cout<<"ChemicalPdesSolver: "<<std::endl;

    //TrianglesMeshReader<2,2> mesh_reader7("mesh/test/data/butterfly");
    TetrahedralMesh<2,2> mesh7;
    mesh7.ConstructRegularSlabMesh(0.1,1.0,1.0);
    //mesh7.Scale(0.2,0.2);

    SchnackenbergCoupledPdeSystem pde7(1e-3, 1e-2, 0.1,0.2, 1.0, 1.0);

    BoundaryConditionsContainer<2,2,2> bcc7;
    ConstBoundaryCondition<2>* p_bc_for_u7 = new ConstBoundaryCondition<2> (2.0);
    ConstBoundaryCondition<2>* p_bc_for_v7 = new ConstBoundaryCondition<2> (0.75);

    for(TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter7 = mesh7.GetBoundaryNodeIteratorBegin(); node_iter7 != mesh7.GetBoundaryNodeIteratorEnd(); node_iter7++){
      bcc7.AddDirichletBoundaryCondition(*node_iter7, p_bc_for_u7, 0);
      bcc7.AddDirichletBoundaryCondition(*node_iter7, p_bc_for_v7, 1);
    }

    LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,2> solver7(&mesh7, &pde7,&bcc7);

    double t_end7 = 10;
    solver7.SetTimes(0, t_end7);
    solver7.SetTimeStep(1e-1);
    solver7.SetSamplingTimeStep(1);
    solver7.SetOutputDirectoryAndPrefix("SchnackenbergButterflyTxt", "resultSlab");
    solver7.SetOutputToTxt(true);
    solver7.SetPrintingTimestepMultiple(1);

    std::vector<double> init_conds(2*mesh7.GetNumNodes());

    for(unsigned int i=0; i<mesh7.GetNumNodes(); i++){
      init_conds[2*i] = fabs(2.0 + RandomNumberGenerator::Instance()->ranf());
      init_conds[2*i+1] = fabs(0.75 + RandomNumberGenerator::Instance()->ranf());
    }

    Vec initial_condition7 = PetscTools::CreateVec(init_conds);
    solver7.SetInitialCondition(initial_condition7);
    //solver7.SolveAndWriteResultsToFile();

    Vec result7 = solver7.Solve();
    ReplicatableVector result_repl7(result7);
    TrianglesMeshWriter<2,2> mesh_writer7("SchnackenbergButterflyTxt", "mesh7", false);
    mesh_writer7.WriteFilesUsingMesh(mesh7);

    PetscTools::Destroy(initial_condition7);
    PetscTools::Destroy(result7);







    std::cout<<"File Complete"<<std::endl;






    }

    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
