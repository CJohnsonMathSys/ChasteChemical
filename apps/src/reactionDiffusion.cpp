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
#include "ExecutableSupport.hpp"

#include "reactionVessel.hpp"
#include "reactionInputProcessing.hpp"



int main(int argc, char *argv[])
{
    ExecutableSupport::StandardStartup(&argc, &argv);
    int exit_code = ExecutableSupport::EXIT_OK;
    try
    {


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
