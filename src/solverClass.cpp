#include "solverClass.hpp"

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "TrianglesMeshWriter.hpp"
#include "BoundaryConditionsContainerImplementation.hpp"
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"

double MyTwoVariablePdeSolver::f(double x, double y){
  return -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y) + sin(2*M_PI*x)*sin(2*M_PI*y);
};

double MyTwoVariablePdeSolver::g(double x, double y){
  return -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y) + sin(M_PI*x)*sin(M_PI*y);
};

c_matrix <double,2*3,2*3>  MyTwoVariablePdeSolver::ComputeMatrixTerm(c_vector<double, 3>& rPhi, c_matrix<double, 2,3>& rGradPhi, ChastePoint<2>& rX, c_vector<double,2>& rU, c_matrix<double,2,2>& rGradU, Element<2,2>* PElement){
  c_matrix<double,2*3,2*3> ret  = zero_matrix<double>(2*3,2*3);
  for(unsigned int i=0; i<3; i++){
    for(unsigned int j=0; j<3; j++){
      for(unsigned int k=0; k<2; k++){
        // stiffness matrix on diagonal blocks
        ret(2*i, 2*j) += rGradPhi(k,i)*rGradPhi(k,j);
        ret(2*i+1, 2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
      }
      // (negative) mass matrix on off-diagonal blocks
      ret(2*i+1, 2*j) = -rPhi(i)*rPhi(j);
      ret(2*i, 2*j+1) = -rPhi(i)*rPhi(j);
    }
  }
  return ret;
};

c_vector <double, 2*3> MyTwoVariablePdeSolver::ComputeVectorTerm(c_vector<double, 3>& rPhi, c_matrix<double, 2, 2+1>& rGradPhi, ChastePoint<2>& rX, c_vector<double,2>& rU, c_matrix<double,2,2>& rGradU, Element<2,2>* pElement){
  c_vector<double,2*3> ret;
  for(unsigned int i=0; i<3; i++){
    ret(2*i) = -f(rX[0], rX[1])*rPhi(i);
    ret(2*i+1) = -g(rX[0],rX[1])*rPhi(i);
  }
  return ret;
};

void MyTwoVariablePdeSolver::SetupLinearSystem(Vec currentSolution, bool computeMatrix){
  SetupGivenLinearSystem(currentSolution, computeMatrix, this -> mpLinearSystem);
};

MyTwoVariablePdeSolver::MyTwoVariablePdeSolver(TetrahedralMesh<2,2>* pMesh, BoundaryConditionsContainer<2,2,2>* pBoundaryConditions)
: AbstractAssemblerSolverHybrid<2,2,2,NORMAL>(pMesh, pBoundaryConditions), AbstractStaticLinearPdeSolver<2,2,2>(pMesh){

};


double ThreeParabolicPdesSolver::g(double t, ChastePoint<2>& rX){
  return t*(rX[0]>0.5);
};

c_matrix<double, 3*3,3*3> ThreeParabolicPdesSolver::ComputeMatrixTerm(c_vector<double,3>& rPhi, c_matrix<double,2,3>& rGradPhi, ChastePoint<2>& rX, c_vector<double,3>& rU, c_matrix<double,3,2>& rGradU, Element<2,2>* pElement){
  c_matrix<double, 9,9> ret = zero_matrix<double>(9,9);

  // get current timestep
  double dt = PdeSimulationTime::GetPdeTimeStep();

  for(unsigned int i=0; i<3; i++){
    for(unsigned int j=0; j<3; j++){
      // mass matrix on the diagonal blocks
      ret(3*i, 3*j) = rPhi(i)*rPhi(j)/dt;
      ret(3*i+1, 3*j+1) = rPhi(i)*rPhi(j)/dt;
      ret(3*i+2, 3*j+2) = rPhi(i)*rPhi(j)/dt;

      // mass matrix on some off-diagonal blocks
      ret(3*i, 3*j+1) = -rPhi(i)*rPhi(j);
      ret(3*i+1, 3*j) = -rPhi(i)*rPhi(j);
      ret(3*i+1, 3*j+2) = -2*rPhi(i)*rPhi(j);

      // stiffness matrix on the diagonal blocks
      for(unsigned int dim=0; dim<2; dim++){
        ret(3*i, 3*j) += rGradPhi(dim,i)*rGradPhi(dim,j);
        ret(3*i+1, 3*j+1) += rGradPhi(dim,i)*rGradPhi(dim,j);
        ret(3*i+2, 3*j+2) += rGradPhi(dim,i)*rGradPhi(dim,j);
      }


    }
  }

  return ret;
};

c_vector<double,3*3> ThreeParabolicPdesSolver::ComputeVectorTerm(c_vector<double,3>& rPhi, c_matrix<double, 2, 3>& rGradPhi, ChastePoint<2>& rX, c_vector<double,3>& rU, c_matrix<double,3,2>& rGradU, Element<2,2>* pElement){
  c_vector<double,3*3> ret;

  // get u,v,w out of provided paramters
  double u = rU(0);
  double v = rU(1);
  double w = rU(2);

  double t = PdeSimulationTime::GetTime();
  double dt = PdeSimulationTime::GetPdeTimeStep();
  double inverse_dt = PdeSimulationTime::GetPdeTimeStepInverse();

  for(unsigned int i=0; i<3; i++){
    ret(3*i) = u* inverse_dt * rPhi(i);
    ret(3*i +1) = v* inverse_dt * rPhi(i);
    ret(3*i+2) = (w* inverse_dt + g(t+dt,rX) * rPhi(i));
  }
  return ret;
};

void ThreeParabolicPdesSolver::SetupLinearSystem(Vec currentSolution, bool computeMatrix){
  SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
};

ThreeParabolicPdesSolver::ThreeParabolicPdesSolver(TetrahedralMesh<2,2>* pMesh, BoundaryConditionsContainer<2,2,3>* pBoundaryConditions)
: AbstractAssemblerSolverHybrid<2,2,3,NORMAL>(pMesh,pBoundaryConditions), AbstractDynamicLinearPdeSolver<2,2,3>(pMesh){
  this->mMatrixIsConstant = true;
};

MyPde::MyPde(){
  mDiffusionTensor(0,0) = 2.0;
  mDiffusionTensor(0,1) = 0.0;
  mDiffusionTensor(1,0) = 0.0;
  mDiffusionTensor(1,1) = 1.0;
};

double MyPde::ComputeConstantInUSourceTerm(const ChastePoint<2>& rX, Element<2,2>* pElement){
  return rX[0]*rX[0] + rX[1]*rX[1];
};

double MyPde::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& rX, Element<2,2>* pElement){
  return 1.0;
};

c_matrix<double,2,2> MyPde::ComputeDiffusionTerm(const ChastePoint<2>& rX){
  return mDiffusionTensor;
};


double MyNonlinearPde::ComputeLinearSourceTerm(const ChastePoint<2>& rX){
  return 1.0;
};

double MyNonlinearPde::ComputeNonlinearSourceTerm(const ChastePoint<2>& rX, double u){
  return 0.0;
};

c_matrix<double, 2,2> MyNonlinearPde::ComputeDiffusionTerm(const ChastePoint<2>& rX, double u){
  return identity_matrix<double>(2)*u;
};

double MyNonlinearPde::ComputeNonlinearSourceTermPrime(const ChastePoint<2>& rX, double u){
  return 0.0;
};

c_matrix<double, 2,2> MyNonlinearPde::ComputeDiffusionTermPrime(const ChastePoint<2>& rX, double u){
  return identity_matrix<double>(2);
};

double MyNeummanFunction(const ChastePoint<2>& rX){
  return rX[1];
};



/*
SchnackenbergCoupledPdeSystem::SchnackenbergCoupledPdeSystem(double d1=1.0,
                              double d2=1.0,
                              double kappa1=1.0,
                              double kappa_1=1.0,
                              double kappa2=1.0,
                              double kappa3=1.0)
    : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 2>(),
      mD1(d1),
      mD2(d2),
      mKappa1(kappa1),
      mKappa_1(kappa_1),
      mKappa2(kappa2),
      mKappa3(kappa3)
{
};
*/
SchnackenbergCoupledPdeSystem::SchnackenbergCoupledPdeSystem(double d1, double d2, double kappa1, double kappa_1, double kappa2, double kappa3)
: AbstractLinearParabolicPdeSystemForCoupledOdeSystem<2, 2, 2>(){
      mD1 = d1;
      mD2 = d2;
      mKappa1 = kappa1;
      mKappa_1 = kappa_1;
      mKappa2 = kappa2;
      mKappa3 = kappa3;
};


double SchnackenbergCoupledPdeSystem::ComputeDuDtCoefficientFunction(const ChastePoint<2>& rX, unsigned int index){
    return 1.0;
};

double SchnackenbergCoupledPdeSystem::ComputeSourceTerm(const ChastePoint<2>& rX, c_vector<double,2>& rU, std::vector<double>& rOdeSolution, unsigned int pdeIndex){
    assert(pdeIndex == 0 || pdeIndex == 1);

    double source_term;
    if (pdeIndex == 0)
    {
        source_term = mKappa1 - mKappa_1*rU(0) + mKappa3*rU(1)*rU(0)*rU(0);
    }
    else // pdeIndex == 1
    {
        source_term = mKappa2 - mKappa3*rU(1)*rU(0)*rU(0);
    }
    return source_term;
}

c_matrix<double, 2, 2> SchnackenbergCoupledPdeSystem::ComputeDiffusionTerm(const ChastePoint<2>& rX, unsigned int pdeIndex, Element<2,2>* pElement){
    assert(pdeIndex == 0 || pdeIndex == 1);

    c_matrix<double, 2, 2> diffusion_term;
    if (pdeIndex == 0)
    {
        diffusion_term = mD1*identity_matrix<double>(2);
    }
    else // pdeIndex == 1
    {
        diffusion_term = mD2*identity_matrix<double>(2);
    }
    return diffusion_term;
}

/*
template<unsigned int SPACE_DIM> c_matrix<double, SPACE_DIM, SPACE_DIM> SchnackenbergCoupledPdeSystem<SPACE_DIM>::ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned int pdeIndex, Element<SPACE_DIM,SPACE_DIM>* pElement){
    assert(pdeIndex == 0 || pdeIndex == 1);

    c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;
    if (pdeIndex == 0)
    {
        diffusion_term = mD1*identity_matrix<double>(SPACE_DIM);
    }
    else // pdeIndex == 1
    {
        diffusion_term = mD2*identity_matrix<double>(SPACE_DIM);
    }
    return diffusion_term;
}
*/
