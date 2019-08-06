# pragma once

//#include <cxxtest/TesSuite.h>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
//#include "PetscSetupAndFinalize.hpp"
#include "TrianglesMeshWriter.hpp"
#include "BoundaryConditionsContainerImplementation.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"

class MyTwoVariablePdeSolver
: public AbstractAssemblerSolverHybrid <2,2,2,NORMAL>,  public AbstractStaticLinearPdeSolver <2,2,2>
{
private:
  double f(double, double);
  double g(double, double);
  c_matrix <double,2*3,2*3> ComputeMatrixTerm(c_vector<double, 3>&, c_matrix<double, 2,3>&, ChastePoint<2>&, c_vector<double,2>&, c_matrix<double,2,2>&, Element<2,2>*);
  c_vector <double, 2*3> ComputeVectorTerm(c_vector<double, 3>&, c_matrix<double, 2, 2+1>&, ChastePoint<2>&, c_vector<double,2>&, c_matrix<double,2,2>&, Element<2,2>*);
  void SetupLinearSystem(Vec, bool);

public:
  MyTwoVariablePdeSolver(TetrahedralMesh<2,2>*, BoundaryConditionsContainer<2,2,2>*);
};



class ThreeParabolicPdesSolver
: public AbstractAssemblerSolverHybrid<2,2,3,NORMAL>, public AbstractDynamicLinearPdeSolver<2,2,3>
{
private:
  double g(double, ChastePoint<2>&);
  c_matrix<double, 3*3,3*3> ComputeMatrixTerm(c_vector<double,3>&, c_matrix<double,2,3>&, ChastePoint<2>&, c_vector<double,3>&, c_matrix<double,3,2>&, Element<2,2>*);
  c_vector<double,3*3> ComputeVectorTerm(c_vector<double,3>&, c_matrix<double, 2, 3>&, ChastePoint<2>&, c_vector<double,3>&, c_matrix<double,3,2>&, Element<2,2>*);
  void SetupLinearSystem(Vec, bool);

public:
  ThreeParabolicPdesSolver(TetrahedralMesh<2,2>*, BoundaryConditionsContainer<2,2,3>*);
};


class MyPde : public AbstractLinearEllipticPde<2,2>
{
private:
  c_matrix<double, 2, 2> mDiffusionTensor;

public:
  MyPde();

  double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>*);
  double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*);
  c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>&);
};


class MyNonlinearPde : public AbstractNonlinearEllipticPde<2>
{
public:
  double ComputeLinearSourceTerm(const ChastePoint<2>& );
  double ComputeNonlinearSourceTerm(const ChastePoint<2>&, double);
  c_matrix<double, 2,2> ComputeDiffusionTerm(const ChastePoint<2>&, double);
  double ComputeNonlinearSourceTermPrime(const ChastePoint<2>&, double);
  c_matrix<double, 2,2> ComputeDiffusionTermPrime(const ChastePoint<2>&, double);
};

double MyNeummanFunction(const ChastePoint<2>&);

/**
 * Two coupled PDEs defining the Schnackenberg reaction-diffusion system
 *
 * u_t = D1*del^2 u + kappa1 - kappa_1*u + kappa3*u^2*v,
 * v_t = D2*del^2 u + kappa2 - kappa3*u^2*v.
 */

class SchnackenbergCoupledPdeSystem : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<2, 2, 2>
{
private:

    double mD1;      /**< Parameter D1 in the Schnackenberg system. */
    double mD2;      /**< Parameter D2 in the Schnackenberg system. */
    double mKappa1;  /**< Parameter kappa1 in the Schnackenberg system. */
    double mKappa_1; /**< Parameter kappa_1 in the Schnackenberg system. */
    double mKappa2;  /**< Parameter kappa2 in the Schnackenberg system. */
    double mKappa3;  /**< Parameter kappa3 in the Schnackenberg system. */

public:
  /*
    SchnackenbergCoupledPdeSystem(double d1=1.0,
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
    }
    */
  SchnackenbergCoupledPdeSystem(double , double , double , double , double , double);
  //SchnackenbergCoupledPdeSystem(double =1.0, double =1.0, double =1.0, double =1.0, double =1.0, double =1.0);
  double ComputeDuDtCoefficientFunction(const ChastePoint<2>&, unsigned);
  double ComputeSourceTerm(const ChastePoint<2>&, c_vector<double,2>&, std::vector<double>&, unsigned int);
  c_matrix<double, 2, 2> ComputeDiffusionTerm(const ChastePoint<2>&, unsigned int, Element<2,2>* = NULL);
};
