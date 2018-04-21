// $Id$
//==============================================================================
//!
//! \file NLKirchhoffLoveShell.C
//!
//! \date Apr 20 2018
//!
//! \author ... and ... / NTNU
//!
//! \brief Class for nonlinear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#include "NLKirchhoffLoveShell.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"


void NLKirchhoffLoveShell::setMode (SIM::SolutionMode mode)
{
  this->KirchhoffLoveShell::setMode(mode);
  if (mode == SIM::STATIC || mode == SIM::RHS_ONLY)
    iS = 1;
}


int NLKirchhoffLoveShell::getIntegrandType () const
{
  return SECOND_DERIVATIVES | UPDATED_NODES;
}


bool NLKirchhoffLoveShell::evalInt (LocalIntegral& elmInt,
                                    const FiniteElement& fe,
                                    const Vec3& X) const
{
  if (elmInt.vec.size() < 2)
    return false;

  // Co-variant basis and Hessian in deformed configuration
  Matrix Gd, Hd;
  Matrix3D Hess;
  Gd.multiplyMat(elmInt.vec.back(),fe.dNdX);
  if (Hess.multiplyMat(elmInt.vec.back(),fe.d2NdX2))
    utl::Hessian(Hess,Hd);
  else
    return false;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eK && iS) // Integrate the stiffness matrix and internal forces
    this->evalKandS(elMat.A[eK-1],elMat.b[iS-1],fe,fe.G,Gd,fe.H,Hd,X);
  else if (iS) // Integrate the internal forces only
  {
    Matrix dummyEK;
    this->evalKandS(dummyEK,elMat.b[iS-1],fe,fe.G,Gd,fe.H,Hd,X);
  }

  if (eS) // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,fe.detJxW);

  return true;
}


bool NLKirchhoffLoveShell::evalKandS (Matrix& EK, Vector& ES,
                                      const FiniteElement& fe,
                                      const Matrix& G0, const Matrix& Gn,
                                      const Matrix& H0, const Matrix& Hn,
                                      const Vec3& X) const
{
  Matrix Dm, Db;
  if (!this->formDmatrix(Dm,Db,fe,X))
    return false;

  // TODO: Replace the below with the nonlinear implementation
  Matrix Bm, Bb;
  if (!this->formBmatrix(Bm,Bb,fe))
    return false;

  Matrix dN_ca, dM_ca;
  dN_ca.multiply(Dm,Bm);
  dM_ca.multiply(Db,Bb);

  Matrix kem, keb;
  kem.multiply(Bm,dN_ca,true);
  keb.multiply(Bb,dM_ca,true);

  if (!EK.empty())
    EK.add(kem,fe.detJxW).add(keb,fe.detJxW);

  return true;
}


bool NLKirchhoffLoveShell::formBmatrix (Matrix& Bm, Matrix& Bb,
                                        const FiniteElement& fe) const
{
  // TODO: Replace this with the nonlinear formulation

  // Covariant basis vectors
  Vec3 g1 = fe.G.getColumn(1);
  Vec3 g2 = fe.G.getColumn(2);

  // Basis vector g3
  Vec3 g3(g1,g2);
  double lg3 = g3.length();
  Vec3 n(g3); n.normalize();

  // Covariant metric gab
  double gab11 = g1*g1;
  double gab12 = g1*g2;
  double gab22 = g2*g2;

  // Contravariant metric gab_con and base vectors g_con
  double invdetgab =  1.0/(gab11*gab22-gab12*gab12);
  double gab_con11 =  invdetgab*gab22;
  double gab_con12 = -invdetgab*gab12;
  double gab_con22 =  invdetgab*gab11;
  Vec3 g1_con = g1*gab_con11 + g2*gab_con12;
  Vec3 g2_con = g1*gab_con12 + g2*gab_con22;

  // Local cartesian coordinates
  Vec3 e1(g1);     e1.normalize();
  Vec3 e2(g2_con); e2.normalize();

  // Transformation matrix T from contravariant to local cartesian basis
  Matrix T(3,3);
  double eg11 = e1*g1_con;
  double eg12 = e1*g2_con;
  double eg21 = e2*g1_con;
  double eg22 = e2*g2_con;
  T(1,1) =      eg11*eg11;
  T(1,2) =      eg12*eg12;
  T(1,3) = 2.0* eg11*eg12;
  T(2,1) =      eg21*eg21;
  T(2,2) =      eg22*eg22;
  T(2,3) = 2.0* eg21*eg22;
  T(3,1) = 2.0* eg11*eg21;
  T(3,2) = 2.0* eg12*eg22;
  T(3,3) = 2.0*(eg11*eg22+eg12*eg21);

  // Strain
  int ndof = fe.N.size()*3;
  Matrix dE_ca(3,ndof); // dE_ca = epsilon cartesian coordinates
  Matrix dK_ca(3,ndof); // dK_ca = kappa cartesian coordinates
  Matrix dE_cu(3,ndof); // dE_cu = epsilon curvelinear coordinate system
  Matrix dK_cu(3,ndof); // dK_cu = kappa curvelinear

  for (int i = 1; i <= ndof; i++)
  {
    double dummy = i;
    double k = ceil(dummy/3);
    int dir = i-3*(k-1);

    dE_cu(1,i) = fe.dNdX(k,1)*g1(dir);
    dE_cu(2,i) = fe.dNdX(k,2)*g2(dir);
    dE_cu(3,i) = 0.5*(fe.dNdX(k,1)*g2(dir) + fe.dNdX(k,2)*g1(dir));

    Vec3 dg1, dg2;
    dg1(dir) = fe.dNdX(k,1);
    dg2(dir) = fe.dNdX(k,2);
    Vec3 dg3 = Vec3(g1,dg2) + Vec3(dg1,g2);

    double g3dg3lg3_3 = g3*dg3/(lg3*lg3*lg3); // eq 5.31 last part

    // eq 5.31 dn = a3,rs (kanskje)
    Vec3 dn = dg3/lg3 - g3*g3dg3lg3_3;

    dK_cu(1,i) = -(fe.d2NdX2(k,1,1)*n(dir) + fe.H.getColumn(1)*dn);
    dK_cu(2,i) = -(fe.d2NdX2(k,2,2)*n(dir) + fe.H.getColumn(2)*dn);
    dK_cu(3,i) = -(fe.d2NdX2(k,1,2)*n(dir) + fe.H.getColumn(3)*dn);
  }

  Bm.multiply(T,dE_cu); // Bm = dE_ca.multiply(T,dE_cu);
  Bb.multiply(T,dK_cu); // Bb = dK_ca.multiply(T,dK_cu);

  return true;
}
