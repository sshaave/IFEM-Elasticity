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
 int ndof = fe.N.size()*3;
      Matrix kem(ndof,ndof);
      Matrix keb(ndof,ndof);
      Matrix fiem(ndof,1);
      Matrix fieb(ndof,1);
      Matrix dN_ca(3,ndof);
      Matrix dM_ca(3,ndof);

      Matrix Dm;
      Matrix Db;
      if (!this->formDmatrix(Dm,Db,fe,X))
        return false;      

      Matrix dK_ca; Matrix dE_ca; Matrix3D ddE_ca; Matrix3D ddK_ca;
      Vec3 E_ca; Vec3 K_ca;
      if (!this->formBmatrix(dE_ca,dK_ca,ddE_ca,ddK_ca,E_ca,K_ca,fe))
          return false;
      Matrix E_ca_m; Matrix K_ca_m;
      E_ca_m(1,1) = E_ca(1); E_ca_m(2,1) = E_ca(2); E_ca_m(3,1) = E_ca(3);
      K_ca_m(1,1) = K_ca(1); K_ca_m(2,1) = K_ca(2); K_ca_m(3,1) = K_ca(3);
      Matrix N_ca; Matrix M_ca;
      N_ca.multiply(Dm,E_ca_m);
      M_ca.multiply(Db,K_ca_m);
      dN_ca.multiply(Dm,dE_ca);
      dM_ca.multiply(Db,dK_ca);
      // Uferdig --->
      //kem.multiply(Bm,dN_ca,true,false,true);
      //keb.multiply(Bb,dM_ca,true,false,true);
      //kem.multiply(fe.detJxW);
      //keb.multiply(fe.detJxW);
      //EK.add(kem).add(keb);
      //   <-----

      return true;
}


bool NLKirchhoffLoveShell::formBmatrix (Matrix& dE_ca, Matrix& dK_ca,
                    Matrix3D ddE_ca, Matrix3D ddK_ca,Vec3& E_ca,
                                            Vec3& K_ca,const FiniteElement& fe) const
{
    // Declaring all variables
    Vec3 g1; Vec3 g2;
    Vec3 g3; double lg3; Vec3 n; Vec3 gab; Vec3 Gab; Vec3 Bv; Vec3 bv; Matrix T;
    if (!this->getAllMetrics(g1,g2,g3,lg3,n,Gab,Bv,T, true, fe))
        return false;
    if (!this->getAllMetrics(g1,g2,g3,lg3,n,gab,bv,T, false, fe))
        return false;

    double lg3_3 = lg3*lg3*lg3;
    double lg3_5 = lg3*lg3*lg3*lg3*lg3;

    // Strain vector referred to curvilinear coor sys
    Vec3 E_cu = 0.5*(gab-Gab);
    // Strain vector referred to cartesian coor sys  -
    Matrix E_cu_m;
    E_cu_m(1,1) = E_cu(1); E_cu_m(2,1) = E_cu(2); E_cu_m(3,1) = E_cu(3);  //*
    Matrix E_ca_m;                                                        //*
    E_ca_m.multiply(E_cu_m,T);                                            //*
    E_ca(1) = E_ca_m(1,1); E_ca(2) = E_ca_m(2,1); E_ca(3) = E_ca_m(3,1);  //*
    // Curvature vector [K11,K22,K12] referred to curvilinear coor sys
    Vec3 K_cu = (Bv -bv);
    // Curvature vector referred to cart coor sys   -
    Matrix K_cu_m;
    K_cu_m(1,1) = K_cu(1); K_cu_m(2,1) = K_cu(2); K_cu_m(3,1) = K_cu(3);  //*
    Matrix K_ca_m;                                                        //*
    K_ca_m.multiply(T,K_cu_m);                                            //*
    K_ca(1) = K_ca_m(1,1); K_ca(2) = K_ca_m(2,1); K_ca(3) = K_ca_m(3,1);  //*
    // Strain
    int ndof = fe.N.size()*3;
   // Matrix dE_ca(3,ndof); // dE_ca = epsilon cartesian coordinates
   // Matrix dK_ca(3,ndof); // dK_ca = kappa cartesian coordinates
    Matrix dE_cu(3,ndof); // dE_cu = epsilon curvelinear coordinate system
    Matrix dK_cu(3,ndof); // dK_cu = kappa curvelinear

    // declaring all variables outside the loop
    double dummy; double k; int dir; Vec3 dg1;  Vec3 dg2; Matrix dg3(3,ndof);
    Matrix g3dg3lg3_3(1,ndof);Matrix g3dg3(1,ndof); Matrix dn(3,ndof); Vec3 ddg3;
    double C; double D; Vec3 ddn; Vec3 ddK_cu;
      for (int i = 1; i <= ndof; i++)
      {
        dummy = i;
        k = ceil(dummy/3);
        dir = i-3*(k-1);

        dE_cu(1,i) = fe.dNdX(k,1)*g1(dir);
        dE_cu(2,i) = fe.dNdX(k,2)*g2(dir);
        dE_cu(3,i) = 0.5*(fe.dNdX(k,1)*g2(dir) + fe.dNdX(k,2)*g1(dir));

        dg1(dir) = fe.dNdX(k,1);
        dg2(dir) = fe.dNdX(k,2);

        dg3(1,i) = dg1(2)*g2(3) - dg1(3)*g2(2) + g1(2)*dg2(3) - g1(3)*dg2(2);  // *
        dg3(2,i) = dg1(3)*g2(1) - dg1(1)*g2(3) + g1(3)*dg2(1) - g1(1)*dg2(3);  // *
        dg3(3,i) = dg1(1)*g2(2) - dg1(2)*g2(1) + g1(1)*dg2(2) - g1(2)*dg2(1);  // *
        g3dg3(1,i) = g3*dg3.getColumn(i);
        g3dg3lg3_3(1,i) = g3dg3(1,i)/(lg3_3);
        Vec3 temp = dg3.getColumn(i)/lg3 - g3*g3dg3lg3_3(1,i);                 //*
        dn(1,i) = temp(1); dn(2,i) = temp(2); dn(3,i) = temp(3);               //*
        dK_cu(1,i) = -(fe.d2NdX2(k,1,1)*n(dir) + fe.H.getColumn(1)*dn);
        dK_cu(2,i) = -(fe.d2NdX2(k,2,2)*n(dir) + fe.H.getColumn(2)*dn);
        dK_cu(3,i) = -(fe.d2NdX2(k,1,2)*n(dir) + fe.H.getColumn(3)*dn);

      } // for int i = 0; i <= ndof; i++

      dE_ca.multiply(T,dE_cu);
      dK_ca.multiply(T,dK_cu);
      // Bm.multiply(T,dE_cu); // Bm = dE_ca.multiply(T,dE_cu);
      // Bb.multiply(T,dK_cu); // Bb = dK_ca.multiply(T,dK_cu);
       Vec3 ddE_cu;
      int kr; int dirr; int ks;
      for (int i = 1; i <= ndof; i++) {
          kr = ceil(i/3);
          dirr = r-3*(kr-1);
          for (int s = 1; s<=i; s++) {
              ks = ceil(s/3);
              dirs = s-3*(ks-1);
              //strain
              if (dirr == dirs) {
                  ddE_cu(1) = fe.dNdX(kr,1)*fe.dNdX(ks,1);
                  ddE_cu(2) = fe.dNdX(kr,2)*fe.dNdX(ks,2);
                  ddE_cu(3) = 0.5*(fe.dNdX(kr,1)*fe.dNdX(ks,2)
                                   + fe.dNdX(kr,2)*fe.dNdX(ks,1));
                } else {
                  ddE_cu = ddE_cu*0;
                }
              ddE_ca(1,r,s) = T.getRow(1)*ddE_cu;
              ddE_ca(2,r,s) = T.getRow(2)*ddE_cu;
              ddE_ca(3,r,s) = T.getRow(3)*ddE_cu;

              // Curvature
              ddg3 = ddg3*0;
              dirt = 6-dirr-dirs;
              ddir = dirr-dirs;
              if (ddir == -1 || ddir = 2) {
                  ddg3(dirt) = fe.dNdX(kr,1)*fe.dNdX(ks,2) - fe.dNdX(ks,1)*fe.dNdX(kr,2);
              } else if (ddir == 1 || ddir == -2) {
                  ddg3(dirt) =  -fe.dNdX(kr,1)*fe.dNdX(ks,2) + fe.dNdX(ks,1)*fe.dNdX(kr,2);
              } // end if
              C = -( ddg3*g3+ dg3(1,r)*dg3(1,s) + dg3(2,r)*dg3(2,s) + dg3(3,r)*dg3(3,s)
                        )/lg3_3;
              D = 3*g3dg3(r)*g3dg3(s)/lg3_5;
              ddn = ddg3/lg3 - dg3.getColumn(r)*g3dg3lg3_3(s)
                  - g3dg3lg3_3(r)*dg3.getColumn(s) + C*g3 + D*g3;
              ddK_cu(1) = -(fe.d2NdX2(kr,1)*dn(dirr,s) + fe.d2NdX2(ks,1)*
                  dn(dirs,r) + fe.H(1,1)*ddn(1) + fe.H(2,1)*ddn(2) + fe.H(3,1)*ddn(3));
              ddK_cu(2) = -(fe.d2NdX2(kr,2)*dn(dirr,s) + fe.d2NdX2(ks,2)*
                  dn(dirs,r) + fe.H(1,2)*ddn(1) + fe.H(2,2)*ddn(2) + fe.H(3,2)*ddn(3));
              ddK_cu(3) = -(fe.d2NdX2(kr,3)*dn(dirr,s) + fe.d2NdX2(ks,3)*
                  dn(dirs,r) + fe.H(1,3)*ddn(1) + fe.H(2,3)*ddn(2) + fe.H(3,3)*ddn(3));
              ddK_ca(1,r,s) = T.getRow(1)*ddK_cu;
              ddK_ca(2,r,s) = T.getRow(2)*ddK_cu;
              ddK_ca(3,r,s) = T.getRow(3)*ddK_cu;
            }
        }
    return true;
}


bool NLKirchhoffLoveShell::getAllMetrics (Vec3& g1, Vec3& g2, Vec3& g3, double& lg3,Vec3& n,Vec3& gab,Vec3& Bv,Matrix& T,
                                        const bool ref, const FiniteElement& fe) const
 {
  // Har ikke definert "dA"
  if (ref) {
    Vec3 g1 = fe.G.getColumn(1);
    Vec3 g2 = fe.G.getColumn(2);
  } else {
    Vec3 g1 = fe.G.getColumn(3); // Må kanskje endre her
    Vec3 g2 = fe.G.getColumn(4);
  }

  // Basis vector g3
  g3(g1,g2);
  lg3 = g3.length();	
  n(g3);
  n.normalize();

  // Covariant metric gab
  gab(1) = g1*g1;
  gab(2) = g1*g2;
  gab(3) = g2*g2;

  if (ref)
    {
    // Contravariant metric gab_con and base vectors g_con
    double invdetgab =  1.0/(gab(1)*gab(3)-gab(2)*gab(2));
    double gab_con11 =  invdetgab*gab22;
    double gab_con12 = -invdetgab*gab12;
    double gab_con22 =  invdetgab*gab11;
    Vec3 g1_con = g1*gab_con11 + g2*gab_con12;
    Vec3 g2_con = g1*gab_con12 + g2*gab_con22;
    // Local cartesian coordinates
    Vec3 e1(g1);     e1.normalize();
    Vec3 e2(g2_con); e2.normalize();
    // Transformation matrix T from contravariant to local cartesian basis
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
  }
  Bv(1) = fe.H.getColumn(1) * n;
  Bv(2) = fe.H.getColumn(2) * n;
  Bv(3) = fe.H.getColumn(3) * n;

  return true;
}



