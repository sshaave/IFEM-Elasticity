// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.C
//!
//! \date Feb 25 2018
//!
//! \author ... and ... / NTNU
//!
//! \brief Class for linear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#include "KirchhoffLoveShell.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "Vec3Oper.h"
#include "IFEM.h"


void KirchhoffLoveShell::printLog () const
{
  IFEM::cout <<"KirchhoffLoveShell: thickness = "<< thickness
             <<", gravity = "<< gravity << std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<KirchhoffLoveShell*>(this)->material = &defaultMat;
  }

  material->printLog();
}


bool KirchhoffLoveShell::evalInt (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Covariant basis vectors
  Vec3 g1 = fe.G.getColumn(1);
  Vec3 g2 = fe.G.getColumn(2);

  // Basis vector g3
  Vec3 g3(g1,g2);
  double lg3 = g3.length();

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

  if (eM) // Integrate the mass matrix
    ; // TODO: Include the mass-matrix terms in elmMat.A[eM-1]

  if (eK) // Integrate the stiffness matrix
    ; // TODO: Include the stiffness-matrix terms in elmMat.A[eK-1]

  if (eS) // Integrate the load vector due to gravitation and other body forces
    ; // TODO: Include the pressure/gravity load terms in elmMat.b[eS-1]

  return true;
}


bool KirchhoffLoveShell::evalBou (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X, const Vec3& normal) const
{
  // TODO (if you want to support Neumann boundary conditions)
  std::cerr <<" *** KirchhoffLoveShell::evalBou not implemented."<< std::endl;
  return false;
}


bool KirchhoffLoveShell::evalSol (Vector& s,
                                  const FiniteElement& fe, const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vector eV;
  if (!primsol.empty() && !primsol.front().empty())
  {
    int ierr = utl::gather(MNPC,3,primsol.front(),eV);
    if (ierr > 0)
    {
      std::cerr <<" *** KirchhoffLoveShell::evalSol: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;
      return false;
    }
  }

  // Evaluate the stress resultant tensor
  return this->evalSol(s,eV,fe,X,true);
}


bool KirchhoffLoveShell::evalSol (Vector& s, const Vector& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  // TODO (if you want to support postprocessing of stress resultants)
  std::cerr <<" *** KirchhoffLoveShell::evalSol not implemented."<< std::endl;
  return false;
}


std::string KirchhoffLoveShell::getField1Name (size_t i,
					       const char* prefix) const
{
  if (i >= 3) return "";

  char name = 'u'+i;
  if (!prefix)
    return std::string(name,1);

  return prefix + std::string(" ") + std::string(name,1);
}


std::string KirchhoffLoveShell::getField2Name (size_t i,
					       const char* prefix) const
{
  if (i >= 6) return "";

  static const char* s[6] = { "n_xx", "n_yy", "n_xy", "m_xx", "m_yy", "m_xy" };

  if (!prefix)
    return s[i];

  return prefix + std::string(" ") + s[i];
}
