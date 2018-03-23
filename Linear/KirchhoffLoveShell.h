// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.h
//!
//! \date Feb 25 2018
//!
//! \author ... and ... / NTNU
//!
//! \brief Class for linear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#ifndef _KIRCHHOFF_LOVE_SHELL_H
#define _KIRCHHOFF_LOVE_SHELL_H

#include "KirchhoffLove.h"


/*!
  \brief Class representing the integrand of thin shell problems.

  \details The formulation is based on Kirchhoff-Love shell theory
  and therefore requires second-derivatives of the basis functions.
*/

class KirchhoffLoveShell : public KirchhoffLove
{
public:
  //! \brief Default constructor.
<<<<<<< HEAD
  KirchhoffLoveShell() : KirchhoffLove(3) {}
=======
  KirchhoffLoveShell();
>>>>>>> 8a2d63d2b6c208393af4928fd655926d28980aac
  //! \brief Empty destructor,
  virtual ~KirchhoffLoveShell() {}

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

<<<<<<< HEAD
=======
  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf) { tracFld = tf; }
  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { fluxFld = tf; }

>>>>>>> 8a2d63d2b6c208393af4928fd655926d28980aac
  using KirchhoffLove::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  bool evalK(Matrix& EK, const FiniteElement& fe, const Vec3& X) const;

  bool formDmatrix (Matrix& Dm, Matrix& Db, const FiniteElement& fe, const Vec3& X, bool invers = false) const;

  void formMassMatrix (Matrix& EM, const Vector& N,
                       const Vec3& X, double detJW) const;

  void formBodyForce (Vector& ES, const Vector& N, size_t iP,
                      const Vec3& X, double detJW) const;

  bool formBmatrix (Matrix& Bm, Matrix& Bb, const FiniteElement& fe) const;

  using KirchhoffLove::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using KirchhoffLove::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress resultant values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  bool evalSol(Vector& s, const Vector& eV,
               const FiniteElement& fe, const Vec3& X,
               bool toLocal = false) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
<<<<<<< HEAD
  virtual size_t getNoFields(int fld = 2) const { return fld < 2 ? 3 : 5; }
=======
  virtual size_t getNoFields(int fld = 2) const { return fld < 2 ? 3 : 6; }
>>>>>>> 8a2d63d2b6c208393af4928fd655926d28980aac
  //! \brief Returns the name of the primary solution field.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;
<<<<<<< HEAD
=======

protected:
  TractionFunc* tracFld; //!< Pointer to implicit boundary traction field
  VecFunc*      fluxFld; //!< Pointer to explicit boundary traction field
>>>>>>> 8a2d63d2b6c208393af4928fd655926d28980aac
};

#endif
