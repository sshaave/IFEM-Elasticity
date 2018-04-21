// $Id$
//==============================================================================
//!
//! \file SIMKLShell.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love shells.
//!
//==============================================================================

#ifndef _SIM_KL_SHELL_H
#define _SIM_KL_SHELL_H

#include "SIMElasticity.h"
#include "SIM2D.h"

class KirchhoffLove;
class ScalarFunc;


/*!
  \brief Driver class for isogeometric FEM analysis of Kirchhoff-Love shells.
*/

class SIMKLShell : public SIMElasticity<SIM2D>
{
protected:
  typedef std::pair<int,int> Ipair; //!< Convenience declaration

  /*!
    \brief Struct defining a nodal point load.
  */
  struct PointLoad
  {
    size_t      patch; //!< Patch index [0,nPatch>
    Ipair       ldof;  //!< Local node/DOF number of the closest node
    double      xi[2]; //!< Parameters of the point (u,v)
    Vec3        X;     //!< Spatial coordinates of the point
    ScalarFunc* p;     //!< Load magnitude
    //! \brief Default constructor.
    PointLoad() : patch(0), ldof(Ipair(0,1)), p(nullptr) { xi[0]=xi[1] = 0.0; }
    //! \brief The destructor deletes the load function.
    ~PointLoad() { delete p; }
  };

  typedef std::vector<PointLoad> PloadVec; //!< Point load container

public:
  //! \brief Default constructor.
  explicit SIMKLShell(bool shell = true);
  //! \brief Destructor.
  virtual ~SIMKLShell();

  //! \brief Returns that no analytical solution is available.
  virtual bool haveAnaSol() const { return false; }

protected:
  //! \brief Returns the actual integrand.
  virtual KirchhoffLove* getProblem(int version = 1);
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand();
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual void preprocessA();
  //! \brief Performs some pre-processing tasks on the FE model.
  virtual bool preprocessB();

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);
  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd);
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase*, const TimeDomain&);

  //! \brief Computes problem-dependent external energy contributions.
  virtual double externalEnergy(const Vectors& u, const TimeDomain& time) const;

protected:
  RealArray tVec;     //!< Shell thickness data
  PloadVec  myLoads;  //!< Nodal point loads
  int       aCode[3]; //!< Analytical BC codes (used by destructor)
};

#endif
