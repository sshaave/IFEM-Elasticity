// $Id$
//==============================================================================
//!
//! \file SIMLinEl3D.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMLinEl.h"
#include "AnalyticSolutions.h"
#include "Vec3Oper.h"


/*!
  \brief Local coordinate system for a cylinder along global z-axis.
*/

class CylinderCS : public LocalSystem
{
public:
  //! \brief Constructor printing a message making user aware of its presense.
  CylinderCS()
  {
    std::cout <<"\nLocal coordinate system: Cylindric"<< std::endl;
  }

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
    double r = hypot(X.x,X.y);
    T(1,1) = X.x/r;
    T(1,2) = X.y/r;
    T(2,1) = -T(1,2);
    T(2,2) = T(1,1);
    T(3,3) = 1.0;
    return T;
  }
};


#ifdef PRINT_CS
#include <fstream>
#endif

/*!
  \brief Local coordinate system for a cylinder along global z-axis,
  closed by a spherical cap.
*/

class CylinderSphereCS : public LocalSystem
{
public:
  //! \brief Constructor printing a message making user aware of its presense.
  CylinderSphereCS(double H = 0.0) : h(H)
  {
    std::cout <<"\nLocal coordinate system: Cylindric with Spherical cap, h="
	      << h << std::endl;
#ifdef PRINT_CS
    sn.open("nodes.dat");
    se.open("elements.dat");
    s1.open("v1.dat");
    s2.open("v2.dat");
    s3.open("v3.dat");
    sn <<"\n*NODES 4\n";
    se <<"\n*ELEMENTS 4\n%NODES #4\n"
       <<"%NO_ID\n%MAP_NODE_INDICES\n%PART_ID 4\n%POINTS\n";
    s1 <<"\n*RESULTS 31\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
    s2 <<"\n*RESULTS 32\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
    s3 <<"\n*RESULTS 33\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
  }

  virtual ~CylinderSphereCS()
  {
    std::cout <<"Finalizing VTF-output of local coordinate systems"<< std::endl;
    s1 <<"\n*GLVIEWVECTOR 2\n%NAME \"v1\"\n%STEP 1\n31\n";
    s2 <<"\n*GLVIEWVECTOR 3\n%NAME \"v2\"\n%STEP 1\n32\n";
    s3 <<"\n*GLVIEWVECTOR 4\n%NAME \"v3\"\n%STEP 1\n33\n";
#endif
  }

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
#ifdef PRINT_CS
    sn << X <<'\n';
    static int iel = 0;
    se << ++iel <<'\n';
#endif
    if (patch == 1) // Cylindric system {-z,theta,r}
    {
      T.zero();
      double r = hypot(X.x,X.y);
      T(1,3) = -1.0;
      T(2,1) = -X.y/r;
      T(2,2) =  X.x/r;
      T(3,1) =  T(2,2);
      T(3,2) = -T(2,1);
#ifdef PRINT_CS
      s1 <<"0 0 -1\n";
      s2 << T(2,1) <<" "<< T(2,2) <<" 0\n";
      s3 << T(3,1) <<" "<< T(3,2) <<" 0\n";
#endif
    }
    else // Spherical system {phi,theta,r}
    {
      Vec3 v3(X.x,X.y,X.z-h);
      v3 /= v3.length();
      double theta = atan2(X.y,X.x);
      double phi = acos(v3.z);
      Vec3 v1(cos(theta)*cos(phi),sin(theta)*cos(phi),-sin(phi));
      Vec3 v2(v3,v1);
      for (int i = 1; i <= 3; i++)
      {
	T(1,i) = v1[i-1];
	T(2,i) = v2[i-1];
	T(3,i) = v3[i-1];
      }
#ifdef PRINT_CS
      s1 << v1 <<'\n';
      s2 << v2 <<'\n';
      s3 << v3 <<'\n';
#endif
    }
    return T;
  }

private:
  double h; //!< Height above global origin of the centre of the sphere
#ifdef PRINT_CS
  mutable std::ofstream sn; //!< VTF output stream for CS nodes
  mutable std::ofstream se; //!< VTF output stream for CS point elements
  mutable std::ofstream s1; //!< VTF output stream for vector v1 of local CS
  mutable std::ofstream s2; //!< VTF output stream for vector v2 of local CS
  mutable std::ofstream s3; //!< VTF output stream for vector v3 of local CS
#endif
};


template<> bool SIMLinEl3D::parseDimSpecific (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ANASOL",6)) {
    int code = -1;
    char* cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"HOLE",4))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
                <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Hole(a,F0,nu,true));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
                <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Lshape(a,F0,nu,true));
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
                <<" F0="<< F0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTS(L,H,F0,true));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines,false);
    }
    else
    {
      std::cerr <<"  ** SIMLinEl3D::parse: Invalid analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getStressSol())
    {
      std::cout <<"Pressure code "<< code <<": Analytical traction"<< std::endl;
      this->setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new TractionField(*mySol->getStressSol());
    }
  }
  else if (!strncasecmp(keyWord,"LOCAL_SYSTEM",12)) {
    size_t i = 12;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    if (!strncasecmp(keyWord+i,"CYLINDRICZ",10))
      this->getIntegrand()->setLocalSystem(new CylinderCS);
    else if (!strncasecmp(keyWord+i,"CYLINDER+SPHERE",15))
    {
      double H = atof(keyWord+i+15);
      this->getIntegrand()->setLocalSystem(new CylinderSphereCS(H));
    }
    else
      std::cerr <<"  ** SIMLinEl3D::parse: Unsupported coordinate system: "
                << keyWord+i <<" (ignored)"<< std::endl;
  }
  else
    return false;

  return false;
}


template<> bool SIMLinEl3D::parseDimSpecific (const TiXmlElement* child)
{
  if (!strcasecmp(child->Value(),"anasol")) {
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "hole") {
      double a = 0.0, F0 = 0.0, nu = 0.0;
      utl::getAttribute(child,"a",a);
      utl::getAttribute(child,"F0",F0);
      utl::getAttribute(child,"nu",nu);
      std::cout <<"\tAnalytical solution: Hole a="<< a <<" F0="<< F0
                <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Hole(a,F0,nu,true));
    }
    else if (type == "lshape") {
      double a = 0.0, F0 = 0.0, nu = 0.0;
      utl::getAttribute(child,"a",a);
      utl::getAttribute(child,"F0",F0);
      utl::getAttribute(child,"nu",nu);
      std::cout <<"\tAnalytical solution: Lshape a="<< a <<" F0="<< F0
                <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Lshape(a,F0,nu,true));
    }
    else if (type == "cants") {
      double L = 0.0, H = 0.0, F0 = 0.0;
      utl::getAttribute(child,"L",L);
      utl::getAttribute(child,"H",H);
      utl::getAttribute(child,"F0",F0);
      std::cout <<"\tAnalytical solution: CanTS L="<< L <<" H="<< H
                <<" F0="<< F0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTS(L,H,F0,true));
    }
    else if (type == "expression") {
      std::cout <<"\tAnalytical solution: Expression"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(child,false);
    }
    else
      std::cerr <<"  ** SIMLinEl3D::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;

    // Define the analytical boundary traction field
    int code = 0;
    utl::getAttribute(child,"code",code);
    if (code > 0 && mySol && mySol->getStressSol()) {
      std::cout <<"\tNeumann code "<< code
                <<": Analytical traction"<< std::endl;
      setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new TractionField(*mySol->getStressSol());
    }
  }
  else if (!strcasecmp(child->Value(),"localsystem") && child->FirstChild()) {
    // Caution: When running adaptively, the below will cause a small memory
    // leak because the coordinate system objects are only deleted by the
    // Elasticity destructor (and not in SIMbase::clearProperties).
    if (!strcasecmp(child->FirstChild()->Value(),"cylindricz"))
      this->getIntegrand()->setLocalSystem(new CylinderCS);
    else if (!strcasecmp(child->FirstChild()->Value(),"cylinder+sphere")) {
      double H = 0.0;
      utl::getAttribute(child,"H",H);
      this->getIntegrand()->setLocalSystem(new CylinderSphereCS(H));
    }
    else
      std::cerr <<"  ** SIMLinEl3D::parse: Unsupported coordinate system: "
                << child->FirstChild()->Value() <<" (ignored)"<< std::endl;
  }
  else
    return false;

  return true;
}