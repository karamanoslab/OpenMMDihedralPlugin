%module dihedralplugin

%include "typemaps.i"

%ignore OpenMM::Force::Force;

namespace OpenMM {
	class Force {
	public:
	   virtual ~Force() ;
	
	   int getForceGroup() const ;
	   void setForceGroup(int group) ;
	   const std::string & getName() const ;
	   void setName(const std::string &name) ;
	   virtual bool usesPeriodicBoundaryConditions() const ;
	};
}


%typemap(in) const OpenMM::RealVec& {
    if (!PySequence_Check($input) || PySequence_Size($input) != 3) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of 3 floats");
        return NULL;
    }
    double x = PyFloat_AsDouble(PySequence_GetItem($input, 0));
    double y = PyFloat_AsDouble(PySequence_GetItem($input, 1));
    double z = PyFloat_AsDouble(PySequence_GetItem($input, 2));
    $1 = new OpenMM::RealVec(x, y, z);
}

%typemap(freearg) const OpenMM::RealVec& {
    delete $1;
}

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include <vector>
#include <string>
#include <stdexcept>
#include "openmm/Force.h"
#include "openmm/Context.h"
#include "DihedralForce.h"
#include "openmm/Vec3.h"
# include "OpenMM.h"
# include "OpenMMAmoeba.h"
# include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}


%pythoncode %{
import openmm as mm
from openmm import unit
%}

/*
 * Add units to function outputs.
*/
%pythonappend DihedralPlugin::DihedralForce::getDihedralParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                               double& ka1, double& phia1, double& epsd, double& ka2, double& phia2, double& eps0, double& kb1, double& phib1, double& eps1, double& kb2, double& phib2, double& eps2) const %{
    val[4] = unit.Quantity(val[4], 1/unit.radian**2)   
    val[5] = unit.Quantity(val[5], unit.radian)        
    val[7] = unit.Quantity(val[7], 1/unit.radian**4)   
    val[8] = unit.Quantity(val[8], unit.radian)        
    val[10] = unit.Quantity(val[10], 1/unit.radian**2) 
    val[11] = unit.Quantity(val[11], unit.radian)     
    val[13] = unit.Quantity(val[13], 1/unit.radian**4) 
    val[14] = unit.Quantity(val[14], unit.radian)      
    
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace DihedralPlugin {

class DihedralForce : public OpenMM::Force {
public:
    DihedralForce();

    int getNumDihedrals() const;

    int addDihedral(int particle1, int particle2, int particle3, int particle4,
                     double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2);

    void setDihedralParameters(int index, int particle1, int particle2, int particle3, int particle4,
                     double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2);

    double calcDihedralAngleFromIndexes(const OpenMM::Context& context, int p1, int p2, int p3, int p4);
    
    double calcDihedralAngleFromCoords(const OpenMM::RealVec& pos1,
                                   const OpenMM::RealVec& pos2,
                                   const OpenMM::RealVec& pos3,
                                   const OpenMM::RealVec& pos4);

    
    void updateParametersInContext(OpenMM::Context& context);
    
    void setUsesPeriodicBoundaryConditions(bool periodic);

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply int& OUTPUT {int& particle1};
    %apply int& OUTPUT {int& particle2};
    %apply int& OUTPUT {int& particle3};
    %apply int& OUTPUT {int& particle4};
    %apply double& OUTPUT {double& ka1};
    %apply double& OUTPUT {double& phia1};
    %apply double& OUTPUT {double& epsd};
    %apply double& OUTPUT {double& ka2};
    %apply double& OUTPUT {double& phia2};
    %apply double& OUTPUT {double& eps0};
    %apply double& OUTPUT {double& kb1};
    %apply double& OUTPUT {double& phib1};
    %apply double& OUTPUT {double& eps1};
    %apply double& OUTPUT {double& kb2};
    %apply double& OUTPUT {double& phib2};
    %apply double& OUTPUT {double& eps2};
    
    void getDihedralParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                               double& ka1, double& phia1, double& epsd, double& ka2, double& phia2, double& eps0, double& kb1, double& phib1, double& eps1, double& kb2, double& phib2, double& eps2) const;
    %clear int& particle1;
    %clear int& particle2;
    %clear int&   particle3;
    %clear int&  particle4;
    %clear double&   ka1;
    %clear double&   phia1;
    %clear double&   epsd;
    %clear double&  ka2;
    %clear double& phia2;
    %clear double&  eps0;
    %clear double&   kb1;
    %clear double&   phib1;
    %clear double&   eps1;
    %clear double&   kb2;
    %clear double& phib2;
    %clear double&  eps2;
   

    /*
     * Add methods for casting a Force to an DihedralForce.
    */
    %extend {
        static DihedralPlugin::DihedralForce& cast(OpenMM::Force& force) {
            return dynamic_cast<DihedralPlugin::DihedralForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<DihedralPlugin::DihedralForce*>(&force) != NULL);
        }
    }
};

}
