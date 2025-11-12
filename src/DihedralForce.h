#ifndef OPENMM_DIHEDRALFORCE_H_
#define OPENMM_DIHEDRALFORCE_H_


#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "windowsExportExample.h"
#include "openmm/reference/RealVec.h"


namespace DihedralPlugin {

/**
 * This class implements an anharmonic bond force of the form E(r)=k*(r-length)^4.  It exists to
 * serve as an example of how to write plugins.
 */

class OPENMM_EXPORT_EXAMPLE DihedralForce : public OpenMM::Force {
public:
    /**
     * Create an ExampleForce.
     */
    DihedralForce();
    /**
     * Get the number of bond stretch terms in the potential function
     */
//     int getNumBonds() const {
//         return bonds.size();
//     }
    
    int getNumDihedrals() const {
         return dihedrals.size(); 
         }
    /**
     * Add a bond term to the force.
     *
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the force constant for the bond, measured in kJ/mol/nm^4
     * @return the index of the bond that was added
     */
    int addDihedral(int p1, int p2, int p3, int p4,
                     double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2);
    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to get parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond, measured in kJ/mol/nm^4
     */
    void getDihedralParameters(int index, int& p1, int& p2, int& p3, int& p4,
                               double& ka1, double& phia1, double& epsd, double& ka2, double& phia2, double& eps0, double& kb1, double& phib1, double& eps1, double& kb2, double& phib2, double& eps2) const;
    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to set parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond, measured in kJ/mol/nm^4
     */
   
    void setDihedralParameters(int index, int p1, int p2, int p3, int p4,
                               double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2);
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     */ 
     
    double calcDihedralAngleFromIndexes(const OpenMM::Context& context, int& p1, int& p2, int& p3, int& p4);
    
    double calcDihedralAngleFromCoords(const OpenMM::RealVec& pos1,
                                   const OpenMM::RealVec& pos2,
                                   const OpenMM::RealVec& pos3,
                                   const OpenMM::RealVec& pos4);

    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    
    void setUsesPeriodicBoundaryConditions(bool periodic);
    
    bool usesPeriodicBoundaryConditions() const;
    
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class DihedralInfo;
    std::vector<DihedralInfo> dihedrals;
    bool usePeriodic; 
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
// class DihedralForce::DihedralInfo {
// public:
//     int particle1, particle2;
//     double length, k;
//     BondInfo() {
//         particle1 = particle2 = -1;
//         length = k = 0.0;
//     }
//     BondInfo(int particle1, int particle2, double length, double k) :
//         particle1(particle1), particle2(particle2), length(length), k(k) {
//     }
// };

class DihedralForce::DihedralInfo {
    public:
        int particle1, particle2, particle3, particle4;
        double  ka1, phia1,  epsd,  ka2,  phia2,  eps0,  kb1,  phib1, eps1,  kb2,  phib2,  eps2;
        DihedralInfo(){
                       particle1= particle2 =particle3=particle4 =-1;
                       ka1 = phia1 = epsd = ka2 = phia2 =  eps0 = kb1 =  phib1 = eps1 =  kb2 =  phib2 =  eps2 =0.0;
         }
        
        DihedralInfo(int particle1, int particle2, int particle3, int particle4,
                     double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2) :
            particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), ka1(ka1), phia1(phia1),  epsd(epsd),  ka2(ka2),  phia2(phia2),  eps0(eps0),  kb1(kb1),  phib1(phib1), eps1(eps1),  kb2(kb2),  phib2(phib2),  eps2(eps2) {}
    
};


} // namespace DihedralPlugin

#endif /*OPENMM_DIHEDRALFORCE_H_*/
