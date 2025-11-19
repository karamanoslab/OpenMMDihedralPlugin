#ifndef OPENMM_DIHEDRALFORCE_H_
#define OPENMM_DIHEDRALFORCE_H_


#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "windowsExportExample.h"
#include "openmm/reference/RealVec.h"


namespace DihedralPlugin {



class OPENMM_EXPORT_EXAMPLE DihedralForce : public OpenMM::Force {
public:

    DihedralForce();
     
    int getNumDihedrals() const {
         return dihedrals.size(); 
         }
    /**

     */
    int addDihedral(int p1, int p2, int p3, int p4,
                     double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2);
    /**
     */
    void getDihedralParameters(int index, int& p1, int& p2, int& p3, int& p4,
                               double& ka1, double& phia1, double& epsd, double& ka2, double& phia2, double& eps0, double& kb1, double& phib1, double& eps1, double& kb2, double& phib2, double& eps2) const;
    /**
     
     */
   
    void setDihedralParameters(int index, int p1, int p2, int p3, int p4,
                               double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2);
    /**
     */ 
     
    double calcDihedralAngleFromIndexes(const OpenMM::Context& context, int& p1, int& p2, int& p3, int& p4);
    
    double calcDihedralAngleFromCoords(const OpenMM::RealVec& pos1,
                                   const OpenMM::RealVec& pos2,
                                   const OpenMM::RealVec& pos3,
                                   const OpenMM::RealVec& pos4);

    void updateParametersInContext(OpenMM::Context& context);
    /**
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
