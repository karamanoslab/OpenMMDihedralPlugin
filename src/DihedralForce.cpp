#include "DihedralForce.h"
#include "DihedralForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <string>
#include "openmm/internal/ContextImpl.h"



using namespace DihedralPlugin;
using namespace OpenMM;
using namespace std;

bool usePeriodic;
DihedralForce::DihedralForce() : usePeriodic(false) {
}

int DihedralForce::addDihedral(int particle1, int particle2, int particle3, int particle4,  double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2) {
    dihedrals.push_back(DihedralInfo(particle1, particle2, particle3, particle4,  ka1, phia1,  epsd,  ka2,  phia2,  eps0,  kb1,  phib1, eps1,  kb2,  phib2,  eps2));
    return dihedrals.size()-1;
}

void DihedralForce::getDihedralParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& ka1, double& phia1, double& epsd, double& ka2, double& phia2, double& eps0, double& kb1, double& phib1, double& eps1, double& kb2, double& phib2, double& eps2) const {
    ASSERT_VALID_INDEX(index, dihedrals);
    particle1 = dihedrals[index].particle1;
    particle2 = dihedrals[index].particle2;
    particle3 = dihedrals[index].particle3;
    particle4 = dihedrals[index].particle4;
   
    ka1= dihedrals[index].ka1;
    phia1= dihedrals[index].phia1;
    epsd= dihedrals[index].epsd;
    ka2= dihedrals[index].ka2;
    phia2= dihedrals[index].phia2;
    eps0= dihedrals[index].eps0;
    kb1= dihedrals[index].kb1;
    phib1= dihedrals[index].phib1;
    
    eps1= dihedrals[index].eps1;
    kb2= dihedrals[index].kb2;
    phib2= dihedrals[index].phib2;
    eps2= dihedrals[index].eps2;
}

void DihedralForce::setDihedralParameters(int index, int particle1, int particle2, int particle3, int particle4,  double ka1, double phia1, double epsd, double ka2, double phia2, double eps0, double kb1, double phib1, double eps1, double kb2, double phib2, double eps2) {
    ASSERT_VALID_INDEX(index, dihedrals);
    dihedrals[index].particle1 = particle1;
    dihedrals[index].particle2 = particle2;
    dihedrals[index].particle3 = particle3;
    dihedrals[index].particle4 = particle4;
        
    dihedrals[index].ka1 = ka1;
    dihedrals[index].phia1 = phia1;
    dihedrals[index].epsd = epsd;
    dihedrals[index].ka2 = ka2;
    dihedrals[index].phia2 = phia2;
    dihedrals[index].eps0 = eps0;
    dihedrals[index].kb1 = kb1;
    dihedrals[index].phib1 = phib1;
    
    dihedrals[index].eps1 = eps1;
    dihedrals[index].kb2 = kb2;
    dihedrals[index].phib2 = phib2;
    dihedrals[index].eps2 = eps2;
    
}


double dot(const RealVec& a, const RealVec& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

RealVec cross(const RealVec& a, const RealVec& b) {
    return RealVec(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}



double DihedralForce::calcDihedralAngleFromIndexes(const OpenMM::Context& context,
                                        int& p1, int& p2, int& p3, int& p4) {
    const auto& positions = context.getState(OpenMM::State::Positions).getPositions();
   
    return calcDihedralAngleFromCoords(positions[p1], positions[p2], positions[p3], positions[p4]);
}


double DihedralForce::calcDihedralAngleFromCoords(const OpenMM::RealVec& pos1,
                                   const OpenMM::RealVec& pos2,
                                   const OpenMM::RealVec& pos3,
                                   const OpenMM::RealVec& pos4) {

//   Input atomic positions

    const double TOLERANCE = 0.05;
    const double SMALL = 0.001;
    const double SMALLER = 0.00001;
    const double PI = M_PI;

	double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
	double sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
	double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
	double c2mag,sc1,sc2,s1,s12,c,pp,ppd,a,a11,a22;
	double s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2;
	
	 
	Vec3 vb1 = pos1 - pos2;
	Vec3 vb2 = pos3 - pos2;
	Vec3 vb3 = pos4 - pos3;
  
	vb1x=vb1[0];
	vb1y=vb1[1];
	vb1z=vb1[2];
	
	vb2x=vb2[0];
	vb2y=vb2[1];
	vb2z=vb2[2];
	
	vb2xm = -vb2x;
	vb2ym = -vb2y;
	vb2zm = -vb2z;
	
	vb3x=vb3[0];
	vb3y=vb3[1];
	vb3z=vb3[2];
	
	sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
	sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
	sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

	rb1 = sqrt(sb1);
	rb3 = sqrt(sb3);

	c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

	// 1st and 2nd angle
	b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
	b1mag = sqrt(b1mag2);
	b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
	b2mag = sqrt(b2mag2);
	b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
	b3mag = sqrt(b3mag2);

	ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
	r12c1 = 1.0 / (b1mag*b2mag);
	c1mag = ctmp * r12c1;

	ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
	r12c2 = 1.0 / (b2mag*b3mag);
	c2mag = ctmp * r12c2;

	// cos and sin of 2 angles and final c
	sin2 = max(1.0 - c1mag*c1mag,0.0);
	sc1 = sqrt(sin2);
	if (sc1 < SMALL) sc1 = SMALL;
	sc1 = 1.0/sc1;

	sin2 = max(1.0 - c2mag*c2mag,0.0);
	sc2 = sqrt(sin2);
	if (sc2 < SMALL) sc2 = SMALL;
	sc2 = 1.0/sc2;

	s1 = sc1 * sc1;
	s2 = sc2 * sc2;
	s12 = sc1 * sc2;
	c = (c0 + c1mag*c2mag) * s12;

	cx = vb1y*vb2z - vb1z*vb2y;
	cy = vb1z*vb2x - vb1x*vb2z;
	cz = vb1x*vb2y - vb1y*vb2x;
	cmag = sqrt(cx*cx + cy*cy + cz*cz);
	dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

	if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
		char str[128];
		printf(" Dihedral Warning: cos_phi= %g\n",
				 c); 
	}

	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;

	phi = acos(c);
	if (dx > 0.0) phi *= -1.0;
	
	return phi;
    
}


void DihedralForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool DihedralForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}

ForceImpl* DihedralForce::createImpl() const {
    return new DihedralForceImpl(*this);
}

void DihedralForce::updateParametersInContext(Context& context) {
    dynamic_cast<DihedralForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}


