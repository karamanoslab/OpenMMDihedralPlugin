/* ------ */

#include "ReferenceDihedralKernels.h"
#include "DihedralForce.h"
#include "math.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/ReferenceForce.h"

using namespace DihedralPlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

inline double dot(const RealVec& a, const RealVec& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline OpenMM::Vec3 cross(const OpenMM::Vec3& a, const OpenMM::Vec3& b) {
    return OpenMM::Vec3(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return data->periodicBoxVectors;
}

void ReferenceCalcDihedralForceKernel::initialize(const System& system, const DihedralForce& force) {
    // Initialize dihedral parameters.
    
    int numDihedrals = force.getNumDihedrals();
    particle1.resize(numDihedrals);
    particle2.resize(numDihedrals);
    particle3.resize(numDihedrals);
    particle4.resize(numDihedrals);
  
    
    ka1.resize(numDihedrals);
    phia1.resize(numDihedrals);
    epsd.resize(numDihedrals);
    ka2.resize(numDihedrals);
    phia2.resize(numDihedrals);
    eps0.resize(numDihedrals);
    kb1.resize(numDihedrals);
    phib1.resize(numDihedrals);
    eps1.resize(numDihedrals);
    kb2.resize(numDihedrals);
    phib2.resize(numDihedrals);
    eps2.resize(numDihedrals);
    
    periodic = force.usesPeriodicBoundaryConditions();
    
    for (int i = 0; i < numDihedrals; i++)
        force.getDihedralParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], ka1[i], phia1[i],  epsd[i],  ka2[i],  phia2[i],  eps0[i],  kb1[i],  phib1[i], eps1[i],  kb2[i],  phib2[i],  eps2[i]  );
}



double ReferenceCalcDihedralForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
	vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    Vec3* boxVectors = extractBoxVectors(context);
    
    int numDihedrals = particle1.size();
    double energy = 0;
    const double TOLERANCE = 0.05;
    const double SMALL = 0.001;
    const double SMALLER = 0.00001;
    const double PI = M_PI;

	double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
	Vec3 f1, f2, f3, f4;
	double sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
	double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
	double c2mag,sc1,sc2,s1,s12,c,pp,ppd,a,a11,a22;
	double a33,a12,a13,a23,sx2,sy2,sz2;
	double s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2;
	double dphia,dphib,dphib2,dphic,dphic2,dphid,dphid2;
	double pa,pb,pb2,pc,pc2,pd,pd2,ea,eb,eb2,ec,ec2,ed,ed2;
	double fea,feb,feb2,fec,fec2,fed,fed2;
    const double kcal_to_kj = 4.184;
	   
	for (int i = 0; i < numDihedrals; i++) {
		int p1 = particle1[i];
		int p2 = particle2[i];
		int p3 = particle3[i];
		int p4 = particle4[i];
		
		double ka = ka1[i], kb = ka2[i];
		double kc = kb1[i], kd = kb2[i];
		double fa = phia1[i], fb = phia2[i];
		double fc = phib1[i], fd = phib2[i];
		double eps0_ = eps0[i], eps1_ = eps1[i], eps2_ = eps2[i], epsd_ = epsd[i];
		
		epsd_ *= kcal_to_kj;
        eps0_ *= kcal_to_kj;
        eps1_ *= kcal_to_kj;
        eps2_ *= kcal_to_kj;

 		double vb1[ReferenceForce::LastDeltaRIndex];
 		double vb2[ReferenceForce::LastDeltaRIndex];
		double vb3[ReferenceForce::LastDeltaRIndex];
		
 		if (periodic) {
 		    ReferenceForce::getDeltaRPeriodic(pos[p2], pos[p1], boxVectors, vb1);
		    ReferenceForce::getDeltaRPeriodic(pos[p2], pos[p3], boxVectors, vb2);
		    ReferenceForce::getDeltaRPeriodic(pos[p3], pos[p4], boxVectors, vb3); 
        } else {
             ReferenceForce::getDeltaR(pos[p2], pos[p1], vb1);
		    ReferenceForce::getDeltaR(pos[p2], pos[p3], vb2);
 		    ReferenceForce::getDeltaR(pos[p3], pos[p4], vb3);
 		}    
	  
	  
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
	
		// error check
		if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
			char str[128];
			printf("Dihedral Warning: atoms: %i %i %i %i cos_phi= %g\n",
					p1, p2, p3, p4, c); 
		}
	
		if (c > 1.0) c = 1.0;
		if (c < -1.0) c = -1.0;
	
		// force & energy
		// p = k ( phi- phi0)^2
		// pd = dp/dc
	
		phi = acos(c);
		if (dx > 0.0) phi *= -1.0;
		si = sin(phi);
		if (fabs(si) < SMALLER) si = SMALLER;
		siinv = 1.0/si;
	
		dphia = phi - fa;
		dphib = phi - fb;
		dphib2 = dphib + 2.0*M_PI;
		dphic = phi - fc;
		dphic2 = dphic - 2.0*M_PI;
		dphid = phi - fd;
		dphid2 = dphid - 2.0*M_PI;
	
		pa = -ka * dphia;
		pb = -kb * dphib * dphib * dphib;;
		pb2 = -kb * dphib2 * dphib2 * dphib2;
		pc = -kc * dphic;
		pc2 = -kc * dphic2;
		pd = -kd * dphid * dphid * dphid;
		pd2 = -kd * dphid2 * dphid2 * dphid2;
	
		ea = pa * dphia - epsd_;
		eb = pb * dphib + eps0_;
		eb2 = pb2 * dphib2 + eps0_;
		ec = pc * dphic + epsd_ + eps1_;
		ec2 = pc2 * dphic2 + epsd_ + eps1_; 
		ed = pd * dphid + eps2_ + eps1_;
		ed2 = pd2 * dphid2 + eps2_ + eps1_;
		
		fea = exp(ea);
		feb = exp(eb);
		feb2 = exp(eb2);
		fec = exp(ec);
		fec2 = exp(ec2);
		fed = exp(ed);
		fed2 = exp(ed2);
	
		pp = fea + feb + feb2 + fec +fec2 + fed + fed2;
		ppd = 2.0*pa*fea + 4.0*pb*feb + 4.0*pb2*feb2 + 2.0*pc*fec
			+ 2.0*pc2*fec2 + 4.0*pd*fed + 4.0*pd2*fed2;
	
		ppd /= pp;
		ppd *= siinv;
	
		energy += -log(pp);
		
		a = ppd;
		c = c * a;
		s12 = s12 * a;
		a11 = c*sb1*s1;
		a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2));
		a33 = c*sb3*s2;
		a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12);
		a13 = -rb1*rb3*s12;
		a23 = r12c2 * (c2mag*c*s2 + c1mag*s12);
	
		sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
		sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
		sz2  = a12*vb1z + a22*vb2z + a23*vb3z;
	
		f1[0] = a11*vb1x + a12*vb2x + a13*vb3x;
		f1[1] = a11*vb1y + a12*vb2y + a13*vb3y;
		f1[2] = a11*vb1z + a12*vb2z + a13*vb3z;
	
		f2[0] = -sx2 - f1[0];
		f2[1] = -sy2 - f1[1];
		f2[2] = -sz2 - f1[2];
	
		f4[0] = a13*vb1x + a23*vb2x + a33*vb3x;
		f4[1] = a13*vb1y + a23*vb2y + a33*vb3y;
		f4[2] = a13*vb1z + a23*vb2z + a33*vb3z;
	
		f3[0] = sx2 - f4[0];
		f3[1] = sy2 - f4[1];
		f3[2] = sz2 - f4[2];
	
		// Apply forces
		force[p1] += f1;
		force[p2] += f2;
		force[p3] += f3;
		force[p4] += f4;
		
	}

    return energy;
}



void ReferenceCalcDihedralForceKernel::copyParametersToContext(ContextImpl& context, const DihedralForce& force) {
    if (force.getNumDihedrals() != particle1.size())
        throw OpenMMException("updateParametersInContext: The number of Dihedrals has changed");
    for (int i = 0; i < force.getNumDihedrals(); i++) {
        int p1, p2, p3, p4;
        force.getDihedralParameters(i, p1, p2, p3, p4, ka1[i], phia1[i],  epsd[i],  ka2[i],  phia2[i],  eps0[i],  kb1[i],  phib1[i], eps1[i],  kb2[i],  phib2[i],  eps2[i]); 
        if (p1 != particle1[i] || p2 != particle2[i] || p3 != particle3[i] || p4 != particle4[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}
