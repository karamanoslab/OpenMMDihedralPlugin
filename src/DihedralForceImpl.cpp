#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "DihedralForceImpl.h"
#include "DihedralKernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>

using namespace DihedralPlugin;
using namespace OpenMM;
using namespace std;

DihedralForceImpl::DihedralForceImpl(const DihedralForce& owner) : owner(owner) {
}

DihedralForceImpl::~DihedralForceImpl() {
}

void DihedralForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcDihedralForceKernel::Name(), context);
    kernel.getAs<CalcDihedralForceKernel>().initialize(context.getSystem(), owner);
}

double DihedralForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcDihedralForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> DihedralForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcDihedralForceKernel::Name());
    return names;
}

vector<tuple<int, int, int, int> > DihedralForceImpl::getDihedralParticles() const {
    int numDihedrals = owner.getNumDihedrals();
    vector<tuple<int, int, int, int> > dihedrals(numDihedrals);
    for (int i = 0; i < numDihedrals; i++) {
        int p1, p2, p3, p4;
        double  ka1, phia1,  epsd,  ka2,  phia2,  eps0,  kb1,  phib1, eps1,  kb2,  phib2,  eps2;
        owner.getDihedralParameters(i, p1, p2, p3, p4, ka1, phia1,  epsd,  ka2,  phia2,  eps0,  kb1,  phib1, eps1,  kb2,  phib2,  eps2);
    }
    return dihedrals;
}

void DihedralForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcDihedralForceKernel>().copyParametersToContext(context, owner);
}
