#include "crpropa/module/PlasmaInstability.h"

PlasmaInstability::PlasmaInstability(double a, double Energyscale, 
                                     double Lengthscale)
                                      : crpropa::Module() {
        setPlasmaIndex(a);
        setPlasmaEnergyScale(Energyscale);
        setPlasmaLengthScale(Lengthscale);
        setDescription("PlasmaInstability");
}

void PlasmaInstability::setPlasmaIndex(double alpha) {
        index = alpha;
}

void PlasmaInstability::setPlasmaEnergyScale(double Escale) {
        tildeE = Escale;
}

void PlasmaInstability::setPlasmaLengthScale(double lambdascale) {
        lambda = lambdascale;
}

double PlasmaInstability::getPlasmaIndex() const {
        return index;
}

double PlasmaInstability::getPlasmaEnergyScale() const {
        return tildeE;
}

double PlasmaInstability::getPlasmaLengthScale() const {
        return lambda;
}

void PlasmaInstability::process(crpropa::Candidate *candidate) const {
        int id = candidate->current.getId();

        // only works for electrons and positrons
        if (fabs(id) != 11)
                return;

        double dx = candidate->getCurrentStep();
        double z = candidate->getRedshift();
        double E = candidate->current.getEnergy();
        double dEdx = coolingPower(E, z);

        if (dEdx < 0)
                dEdx = 0;
        double Enew = E - dEdx * dx;
        candidate->current.setEnergy(Enew);
        candidate->limitNextStep(0.1 * E / dEdx);
}

double PlasmaInstability::coolingPower(double E, double z) const {
        double dEdx = 0;
        E /= (1 + z);
        return (tildeE / lambda) * pow( E  / tildeE, 1-index);
}

std::string PlasmaInstability::getDescription() const {
        std::stringstream s;
        s << "PlasmaInstability";

        return s.str();
}