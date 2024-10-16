#ifndef CRPROPA_PLASMAINSTABILITY_H
#define CRPROPA_PLASMAINSTABILITY_H


#include <crpropa/Cosmology.h>
#include <crpropa/Module.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Units.h>


class PlasmaInstability : public crpropa::Module {
        protected:
                double index;
                double tildeE;
                double lambda;

        public:
                PlasmaInstability(double index, double tildeE, 
                double lambda);
                void setPlasmaIndex(double indexplasma);
                void setPlasmaEnergyScale(double tildeEplasma);
                void setPlasmaLengthScale(double lambdaplasma);
                double getPlasmaIndex() const;
                double getPlasmaEnergyScale() const;
                double getPlasmaLengthScale() const;
                double coolingPower(double energy, double redshift) 
                const;
                std::string getDescription() const;
                void process(crpropa::Candidate *candidate) const;
};


#endif