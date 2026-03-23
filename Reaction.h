#ifndef REACTION_H
#define REACTION_H
#include "Entite.h"
#include <vector>
#include <string>
#include <map>

using namespace std;

class Reaction {
    public: // On déclare les attributs et le constructeur
    vector<Entite*> reactifs;
    vector<Entite*> produits;
    
    double E_a;
    
    double kforward;
    double kbackward;
    
    Reaction(vector<Entite*> r, vector<Entite*> p, double E);
    
    double vitesse(bool, double V=1);
    
    double DeltaG();
    
};

#endif // REACTION_H



