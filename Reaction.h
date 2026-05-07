#ifndef REACTION_H
#define REACTION_H
#include "Entite.h"
#include <vector>
#include <string>
#include <map>

using namespace std;

class Reaction {
    public: // On déclare les attributs et le constructeur
    string name;
    vector<Entite*> reactifs;
    vector<Entite*> produits;
    
    double E_a;
    
    double kforward;
    double kbackward;
    
    Reaction(){};
    Reaction(string n);
    Reaction(string n, vector<Entite*> r, vector<Entite*> p, double E);
    
    void updateReactionRates();
    void addReactant(Entite*);
    void addProduct(Entite*);
    
    double vitesse(bool,bool, double V, const vector<double>& y, const vector<Entite*> entites);
    
    double DeltaG();
    
};

#endif // REACTION_H



