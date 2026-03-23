#ifndef ENTITE_H
#define ENTITE_H
#include <vector>
#include <string>
#include <map>

using namespace std;

class Entite {
    public :
    string name;
    double effectif;
    double energie_libre;
    double concentration_ext; // dépend du Vext, fixe pour chaque entité
    
    Entite(string n, double e, double e_l, double ce=0);

};

#endif // ENTITE_H



