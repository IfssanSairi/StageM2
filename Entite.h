#ifndef ENTITE_H
#define ENTITE_H
#include <vector>
#include <string>
#include <map>

using namespace std;

class Entite {
    public :
    string name = "";
    double effectif=0;
    double energie_libre=0;
    double concentration_ext=0; // dépend du Vext, fixe pour chaque entité
    Entite (){};
    
    Entite(string n);
    
    Entite(string n, double e, double e_l, double ce=0);

};

#endif // ENTITE_H



