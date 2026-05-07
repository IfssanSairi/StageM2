#include "Reaction.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

// Définition des constructeurs

Reaction :: Reaction (string n, vector<Entite*> r, vector<Entite*> p, double E){
    name = n;
    reactifs = r;
    produits = p;
    E_a = E;
    
    //kforward et kbackward sont bien des attributs, mais non explicités en arguments dans le constructeur
    
    kforward =0; // valeur par défaut
    kbackward =0; // valeur par défaut
}

Reaction :: Reaction(string n){
    name =n;
}

void Reaction::updateReactionRates() {
    cout << "Reaction::updateReactionRates()" << endl;
    cout << "k+ k- = " << kforward << " " << kbackward << endl;
    double deltaG = DeltaG();

    if (deltaG <= 0) {
        kforward  = exp(-E_a);
        kbackward = exp(-E_a - abs(deltaG));
    } else {
        kforward  = exp(-E_a - abs(deltaG));
        kbackward = exp(-E_a);
    }
    cout << "bis k+ k- = " << kforward << " " << kbackward << endl;

}

void Reaction::addReactant(Entite* e){
    reactifs.push_back(e);
}

void Reaction::addProduct(Entite* e){
    produits.push_back(e);
    
}


double Reaction::vitesse(bool isForward, bool Gillespie, double V, const vector<double>& y, const vector<Entite*> entites){
    double v = isForward?kforward:kbackward;
    vector<Entite *> *vecEnt=isForward?&reactifs:&produits; // Si réaction forward, vecteur *vecEnt devient le vecteur des réactifs, sinon dans celui des produits
    
    Entite* lastreac = NULL;
    int compteur =0;
    
    if (Gillespie==true){
        // pourquoi *vecEnt et pas simplement vecEnt
        for (const auto& er : *vecEnt){
                
                if (er==lastreac){ // we assume that identical reactants are consecutive
                    compteur++;}
                else {
                    compteur=0;
                }
                
                if (er-> effectif <=0){
                    v*=0;
                }
                else {
                    v*= (er->effectif - compteur)/V;
                }
            lastreac=er; // on met à jour lastreac
        }
        
        return v*V;
    }
    else {
        for (const auto& er : *vecEnt){
            auto it = find_if(entites.begin(), entites.end(),
                              [&](const Entite* e){ return e == er; });
            size_t j = it - entites.begin();
            v *= y[j];
        }
        return v;
    }
}
    

double Reaction :: DeltaG(){
    double G_initial=0;
    double G_final=0;
    for (const auto& er : reactifs){
        G_initial += er->energie_libre; // déférencer le pointeur
        
    }
    for (const auto& er : produits){
        G_final += er->energie_libre;
        
    }
    
    return (G_final - G_initial);
}










