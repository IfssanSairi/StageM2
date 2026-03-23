#include "Reaction.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>

// Définition des constructeurs

Reaction :: Reaction (vector<Entite*> r, vector<Entite*> p, double E){
    
    reactifs = r;
    produits = p;
    E_a = E;
    
    //kforward et kbackward sont bien des attributs, mais non explicités en arguments dans le constructeur
    
    if (DeltaG() <= 0){
        kforward = exp(-E_a);
        kbackward = exp(-E_a - abs(DeltaG())) ;
    }
    else {
        kforward = exp(-E_a - abs(DeltaG()));
        kbackward = exp(-E_a);
    }
}


double Reaction::vitesse(bool isForward, double V){
    double v = isForward?kforward:kbackward;
    Entite* lastreac = NULL;
    int compteur =0;
    vector<Entite *> *vecEnt=isForward?&reactifs:&produits; // Si réaction forward, vecteur *vecEnt devient le vecteur des réactifs, sinon dans celui des produits
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
        }
    return v*V;
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










