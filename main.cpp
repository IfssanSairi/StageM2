#include <iostream>
#include <sstream>
#include <fstream> // to deal with files
#include "Reaction.h"
#include "Entite.h"
#include <cmath>
#include <random>
#include <vector>
#include <map>
#include <string>
#include <list>
#include <getopt.h>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

using namespace std;


void trimWithBadCharacters(std::string & str)
{
    while (str.find('\r') != str.npos)
        str.erase( str.find('\r'));
}


bool isValidInteger(string& s)
{
    // better way : trim input character with spaces, new lines...
    if (s.find('\r') != s.npos)
    {
        int k = s.find('\r');
        s = s.erase(k, 1);
    }
    
    return( strspn( s.c_str(), "-0123456789" ) == s.size() );
}

struct Reseau // classe sans méthodes pour obtenir les entités et les réactions du réseau
{
    vector<Entite*> entites ; // pointeurs vers objets de type Entite
    vector<Reaction*> reactions; // pointeurs vers objets de type Reaction
    // Matrice stoechiométrique pour le code déterministe
    vector<vector<int>> M;
};


int findEntityByName(vector<Entite*> entites, string target_name) // fonction permettant de vérifier la correspondance entre entités de la partie matrice et de la partie entités
{
    int index = -1;
    for (size_t k =0; k < entites.size(); k++){
        if (entites[k]->name==target_name){ // entites est bien un vecteur qui contient des pointeurs
            index = k;
            break;
        }
    }
return index;
}

int findReactionByName(vector<Reaction*> reactions, string target_name)
{
    int index = -1;
    for (size_t k =0; k < reactions.size(); k++){
        cout << "candidate: " << reactions[k]->name << endl;
        if (reactions[k]->name.find('\r') != reactions[k]->name.npos)
            cout << "WARNING" << endl;
        if (reactions[k]->name==target_name){ // entites est bien un vecteur qui contient des pointeurs
            index = k;
            break;
        }
    }
return index;
}

void initialiseStoichiometricMatrix(Reseau & reseau, vector<string> lines) // reseau est passé en référence ici (modifier la valeur de reseau dans la fonction modifie aussi sa valeur dans la fonction principale)
{
        cout << "---initialiseMatrix---"<< endl;
        for (auto & el : lines)
            cout << el << endl;
    
    vector<vector<int>> matrix (lines.size() - 1, vector<int> (0)); // nb lignes = nb entités donc on skip le header

    vector<Reaction*> reac;
    vector<Entite*> ent;
    
    for (size_t j=0; j < lines.size();j++){
        vector<int> row;
        if (j==0){ // 1ère ligne avec les noms des réactions (R1, R2, R3...)
            string element;
            stringstream ss(lines[j]);
            while (getline(ss, element, ','))
            {
                if (!element.empty()) // si ce header n'est pas vide ?
                {
                    cout << "Init reaction with name : " << element << endl;
                    trimWithBadCharacters(element);
                    Reaction * newreac = new Reaction(element); // initialisation du constructeur // pointeur vers un objet Reaction construit avec le constructeur qui rentre seulement les noms
                    reac.push_back(newreac); // on ajoute ce pointeur dans le vecteur qui convient (aucune info sur la réaction pour l'instant)
                }
                    
            }
        } // end j == 0
        else // entities car à la ligne suivante on arrive sur les chaînes de caractères des entités
        {
            int compteur=0;
            string element;
            stringstream ss(lines[j]);
            while (getline(ss, element, ','))
            {
                if (compteur==0){ // colonne avec les noms des entités
                    trimWithBadCharacters(element);
                    Entite * newentity = new Entite(element); // constructeur
                    ent.push_back(newentity);
                    cout << "init entity with name " << element << endl;
                }
                else // autres colonnes
                {
                    if (!isValidInteger(element))
                        throw runtime_error("invalid integer entry.");
                    int s = atoi(element.c_str());
                    row.push_back(s);
                    size_t index_reac = compteur - 1; // si on a des nb stoechio qui dépassent les colonnes réactions
                    if (index_reac>=reac.size())
                        throw runtime_error("Stoichiometric index out of reaction vector.");
                    for (int k=0; k<abs(s); k++)
                    {
                        if (s>0)
                        {
                            reac[index_reac]->addProduct(ent.back());
                        }
                        else
                        {
                            reac[index_reac]->addReactant(ent.back());
                        }
                    }
                }
                compteur++;
            }
            matrix[j-1]=row; // on remplit la matrice des rows
        }
    }
    
    reseau.entites = ent;
    reseau.reactions = reac;
    reseau.M = matrix;

    //debugging
    cout << "\nNetwork init as:" << endl;
    for (auto & r : reac)
    {
        cout << r->name << ": ";
        for (auto & e : r->reactifs) // rappel : les réactifs sont des pointeurs vers les entités
            cout << " " << e->name; //
        cout << " --> ";
        for (auto & e : r->produits)
            cout << " " << e->name;
        cout << endl;
    }
    cout << endl;
    
}


void initialiseEntities(Reseau& reseau, vector<string> lines)
{
    cout << "---initialiseEntities---"<< endl;
    for (auto & el : lines)
        cout << el << endl;

    for (size_t j=1; j < lines.size(); j++){ // on commence à 1 pour exclure le header
        
        // read line that is comma separated
        string line = lines[j];
        stringstream ss(line);
        string name, number, free_energy, concentration_ext;

        // Read model and skip unwanted columns
        getline(ss, name, ',');
        getline(ss, number, ',');
        getline(ss, free_energy, ',');
        getline(ss, concentration_ext, ',');
        
        int entindex = findEntityByName(reseau.entites, name); // est-ce qu'il y a bien correspondance entre les entités de la matrice et les entités de la partie entités
        if (entindex<0)
        {
            cout << "Entity not found !! --> " << name << endl;
            throw runtime_error("Entity not found by name");
        }
        
        // Convert variables string to double
        double effectif = stod(number);
        double energie_libre = stod(free_energy);
        double con_ext = stod(concentration_ext);
        
        //On retrouve l'entité correspondante par son nom dans la partie Entités et on récupère ses attributs
        reseau.entites[entindex]->effectif = effectif;
        reseau.entites[entindex]->energie_libre = energie_libre;
        reseau.entites[entindex]->concentration_ext = con_ext;
    }

}

void initialiseReactions(Reseau& reseau, vector<string> lines)
{
    cout << "---initialiseReactions---"<< endl;
    //for (auto & el : lines)
    //    cout << el << endl;

    for (size_t j=1; j < lines.size(); j++){ // on commence à 1 pour exclure le header
        
        // read line that is comma separated
        string line = lines[j];
        stringstream ss(line);
        string name, Ea;

        // Read model and skip unwanted columns
        getline(ss, name, ',');
        getline(ss, Ea, ',');
        
        int entindex = findReactionByName(reseau.reactions, name);
        if (entindex<0)
        {
            cout << "Reaction not found !! --> " << name << endl;
            throw runtime_error("Reaction not found by name");
        }
        
        // Convert variables string to double
        double activation_energy = stod(Ea);
        
        //On retrouve la réaction correspondante par son nom dans la partie Reactions et on récupère ses attributs
        reseau.reactions[entindex]->E_a = activation_energy;
        cout << "calling updateReactionRates()" << endl;
        reseau.reactions[entindex]->updateReactionRates();
        
    }

}

Reseau initialiseReactionNetwork(string inputnetwork){
    
    Reseau reseau; // On crée un objet réseau, qui a comme attributs des vecteurs de pointeurs vers les entités et les réactions
    
    ifstream csv_file(inputnetwork);
    vector<string> lines;
    string line;
    
    // index flags for matrix, entities and reactions
    size_t flag_matrix = 0;
    size_t flag_reactions = 0;
    size_t flag_entites = 0;
    
    int c=-1;
    while (getline(csv_file, line)) {
        c++; // c est incrémenté en début de boucle, comme ça qu'on "avance" dans le fichier
        stringstream ss(line);
        
        if (line.find("MATRIX") != line.npos)
            flag_matrix = c;
        else if (line.find("REACTIONS") != line.npos)
            flag_reactions = c;
        else if (line.find("ENTITES") != line.npos)
            flag_entites = c;
        
        lines.push_back(line);
        
    }
    
    vector<string> lines_matrix, lines_entites, lines_reactions; // les lignes qui correspondent à chaque partie
    
    for (size_t i=0; i<lines.size(); i++)
        {
            if (i<flag_entites && i>flag_matrix)
                lines_matrix.push_back(lines[i]);
            else if (i>flag_entites && i<flag_reactions)
                lines_entites.push_back(lines[i]);
            else if (i>flag_reactions)
                lines_reactions.push_back(lines[i]);
        }
    // init stoichiometric matrix of the reaction network
    try{
        initialiseStoichiometricMatrix(reseau, lines_matrix);
    }
    catch( const runtime_error& error ){
        cout << error.what() << endl;
    }
    // init entities of the reaction network
    try{
        initialiseEntities(reseau, lines_entites);
    }
    catch( const runtime_error& error ){
        cout << error.what() << endl;
    }
    
    // init reactions of the reaction network
    try{
        initialiseReactions(reseau, lines_reactions);
    }
    catch( const runtime_error& error ){
        cout << error.what() << endl;
    }
    
    return reseau;
}



void printReactionNetwork(Reseau * reseau)
{
    cout << "---- Printing Reaction Network ----" << endl;

    // print reactions in the form 'X  Y --> Z'
    cout << "Reactions stoichiometry:" << endl;
    for (auto & r : reseau->reactions)
    {
        cout << r->name << ": ";
        for (auto & e : r->reactifs)
            cout << " " << e->name;
        cout << " --> ";
        for (auto & e : r->produits)
            cout << " " << e->name;
        cout << endl; 
    }
    cout << endl;


    cout << "Entity properties" << endl;
    for (auto & e : reseau->entites)
    {
        cout << e->name << ".\tNumber : " << e->effectif << "\t. free energy : " << e->energie_libre   << ".\tc_ext= " << e->concentration_ext << endl; 
    }

    cout << "\nReactions properties" << endl;
    for (auto & r : reseau->reactions)
    {
        cout << r->name << ".\tEa : " << r->E_a << ".\tk+ : " << r->kforward << ".\tk- : " << r-> kbackward <<  endl; // PROBLEME : k valent 0 donc aucune réaction ne se produit
    }

}


// Variables générales (déterministe et stochastique)

double V =300; // Volume total
double p_renouvelé = 0.001; // part de volume renouvelé à l'entrée et à la sortie du système

Reseau reseau; // là on définit juste le réseau, on le remplit pas

// On définit la fonction pour la résolution des équations différentielles

void cycle(const vector<double>& y , vector<double> &dydt, double t){
    vector<double> v (reseau.reactions.size(), 0);
    size_t i;
    for (i=0; i< reseau.reactions.size(); i++){
        v[i]=reseau.reactions[i]->vitesse(true,false, V, y, reseau.entites) - reseau.reactions[i]->vitesse(false, false, V, y, reseau.entites);
    }
    for (size_t k=0; k < reseau.entites.size(); k++){
        dydt[k]=0;
        for (size_t j=0; j < reseau.reactions.size(); j++){
            dydt[k]+= reseau.M[k][j] * v[j];
        }
        dydt[k] += reseau.entites[k]->concentration_ext*p_renouvelé - p_renouvelé*(y[k]);
    }
}

// Autocatalytic cycle AB

int main(int argc,char* argv[]) { // for arguments
    
    //Reseau reseau;
    bool verbose = true;
    
    const struct option longopts[] =
    {
        {"reseau",   required_argument,  0, 'r'},// on met le nom du fichier qui correspond au réseau en entrée en argument
        {"tmax",     required_argument,  0, 't'},
        {"gillespie",required_argument,  0, 'g'}, // on ajoute l'option gillespie
        {"verbose",  no_argument,        0, 'v'},
        {"help",     no_argument,        0, 'h'},
        {0,0,0,0},
    };
    
    int index;
    int iarg=0;
    
    double tmax = 1000; // par défaut c'est ce temps maximal
    bool Gillespie = true; // par défaut
    
    //turn off getopt error message
    opterr=1;
    
    while(iarg != -1)
    {
        iarg = getopt_long(argc, argv, "r:t:hv", longopts, &index);
        
        switch (iarg)
        {
            case 'h':
                //std::cout << "You hit help" << std::endl;
                // printHelp(); // possibulité de coder une fonction qui explique comment se servir du programme.
                break;
                
            case 'v':
                verbose = true;
                break;
                
            case 't':
                tmax = stod(optarg);
                break;
                
            case 'r':
                cout << "You hit reaction : " << optarg << endl; // si on a bien spécifié un réseau, on lance l'initialisation du réseau
                try
            {
                reseau = initialiseReactionNetwork(string(optarg)); // c'est vraiment ici qu'on remplit le réseau de réaction
            }
                catch( const runtime_error& error )
            {
                cout << error.what() << endl;
            }
                break;
            case 'g':
            {
                string val = optarg;
                
                if (val == "true")
                    Gillespie = true;
                else if (val == "false")
                    Gillespie = false;
                else
                    throw runtime_error("Invalid value for --gillespie");
                
                break;
            }
        }
    }
    
    
    if (verbose)
        printReactionNetwork(&reseau);
    
    // Definition des vecteurs y et dydt
    
    vector<double> y(reseau.entites.size());
    vector<double> dydt(reseau.entites.size());
    
    if (Gillespie){
        
        
        //double freq_fix_mut=0.0;
        int nRuns = 1;
        vector<int> runs;
        vector <double> temps;
        vector <vector<double>> etats;
        vector<vector<vector<double>>> all_etats;
        vector<vector<double>> all_temps;
        vector<double> x ; // vector x to keep entities concentrations
        vector <vector<double>> propensions;
        
        vector<Entite*> entites_init = reseau.entites;
        
        // Par défaut, rand() fixe la graine (même résultat à chaque simulation)
        // Donc srand pour ne pas fixer la graine
        
        srand(time(NULL));
        
        // Déclaration et définition des variables propres à l'algo Gillespie
        
        double atot, t, tau, somme;
        unsigned mu; // unsigned pour entiers non négatifs
        
        // Sauvegarde des conditions initiales (PAS des pointeurs !) pour les effectifs (seulement ce dont on a besoin)
        vector<double> effectifs_init;
        
        for (auto e : reseau.entites){
            effectifs_init.push_back(e->effectif);
        }
        
        for (int i=0;i < nRuns; i++){
            
            // Initialisation des vecteurs pour les prochains runs
            
            for (size_t i = 0; i < reseau.entites.size(); i++){
                reseau.entites[i]->effectif = effectifs_init[i];
            }
            
            x.clear();
            propensions.clear();
            etats.clear();
            temps.clear();
            
            vector <double> r = {0.0, 0.0}; // initialisation pour les nombres aléatoires
            
            // Définition des vecteurs pour stocker les résultats de la dynamique
            
            temps = {0.0}; // initialisation du vecteur temps
            
            for (const auto& e : reseau.entites){
                x.push_back(e->effectif/V);
            }
            
            etats = {x}; // vecteur double états pour stocker tous les x
            
            
            // Initialisation du vecteur des propensions de réactions du cycle
            int nb_events = 2*reseau.reactions.size() + 2*reseau.entites.size();
            vector <double> a (nb_events);
            
            // Algo Gillespie
            // Initialisation
            
            t = 0.0;
            
            while (t>=0 && t<tmax){
                
                // Calcul des propensions de réaction (vitesses des réactions)
                // On distingue sens direct et sens inverse car on ne peut pas avoir de propensions négatives
                size_t i;
                for (i=0; i < reseau.reactions.size(); i++){ // car on cherche chaque objet dans la liste reactions
                    //cout << i <<endl;
                    a[2*i]=reseau.reactions[i]->vitesse(true, true, V, y, reseau.entites); // On divise bien par le volume dans les propensions (voir méthode vitesse dans classe Reaction)
                    a[2*i+1]=reseau.reactions[i]->vitesse(false,true, V, y, reseau.entites);
                }
                
                // A la fin de cette boucle, i = reaction.size() - 1 cad ici 2
                
                // Calcul des propensions de création(ou entrée) et de destruction(ou sortie)
                
                for (size_t j=0; j < reseau.entites.size(); j++){
                    a[2*i+2*j]= reseau.entites[j]->concentration_ext*V*p_renouvelé; // on parcourt tout le tableau donc on doit repartir à partir de 2*i cad 2*(reactions.size()-1)
                    a[2*i+2*j+1]= V*p_renouvelé*max(reseau.entites[j]->effectif/V,0.0); // les Vtot s'annulent
                    
                }
                
                atot=0.0;
                for(unsigned m=0; m < a.size(); ++m){
                    atot += a[m]; // propension totale à réagir
                }
                
                // Génération de r1 et r2 dans une loi uniforme
                
                r.at(0)=(double)rand()/(double)RAND_MAX; // r1
                r.at(1)=(double)rand()/(double)RAND_MAX; // r2
                
                // Temps jusqu'à la prochaine simulation (distribution exponentielle) :
                
                tau = (-log(r[0])) / atot;
                t=t+tau;
                
                // Détermination de l'évènement :
                
                somme = 0.0;  // Initialisation de la somme
                
                for(unsigned j=0; j < a.size(); j++){
                    somme = somme + a[j]; // on fait la somme mais ce qui compte c'est bien la longueur de l'intervalle que prend la réaction sur le segment total
                    
                    //  && somme - a[j] < atot*r[1]
                    
                    // puisque r tiré dans loi uniforme (0,1), en multipliant par atot, on "créé" l'intervalle (0,atot)
                    if(somme >= atot*r[1]){ // Condition pour choisir la réaction
                        mu=j; // mu représente l'indice de la réaction
                        break;
                    } // on sort de la boucle dès qu'on a trouvé la réaction à se produire
                }
                // Une fois la réaction choisie, on change le vecteur x du nombre d'entités
                
                if (mu < 2*reseau.reactions.size()){ // Si on se trouve dans la partie du segment concernant les réactions
                    int dir=(mu%2==0)?1:-1; // mu pair, dir=1 sinon dir=-1
                    int idx = mu/2;
                    
                    for (const auto& reac : reseau.reactions[idx]->reactifs){ // Rappel : pointeurs vers objets Entite dans réactif
                        reac->effectif+=-dir; // si pair, sens direct, donc consommation des réactifs
                        
                    }
                    
                    for (const auto& prod : reseau.reactions[idx]->produits){
                        prod->effectif+=dir; // si impair, sens indirect, donc production des produits
                    }
                    
                }
                // Si on se trouve dans l'autre partie du segment concernant les créations et destructions
                // On commence donc à mu = 2*reactions.size()
                else {
                    int dir=(mu%2==0)?1:-1;// si pair, création sinon destruction
                    int idx = mu - 2*reseau.reactions.size(); // on reprend à 0 après les réactions
                    idx = idx/2; // on retrouve l'indice de l'entité
                    reseau.entites[idx]->effectif+=dir;
                }
                
                for (size_t i=0; i < reseau.entites.size(); i++){
                    x[i]=reseau.entites[i]->effectif/V;
                    
                }
                
                // Enregistrement de l'état :
                temps.push_back(t);
                etats.push_back(x);
                propensions.push_back(a);
            } // fin boucle while
            all_etats.push_back(etats);
            all_temps.push_back(temps);
            
        } // fin boucle Run
        
        int idx_AB = findEntityByName(reseau.entites, "AB");
        int idx_CB = findEntityByName(reseau.entites, "CB");
        int idx_DB = findEntityByName(reseau.entites, "DB");
        int idx_EB = findEntityByName(reseau.entites, "EB");
        int idx_FB = findEntityByName(reseau.entites, "FB");
        int idx_GB = findEntityByName(reseau.entites, "GB");
        int idx_HB = findEntityByName(reseau.entites, "HB");
        
        vector<int> indices = {idx_AB,idx_CB,idx_DB,idx_EB,idx_FB, idx_GB, idx_HB};
        
        map<vector<bool>, int> count_config;
        
        for (size_t run = 0; run < all_etats.size(); run++) {
            const auto& etat_final = all_etats[run].back(); // on prend le dernier vecteur du vecteur de tous les états pour un run donné
            vector<bool> config;
            config.reserve(indices.size()); // vecteur de booléens qui renseigne sur les états finaux
            
            for (int idx : indices)
            {
                config.push_back(etat_final[idx] > 0); // l'inégalité transforme la valeur en booléen
            }
            
            count_config[config]++;
            
        }
        
        // Création d'un fichier csv
        
        ofstream file("gillespie.csv");

        file << "Run,Temps";
        for (size_t i = 0; i < reseau.entites.size(); ++i)
            file << "," << reseau.entites[i]->name;

        file << '\n';
        
        //file << "Run,Temps,A,B,C,AB,ABA,ABAB,ABC,ABCB,CB,CBC,CBCB \n";
        //file << "Run,Temps,A,B,C,D,E,F,G,H,AB,ABA,ABAB,ABC,ABCB,CB,CBC,CBCB,ABD,ABDB,DB,DBD,DBDB,CBE,CBEB,EB,EBE,EBEB,CBF,CBFB,FB,FBF,FBFB,DBG,DBGB,GB,GBG,GBGB,DBH,DBHB,HB,HBH,HBHB \n"; // en-tête
        
        if (!file.is_open()) {
            cerr << "Failed to open file!" << endl;
            return 1;
        }
        
        //vector<vector<double>> data ;
        
        //for (size_t k=0; k< etats.size(); k++){
        //vector <double> row;
        //row.push_back(temps[k]);
        //for (double x : etats[k]){
        //row.push_back(x);
        //}
        //data.push_back(row);
        //}
        
        //for (const auto& row : data) {
        //for (size_t i = 0; i < row.size(); ++i) {
        //file << row[i];
        //if (i != row.size() - 1) file << ",";
        //}
        //file << "\n";
        //}
        
        for (size_t run = 0; run < all_etats.size(); run++){
            for (size_t k = 0; k < all_etats[run].size(); k++){
                file << run << "," << all_temps[run][k];
                for (double v : all_etats[run][k]){
                    file << "," << v;
                }
                file << "\n";
            }
        }
        
        file.close();
        cout << "CSV file created successfully." << endl;
        
        // Affichage Texte des résultats
        
        cout << "Temps : " << "Propensions : " << " Etats : \n ";
        
        
        for (size_t k=0; k < etats.size(); k++){
            cout << temps[k] << " ";
            
            cout << "{";
            for (double h : propensions[k]){
                
                cout << h << ",";
            }
            cout << "} ";
            
            cout << "{";
            
            for (double x : etats[k]){
                
                cout << x << ",";
            }
            cout << "} \n";
            
        }
        
        for (const auto& config : count_config) {
            cout << "Configuration (AB,CB,DB,EB,FB,GB,HB):";
            
            for (bool b : config.first) {
                cout << b << ",";
            }
            
            cout << " Occurrence: " << config.second << "\n";
        }
        
        
        //for (size_t i=0; i < config_count.size(); i++){
        //if (config_count[i]!=0){
        //cout << "Configuration n°:" << i << "\n";
        //cout << "Nombre de runs dans cette configuration:" << config_count[i] << "\n";
        //}
        // }
        
        
        //cout << "Fréquence de fixation de la mutation :" << (double)freq_fix_mut / nRuns << endl ;
        
        //if (runs.size()!=0){
        //for (size_t i=0; i < runs.size(); i++){
        //cout << runs[i] << ",";
        
        //}
        //cout << endl;
        //}
    } // fin de la condition if (Gillespie)
    // Simulation déterministe:
    else {
        
        for(size_t i=0;i<reseau.entites.size();i++){
            y[i] = reseau.entites[i]->effectif/V; // concentration
        }
        
        ofstream out("resultats.csv");
        out << "Temps";
        for (size_t i = 0; i < reseau.entites.size(); i++) out << "," << reseau.entites[i]->name;
        out << "\n";
        
        auto observer = [&](const vector<double>& y, double t) {
            out << t;
            for (double val : y) out << "," << val;
            out << "\n";
        };
        
        // Résolution avec la fonction integrate
        double dt = 0.1; // définition du pas de temps
        
        integrate_const(runge_kutta4<vector<double>>(),cycle,y,0.0,tmax,dt,observer);
        
        for (size_t k=0; k < reseau.entites.size(); k++){
            for (size_t j=0; j < reseau.reactions.size(); j++){
                cout << reseau.M[k][j] << " \n"[j == reseau.reactions.size()-1];
            }
        }
        return 0;
    }
}
