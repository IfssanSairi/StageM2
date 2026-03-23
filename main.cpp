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


using namespace std;


bool isValidInteger(const std::string& s)
{
    return( strspn( s.c_str(), "-0123456789" ) == s.size() );
}

struct Reseau
{
    vector<Entite*> entites ;
    vector<Reaction*> reactions;
};


int findEntityByName(vector<Entite*> entites, string target_name)
{
    int index = -1;
    for (size_t k =0; k < entites.size(); k++){
        if (entites[k]->name==target_name){
            index = k;
            break;
        }
    }
return index;
}


void initialiseStoichiometricMatrix(Reseau & reseau, vector<string> lines)
{
    
    vector<Reaction*> reac;
    vector<Entite*> ent;
    
    for (size_t j=0; j < lines.size();j++){
        if (j==0){ // save reactions
            string element;
            stringstream ss(lines[j]);
            while (getline(ss, element, ','))
            {
                if (!element.empty())
                {
                    Reaction * newreac = new Reaction(); // initialisation du constructeur
                    reac.push_back(newreac);
                }
                    
            }
        } // end j == 0
        else // entities
        {
            int compteur=0;
            string element;
            stringstream ss(lines[j]);
            while (getline(ss, element, ','))
            {
                if (compteur==0){
                    Entite * newentity = new Entite(element);
                    ent.push_back(newentity);
                }
                else
                {
                    if (!isValidInteger(element))
                        throw runtime_error("invalid integer entry.");
                    int s = atoi(element.c_str());
                    size_t index_reac = compteur - 1;
                    if (index_reac>=reac.size())
                        throw std::runtime_error("Stoichiometric index out of reaction vector.");
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
        }
    }
    
    reseau.entites = ent;
    reseau.reactions = reac;
    
}


void initialiseEntities(Reseau& reseau, vector<string> lines)
{

    for (size_t j=0; j < lines.size(); j++){
        
        // read line that is comma separated
        string line = lines[j];
        stringstream ss(line);
        string name, number, free_energy, concentration_ext;

        // Read model and skip unwanted columns
        getline(ss, name, ',');
        getline(ss, number, ',');
        getline(ss, free_energy, ',');
        getline(ss, concentration_ext, ',');
        
        int entindex = findEntityByName(reseau.entites, name);
        if (entindex<0)
        {
            cout << "Entity nout found !! --> " << name << endl;
            throw std::runtime_error("Entity not found by name");
        }
        
        // Convert variables string to double
        double effectif = stod(number);
        double energie_libre = stod(free_energy);
        double con_ext = stod(concentration_ext);
        
        //
        reseau.entites[entindex]->effectif = effectif;
        reseau.entites[entindex]->energie_libre = energie_libre;
        reseau.entites[entindex]->concentration_ext = con_ext;
    }
    

}

Reseau initialiseReactionNetwork(string inputnetwork){
    
    Reseau reseau;
    
    ifstream csv_file(inputnetwork);
    vector<string> lines;
    string line;
    
    // index flags for matrix, entities and reactions
    size_t flag_matrix = 0;
    size_t flag_reactions = 0;
    size_t flag_entites = 0;
    
    // Read rows with entites data
    int c=-1;
    while (getline(csv_file, line)) {
        c++;
        stringstream ss(line);
        
        if (line.find("MATRIX") != line.npos)
            flag_matrix = c;
        else if (line.find("REACTIONS") != line.npos)
            flag_reactions = c;
        else if (line.find("ENTITES") != line.npos)
            flag_entites = c;
        
        // Read model and skip unwanted columns
        // getline(ss, element, ',');
        lines.push_back(line);
        
        // Add entite to the vector
    }
    
    vector<string> lines_matrix, lines_entites, lines_reactions;
    for (size_t i=0; i<lines.size(); i++)
    {
        if (i<flag_entites && i>flag_matrix)
            lines_matrix.push_back(lines[i]);
        else if (i>flag_entites && i<flag_reactions)
            lines_entites.push_back(lines[i]);
        else
            lines_reactions.push_back(lines[i]);
    }
    
    // init stoichiometric matrix of the reaction network
    try{
        initialiseStoichiometricMatrix(reseau, lines_matrix);
    }
    catch( const std::runtime_error& error ){
        cout << error.what() << endl;
    }
    
    try{
        initialiseEntities(reseau, lines_entites);
    }
    catch( const std::runtime_error& error ){
        cout << error.what() << endl;
    }
    //initialiseReactions();
    
    //for (size_t k=0; k<lines.size(); k++)
     //   cout << lines[k] << endl;
    
   
    
    /*
     define printReactionNetwerk(reac, ent);
    cout << "added " << reac.size() << " reactions" << endl;
    for (size_t i=0; i< ent.size(); i++){
        cout << ent[i]->name << ",";
    }
    cout << endl;
    
    for (size_t i=0; i< reac.size(); i++){
        cout << "reac #" << i << endl;
        cout << "reactifs :" << endl;
        for (size_t j=0; j<reac[i]->reactifs.size(); j++)
            cout << "\t" << reac[i]->reactifs[j]->name << endl;
        cout << "produits :" << i << endl;
        for (size_t j=0; j<reac[i]->produits.size(); j++)
            cout << "\t" << reac[i]->produits[j]->name << endl;
    }
    cout << endl;
    */
    
    return reseau;
}

// Autocatalytic cycle AB

int main(int argc,char* argv[]) { // for arguments
    
    Reseau reseau;
    
    const struct option longopts[] =
      {
        {"reseau",   required_argument,   0, 'r'},
        {"help",      no_argument,        0, 'h'},
        {0,0,0,0},
      };

      int index;
      int iarg=0;

      //turn off getopt error message
      opterr=1;

      while(iarg != -1)
      {
        iarg = getopt_long(argc, argv, "r:h", longopts, &index);

        switch (iarg)
        {
          case 'h':
            std::cout << "You hit help" << std::endl;
            break;

          case 'r':
            std::cout << "You hit reaction : " << optarg << std::endl;
            try
            {
                reseau = initialiseReactionNetwork(string(optarg));
            }
            catch( const std::runtime_error& error )
            {
                cout << error.what() << endl;
            }
            break;

        }
      }

  
    
    
    // File_entites to read
    vector<string> file_entites = {"entites - Feuille 1.csv"};

    // List to store all entites data
    vector<Entite> entites;


    for (const auto& filename : file_entites) {
    ifstream csv_file(filename);
    string line;

    // Skip the header line
    getline(csv_file, line);

    // Read rows with entites data
    while (getline(csv_file, line)) {
        stringstream ss(line);
        string name, number, free_energy, concentration_ext;

        // Read model and skip unwanted columns
        getline(ss, name, ',');
        getline(ss, number, ',');
        getline(ss, free_energy, ',');
        getline(ss, concentration_ext, ',');

        // Convert variables string to double
        double effectif = stod(number);
        double energie_libre = stod(free_energy);
        double con_ext = stod(concentration_ext);

        // Add entite to the vector
        entites.push_back(Entite (name, effectif, energie_libre, con_ext));
    }
}

    //double freq_fix_mut=0.0;
    int nRuns = 1;
    vector<int> runs;
    vector <double> temps;
    vector <vector<double>> etats;
    vector<vector<vector<double>>> all_etats;
    vector<vector<double>> all_temps;
    vector<double> x ; // vector x to keep entities concentrations
    vector <vector<double>> propensions;
    
    vector<Entite> entites_init = entites;
    
    // Par défaut, rand() fixe la graine (même résultat à chaque simulation)
    // Donc srand pour ne pas fixer la graine
    
    srand(time(NULL));
    
    // Déclaration et définition des variables propres à l'algo Gillespie
    
    double V =300; // Volume total
    double p_renouvelé = 0.001; // part de volume renouvelé à l'entrée et à la sortie du système
    
    
    double tmax = 1000; // par défaut c'est ce temps maximal
    if (argc > 1){
        tmax=stod(argv[1]);
        
    }
    
    double atot, t, tau, somme;
    unsigned mu; // unsigned pour entiers non négatifs
        
        for (int i=0;i < nRuns; i++){
            
            // Initialisation des vecteurs pour les prochains runs
            
            entites = entites_init;
            
            x.clear();
            propensions.clear();
            etats.clear();
            
            vector <double> r = {0.0, 0.0}; // initialisation pour les nombres aléatoires
            
            //declarer des pointeurs vers les entites
            Entite* A=&entites[0];
            Entite* B=&entites[1];
            Entite* C=&entites[2];
            Entite* AB=&entites[3];
            Entite* ABA=&entites[4];
            Entite* ABAB=&entites[5];
            Entite* ABC=&entites[6];
            Entite* ABCB=&entites[7];
            Entite* CB=&entites[8];
            Entite* CBC=&entites[9];
            Entite* CBCB=&entites[10];
            
            // Création des objets de la classe Réaction
            // Objets stockés dans un vecteur
            
            // Reaction AB+A=ABA
            // Reaction ABA+B=ABAB
            // Reaction ABAB=2AB
            
            vector<Reaction> reactions = {
                Reaction({A,AB}, {ABA}, 0.0),
                Reaction({ABA, B}, {ABAB}, 0.0),
                Reaction({ABAB}, {AB, AB}, 0.0),
                Reaction({AB, C}, {ABC}, 5.0), // réactions de mutation
                Reaction({ABC, B}, {ABCB}, 0.0),
                Reaction({ABCB}, {AB, CB}, 0.0),
                Reaction({CB, C}, {CBC}, 0.0),// réactions du cycle mutant
                Reaction({CBC,B}, {CBCB}, 0.0),
                Reaction({CBCB}, {CB,CB}, 0.0)
            };
            
            // Définition des vecteurs pour stocker les résultats de la dynamique
            
            temps = {0.0}; // initialisation du vecteur temps
            
            for (const auto& e : entites){
                x.push_back(e.effectif/V);
            }
            
            etats = {x}; // vecteur double états pour stocker tous les x
            
            
            // Initialisation du vecteur des propensions de réactions du cycle
            int nb_events = 2*reactions.size() + 2*entites.size();
            vector <double> a (nb_events);
            
            // Algo Gillespie
            // Initialisation
            
            t = 0.0;
            
            while (t>=0 && t<tmax){
                
                // Calcul des propensions de réaction (vitesses des réactions)
                // On distingue sens direct et sens inverse car on ne peut pas avoir de propensions négatives
                size_t i;
                for (i=0; i < reactions.size(); i++){ // car on cherche chaque objet dans la liste reactions
                    //cout << i <<endl;
                    a[2*i]=reactions[i].vitesse(true, V); // On divise bien par le volume dans les propensions (voir méthode vitesse dans classe Reaction)
                    a[2*i+1]=reactions[i].vitesse(false, V);
                }
                
                // A la fin de cette boucle, i = reaction.size() - 1 cad ici 2
                
                // Calcul des propensions de création(ou entrée) et de destruction(ou sortie)
                
                for (size_t j=0; j < entites.size(); j++){
                    a[2*i+2*j]= entites[j].concentration_ext*V*p_renouvelé; // on parcourt tout le tableau donc on doit repartir à partir de 2*i cad 2*(reactions.size()-1)
                    a[2*i+2*j+1]= V*p_renouvelé*entites[j].effectif/V; // les Vtot s'annulent
                    
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
                
                if (mu < 2*reactions.size()){ // Si on se trouve dans la partie du segment concernant les réactions
                    int dir=(mu%2==0)?1:-1; // mu pair, dir=1 sinon dir=-1
                    int idx = mu/2;
                    
                    for (const auto& reac : reactions[idx].reactifs){ // Rappel : pointeurs vers objets Entite dans réactif
                        reac->effectif+=-dir; // si pair, sens direct, donc consommation des réactifs
                        
                    }
                    
                    for (const auto& prod : reactions[idx].produits){
                        prod->effectif+=dir; // si impair, sens indirect, donc production des produits
                    }
                    
                }
                // Si on se trouve dans l'autre partie du segment concernant les créations et destructions
                // On commence donc à mu = 2*reactions.size()
                else {
                    int dir=(mu%2==0)?1:-1;// si pair, création sinon destruction
                    int idx = mu - 2*reactions.size(); // on reprend à 0 après les réactions
                    idx = idx/2; // on retrouve l'indice de l'entité
                    entites[idx].effectif+=dir;
                }
                
                for (size_t i=0; i < entites.size(); i++){
                    x[i]=entites[i].effectif/V;
                    
                }
                
                // Enregistrement de l'état :
                temps.push_back(t);
                etats.push_back(x);
                propensions.push_back(a);
            } // fin boucle while
            all_etats.push_back(etats);
            all_temps.push_back(temps);
            
        } // fin boucle Run
    
    //for (size_t run = 0; run < all_etats.size(); run++) {
        //const auto& etat_final = all_etats[run].back(); // on prend le dernier vecteur du vecteur de tous les états pour un run donné
        //if (etat_final[3] ==0 && etat_final[8] !=0) {
            //freq_fix_mut++;
            //temps_fix+= all_temps[run][];
        //}
    //}

    // Création d'un fichier csv
    
    ofstream file("gillespie.csv");

    
    file << "Run,Temps,A,B,C,AB,ABA,ABAB,ABC,ABCB,CB,CBC,CBCB \n"; // en-tête
    
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
    //cout << "Fréquence de fixation de la mutation :" << (double)freq_fix_mut / nRuns << endl ;
    //cout << "Temps moyen de fixation sur les runs : "<< (double)temps_fix / (freq_fix_mut) << endl ;
    
    //if (runs.size()!=0){
        //for (size_t i=0; i < runs.size(); i++){
            //cout << runs[i] << ",";
            
        //}
        //cout << endl;
    //}
    return 0;
}
