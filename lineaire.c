#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<gmp.h>
#include"lineaire.h"
#include "macro.h"


uint reduit(long** Registres_b, long* Set, uint ligne, uint ligne_l)
//Prend en entrée une matrice Registres_b, chaques valeurs étant un bit rangé dans des registres.
//Prend aussi un vecteur de longeur égale à la largeur de la matrice, initialisé à 0.
//
//Reduit la matrice de sort que:
//                              - La matrice est triangulaire supérieure.
//                              - Pour une colonne donnée, si c'est un 1 sur la diagonale, la colonne est nulle partout ailleurs.
//                              - Pour une ligne donnée, si c'est un 0 sur la diagonale, la ligne est nulle partout.
//
//Modifie le vecteur de sorte que un 1 signifie que, pour la colonne correspondante, la valeur sur la diagonale est 0.
//Autrement dit, cette coordonée doit être choisie égale à 0 ou 1 dans la solution du système.
//
//Renvoie count tel que 2^count soit le nombre de solutions possibles.

{
    uint count = 0;                                          //Donne le nombre de solution (2^count) de la matrice

    //------------------------------------------------------------------------------------------------------

    for(uint j=0;j<ligne;j++)                                 //On parcourt les differentes colonnes
    {
        uint i;
        for(i=0;i<ligne;i++)                                 //On parcourt les differentes lignes
        {                
            if(COORD(i,j)==1)                               //Si on trouve un 1, la ligne i sera le pivot:
            {
                if(i>=j){break;}                            //          -Si le 1 est sous la diagonale
                if(COORD(i,i)==0){break;}                   //          -Sinon, si la valeur de cette ligne à la diagonale est 0
            }
        }
        if(i==ligne){SETMODIF(j);count++;continue;}          //Si aucun pivot est trouvé, on met un 1 dans set à cette coordonée et on passe à la colonne suivante
                                                            
        if(i!=j)                                            //Si le pivot n'est pas sur la diagonale, on échange les lignes
        {
            long a;
            for(uint k=0;k<ligne_l;k++)
            {
                a = Registres_b[i][k];
                Registres_b[i][k] = Registres_b[j][k];
                Registres_b[j][k] = a;
            }
        }

        for(i=0;i<ligne;i++)                                 //On xor le pivot avec chaqu'une des lignes
        {
            if(i!=j && COORD(i,j)==1)
            {
                for(uint k=0;k<ligne_l;k++)
                {
                    Registres_b[i][k] = Registres_b[i][k]^Registres_b[j][k];
                }
            }
        }
    }
    
    //------------------------------------------------------------------------------------------------------

    return count;                                           //On sort de la fonction en renvoyant count
}


void resoud(long** Registres_b, long* Set, long* Sol, uint ligne, uint ligne_l, long seed)
//Prend en entrée la matrice réduite par "reduit", le vecteur "Set" donné par "reduit",
//un vecteur Sol etant la solution du système, un vecteur Y_dec qui contient la decomposition de Y,
//count représentant le nombre de solutions (2^count) et seed un graine de la solution.
//
//Cette fonction vas fixer les differentes valeurs à choisir, définies par "Set", grâce à "seed",
//Puis vas calculer les valeurs restantes avec la matrice réduite.
//
//A la fin, Set sera tel que le i-ème bit de Set est égale a 1 si Qi est à prendre, 0 sinon.
//Note !!! Le Set est actuellement 0 partout et 1 en les valeurs qu'on doit choisir, on remaque que si on met
//         toutes les valeurs à choisir à 0, on trouve la solution où tout est nul.

{   
    uint seed_curs = 0;                                      //seed_curs indique quel bit de la graine regardé

    for(uint i=0;i<(ligne_l);i++){Sol[i]=Set[i];}            //Sol copie Set


    //------------------------------------------------------------------------------------------------------

    for(uint j=0;j<ligne;j++)                                 //On parcourt Set en largeur
    {
        if(SETBIT(j)==1)                                        //Si le bit est 1, on doit choisir cette valeur
        {                                                    //On prend le bit correspondant à la position de "curs" sur la graine
            if(((seed>>seed_curs)&1)==0){SOLMODIF(j);}      //Si c'est un 0, Qj n'est pas pris dans la solution et rien ne change
            else
            {
                for(uint i=0;i<j;i++)                        //Si c'est un 1, Qj est pris et les lignes de la matrice sont "équilibré",
                //                                          //ie. pour toutes les lignes i ou la colonnes j est 1,                
                {//                                         //on inverse le i-ème bit de Sol.
                    if(COORD(i,j)==1){SOLMODIF(i);}
                }
            }
            if(seed_curs<31){seed_curs++;}                                    //On a utilise le bit de la graine, on passe au bit suivant.
        }
    }
}


uint test_sol(uint* Base_fact, uint** Registres_fact, mpz_t* Liste_A, long* Sol, uint* exp, mpz_t p, mpz_t X, mpz_t Y, mpz_t N, uint ligne, uint ligne_l)
//Prend en entrée principalement Sol, Registre_fact, Liste A et N
//Calcul, à partir de la solution fournie, X et Y tel que X² = Y² [N] et test si pgcd(X-Y,N) ou pgcd(Y+X,N) fournit un facteur non trivial de N.
//X est la racine du produit des Qi tel que i est choisi dans la solution
//Y est le produit des Ai tel que i est choisi dans la solution 
//Revoit 1 si le pgcd est non-trivial, 0 sinon
{
    for(uint j=0;j<ligne;j++) // Parcours de Sol pour calculer X et Y
    {
        if (SOLBIT(j)==1) // Si =1 on garde le Q
        {
            mpz_mul(Y,Y,Liste_A[j]); // On multiplie Y par le nouveau A
            mpz_mod(Y,Y,N); // Réduction modulo N

            for(uint i=1;i<ligne;i++) // Pour chaque premier pi, on ajoute la valuation pi-adique de Q(j) à exp[i]
                {exp[i]+=Registres_fact[i][j];}
        }
    }

    for(uint i=0;i<ligne;i++)
    {
        exp[i]/=2; // Il faut prendre la racine carrée du produit des Q
        mpz_set_ui(p,Base_fact[i]); // p=pi
        mpz_powm_ui(p,p,exp[i],N); // p=pi^(exp[i]) modulo N
        mpz_mul(X,X,p); // On multiplie X par la bonne puissance de pi
        mpz_mod(X,X,N); // Réduction modulo N
    }
    gmp_printf("X: %Zd Y: %Zd\n",X,Y);
    mpz_t A;
    mpz_init(A);
    mpz_add(A,X,Y);
    mpz_mod(A,A,N);
    gmp_printf("somme: %Zd  " ,A);
    mpz_gcd(A,A,N);
    gmp_printf("gcd: %Zd\n" ,A);

    if((mpz_cmp_ui(A,1)!=0)&&(mpz_cmp(A,N)!=0)) // Test gcd non-trivial
    {
        gmp_printf("Un diviseur non-trivial de %Zd est %Zd\n",N,A);
        return 1;
    }

    mpz_sub(A,X,Y);
    mpz_mod(A,A,N);
    gmp_printf("difference: %Zd  " ,A);
    mpz_gcd(A,A,N);
    gmp_printf("gcd: %Zd\n" ,A);

    if((mpz_cmp_ui(A,1)!=0)&&(mpz_cmp(A,N)!=0)) // Test gcd non-trivial
    {
        gmp_printf("Un diviseur non-trivial de %Zd est %Zd\n",N,A);
        return 1;
    }
    mpz_clear(A);

    return 0;
}