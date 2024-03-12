#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<gmp.h>
#include "gen.h"
#include "premier.h"
#include "lineaire.h"
#include "macro.h"

int main(int argc, char *argv[])
{
    //Argument -pl (premier liste) suivi d'un int N (<~4 000 000 000) lance la fonction premier(N)
    if(strcmp(argv[1],"-pl")==0)
    {
        printf("%d\n",premier((uint)atoi(argv[2])));
		return 1;
    }

    //Argument -pp (premier print)  suivi d'un entier N puis d'un entier longueur (<~2 000 000),
	//affiche les premiers nombre de List_Premier, tel que le symbole de legendre de N sur p soit 1
    if(strcmp(argv[1],"-pp")==0)
    {
        uint longueur = atoi(argv[3]);
        uint List[longueur];
		mpz_t N;
		mpz_init(N);
		mpz_set_ui(N,atoi(argv[2]));
        if(premier_rec(List,longueur,N)==1)
        {
            for(uint i=0;i<longueur;i++){printf("%d\n",List[i]);}
        }
		return 1;
	}


	//Argument -fact suivi d'un entier N a factorisé, d'un entier ligne, représentant le nombre de premiers dans la base de factorisation,
	//puis d'un entier facultatif k (par defaut 1) par lequel N est multiplié (cas de periode faible)
	if(strcmp(argv[1],"-fact")==0)
	{
		//On récupère la valeur de N à factoriser
    	srand(time(NULL));
    	clock_t start = clock();

		mpz_t N;
		mpz_init(N);
		mpz_set_str(N,argv[2],10);

		//On détermine la taille l de la base de factorisation dont on aura besoin

		uint ligne; // ligne est la taille d'une ligne en bits
		ligne=(uint)atoi(argv[3]);
		uint ligne_l=((ligne-1)/64) + 1; //ligne_l est la taille en long (un long est 64 bits)

		//Si k est fourni on multiplie N par k
		uint k;
		if(argc>4){k = (uint)atoi(argv[4]);}
		else{k=1;}
		mpz_mul_ui(N,N,k);

		//Création de Base_fact

		uint* Base_fact;
		Base_fact=malloc(ligne*sizeof(int*));
		if(premier_rec(Base_fact,ligne,N)==0){return 0;};

		//Création des registres Registres_fact et Registres_b qui contiendront respectivement les exposants de factorisation des Qi et leur réduction mod 2

		long** Registres_b;
		Registres_b=malloc(ligne*sizeof(long*)); // Il y aura ligne registres, un par élément de Base_fact
		for(uint i=0;i<ligne;i++)
			{Registres_b[i]=calloc(ligne_l,sizeof(long));} // Chaque registre contiendra ligne_l entier long, chacun contenant 64=ligne/ligne_l bits d'information

		uint** Registres_fact;
		Registres_fact=malloc(ligne*sizeof(uint*));
		for(uint i=0;i<ligne;i++)
			{Registres_fact[i]=calloc(ligne,sizeof(uint));}

		//Création des listes de A et Q, Liste_A et Liste_Q

		mpz_t* Liste_A;
		mpz_t* Liste_Q;
		Liste_A = malloc(sizeof(mpz_t)*ligne);
		Liste_Q = malloc(sizeof(mpz_t)*ligne);
		for(uint i=0;i<ligne;i++)
		{
			mpz_init(Liste_A[i]);
			mpz_init(Liste_Q[i]);
		}

		//Remplissage des listes et registres

		gen(Base_fact, ligne, ligne, N, Liste_Q, Liste_A, Registres_b, Registres_fact);
		mpz_divexact_ui(N,N,k); //On recupère le N original

		//Partie linéaire pour déterminer les Qi à multiplier pour obtenir un carré

		long* Set;
		Set=calloc(ligne_l,sizeof(long));

		long* Sol;
		Sol=malloc(sizeof(long)*ligne_l);	

		//2^count est le nombre de solutions du système linéaire.
		uint count;
		// On réduit le système Registres_b, et count est égal au nombre de solutions
		count=reduit(Registres_b,Set,ligne,ligne_l); 
		printf("reduit done\ncount= %d\n",count);

		mpz_t X; //X sera la racine du produit des Qi
		mpz_t Y; //Y sera le produit des Ai
		mpz_init(X);
		mpz_init(Y);

		uint* exp; // Exposants des premiers dans X
		exp=malloc(ligne*sizeof(int));

		mpz_t p;
		mpz_init(p);

		for(uint k=1;(k<100) && (k<((long)1<<count));k++)
		{
			//Seed est random pour choisir au hasard une solution du système
			long seed = rand();
	        seed = (seed<<21)^rand();
	        seed = (seed<<21)^rand();

			//Appel de resoud
	        resoud(Registres_b,Set,Sol,ligne,ligne_l,seed);

			//On réinitialise exp X et Y.
			for(uint i=0;i<ligne;i++){exp[i]=0;}
	        mpz_set_ui(X,1);
	    	mpz_set_ui(Y,1);

			//Appel de test_sol
			if(test_sol(Base_fact, Registres_fact, Liste_A, Sol, exp, p, X, Y, N, ligne, ligne_l) == 1)
			{
				printf("seed = %d\n",k);
				break;
			}

	    }
	    
		clock_t end = clock();

	    printf("%f secondes\n",((double)end - start) / CLOCKS_PER_SEC);
	    printf("End\n");
		// Clear et free

		mpz_clear(p);

		free(exp);

		mpz_clear(X);
		mpz_clear(Y);

		free(Sol);
		free(Set);

		for(uint i=0;i<ligne;i++)
			{free(Registres_fact[i]);
			free(Registres_b[i]);
			mpz_clear(Liste_A[i]);
			mpz_clear(Liste_Q[i]);}

		free(Registres_fact);
		free(Registres_b);

		free(Liste_A);
		free(Liste_Q);

		mpz_clear(N);

	    return 1;
	}

	//Argument -try suivi d'un entier N a factorisé, d'un entier ligne, représentant le nombre de premiers dans la base de factorisation,
	//puis d'un entier facultatif k (par defaut 1) par lequel N est multiplié (cas de periode faible)
	//mesure le temps pour les 20 premières colonnes et estime le temps de résolution total
	if(strcmp(argv[1],"-try")==0)
	{
		//On récupère la valeur de N à factoriser
    	clock_t start = clock();

		mpz_t N;
		mpz_init(N);
		mpz_set_str(N,argv[2],10);

		//On détermine la taille l de la base de factorisation dont on aura besoin

		uint ligne; // ligne est la taille d'une ligne en bits
		ligne=(uint)atoi(argv[3]);

		//Si k est fourni on multiplie N par k
		uint k;
		if(argc>4){k = (uint)atoi(argv[4]);}
		else{k=1;}
		mpz_mul_ui(N,N,k);

		//Création de Base_fact

		uint* Base_fact;
		Base_fact=malloc(ligne*sizeof(int*));
		if(premier_rec(Base_fact,ligne,N)==0){return 0;};

		//Création des registres Registres_fact et Registres_b qui contiendront respectivement les exposants de factorisation des Qi et leur réduction mod 2

		long** Registres_b;
		Registres_b=malloc(ligne*sizeof(long*)); // Il y aura ligne registres, un par élément de Base_fact
		for(uint i=0;i<ligne;i++)
			{Registres_b[i]=calloc(1,sizeof(long));} // Chaque registre contiendra 1 entier long, contenant 64 bits d'information

		uint** Registres_fact;
		Registres_fact=malloc(ligne*sizeof(uint*));
		for(uint i=0;i<ligne;i++)
			{Registres_fact[i]=calloc(ligne,sizeof(uint));}

		//Création des listes de A et Q, Liste_A et Liste_Q

		mpz_t* Liste_A;
		mpz_t* Liste_Q;
		Liste_A = malloc(sizeof(mpz_t)*ligne);
		Liste_Q = malloc(sizeof(mpz_t)*ligne);
		for(uint i=0;i<ligne;i++)
		{
			mpz_init(Liste_A[i]);
			mpz_init(Liste_Q[i]);
		}

		gen(Base_fact, ligne, 20, N, Liste_Q, Liste_A, Registres_b, Registres_fact);
		clock_t end = clock();

		printf("%f secondes\n",((double)end - start) / CLOCKS_PER_SEC);
		printf("temps estimé: %f\n",(ligne/20)*((double)end - start) / CLOCKS_PER_SEC);
		printf("End\n");

		// Clear et free
		for(uint i=0;i<ligne;i++)
		{
			free(Registres_fact[i]);
			free(Registres_b[i]);
			mpz_clear(Liste_A[i]);
			mpz_clear(Liste_Q[i]);
		}
		free(Registres_fact);
		free(Registres_b);
		free(Liste_A);
		free(Liste_Q);
		mpz_clear(N);

		return 1;
	}

	return 0;
}
