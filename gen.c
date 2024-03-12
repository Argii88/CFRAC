#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include "gen.h"
#include "macro.h"

	void gen(uint* Base_fact, uint ligne, uint colonne, mpz_t N, mpz_t* Liste_Q, mpz_t* Liste_A, long** Registres_b, uint** Registres_fact)
		{
		//Arguments :
		//	- Base_fact la base de factorisation de petits nombres premiers
		//	- ligne la taille de Base_fact
		//  - colonne le nombre de Qn/An à calculer (généralement ligne = colonne)
		//	- N le grand entier à factoriser
		//	- Liste_Q la liste des valeurs de Qi à conserver
		//	- Liste_A la liste des valeurs de Ai correspondant aux Qi
		//	- Registres_b une matrice d'entier, chacun ayant 64 bits, représentant les exposants mod 2 des premiers de Base_fact dans l'entier Q
		//	- Registres_fact une matrice d'entiers non-signés représentant les exposants des premiers de Base_fact dans l'entier Q

		// On suit l'article de Pomerance/Wagstaff pour le calcul des Qi et Ai

		//------------------------------------------------------------------------------------
		//Initialisation

		mpz_t g, Q1, Q2, q, r1, r2, A1, A2, G, dr;
		
		mpz_inits(g, Q1, Q2, q, r1, r2, A1, A2, G, dr, NULL);
		
		mpz_sqrt(g,N); // g=floor(sqrt(N))

		mpz_set(Q1,N); // Q(-1)=N
		mpz_set_ui(Q2,1); // Q(0)=1
		mpz_set(q,g); // q(0)=g
		mpz_set(r1,g); // r(-1)=g
		//r2=r(0)=0 déjà initialisé
		mpz_set_ui(A1,1); // A(-1)=1
		mpz_set(A2,g); // A(0)=g
		
		//------------------------------------------------------------------------------------
		// Calcul itératif des Qi

		mpz_mul_2exp(g,g,1); // On n'a plus besoin de g, on aura besoin de const 2g

		uint j=0; // Nombre de Qi déjà enregistrés
		unsigned long neg = 0; //si Q est negatif
		uint nb = 0;

		while (j<colonne) // généralement colonne = ligne pour avoir une relation linéaire non-triviale, et avoir une matrice carrée
		{
			//---------------------------------------------------------------------------
			// On commence par calculer les nouvelles valeurs des suites
			// neg alterne entre 0 et 1 car le signe des Qn alterne aussi entre + et -
			neg^=1;
			nb++;

			mpz_sub(dr,r2,r1);
			mpz_addmul(Q1,q,dr);
			mpz_swap(Q1,Q2); // Q2 = Q1+q.(r2-r1)

			mpz_sub(G,g,r2); // G = 2g-r2
    		
			mpz_fdiv_q(q,G,Q2); // q = G//Q2negation

			mpz_mul(r1,q,Q2);
			mpz_sub(r1,G,r1);
			mpz_swap(r1,r2); // r2 = G-q.Q2

			mpz_addmul(A1,q,A2);
			mpz_mod(A1,A1,N);
			mpz_swap(A1,A2); // A2 = q.A2+A1 mod N
			// Q2=Qi, A2=Ai, q=qi, r2=ri, G=Gi

			//---------------------------------------------------------------------------
			// Trial Division

			// indique si le Qn traité se factorise
			char incr;

			//On lance trail_div
			incr=trial_div(Registres_fact, Registres_b, j, Q2, Base_fact, ligne);

			if (incr==1) // Si Q2 a été complètement factorisé
			{
				Registres_b[0][j>>6]|=(neg<<(j&63));// On met la negation
				mpz_set(Liste_Q[j],Q2); // On ajoute Q2 à la liste des Q
				mpz_set(Liste_A[j],A1); // On ajoute A2 correspondant
				gmp_printf("Colonne %u -- %ld -- Q: %Zd, A: %Zd\n",j,colonne,Q2,A1);
				j++; // On a un nouveau Q factorisé
			}
		}



		//-------------------------------------------------------------------------------------
		// Clear
		printf("nb Q testé:%u\n",nb);
		mpz_clears(g, Q1, Q2, q, r1, r2, A1, A2, G, dr,NULL);
		}

	char trial_div(uint** Registres_fact, long** Registres_b, uint i, mpz_t Q, uint* Base_fact, uint ligne)
		{
		// La procédure calcule les multiplicités des éléments de Base_fact dans Q,
		// et les rentre dans les ièmes cases des registres
		// Si Q est complètement factorisé, la fonction renvoie 1, et 0 sinon.

		// Arguments:
		//	- Registres_fact et Registres_b les registres à remplir
		//	- Q le nombre à factoriser tel que Q=Qi
		//	- Base_fact la base de nombres premiers
		//	- ligne la taille de Base_fact

		uint ind=1; //indice de parcours de Base_fact
		uint mult; //compteur de la multiplicité du nombre premier dans Q

		mpz_t T; // mpz que l'on va diviser au lieu de le faire en place sur Q
		mpz_init(T);
		mpz_set(T,Q); // T prend la valeur Q

		while(ind<ligne)
			{
			// Calcul de la multiplicité de p_ind dans Q

			mult=0;
			while(mpz_divisible_ui_p(T,Base_fact[ind])!=0)
			{
				mpz_divexact_ui(T,T,Base_fact[ind]);  // on divise et on augmente la multiplicité de 1.
				mult+=1; 
			}

			Registres_fact[ind][i]=mult;		// On met mult dans Registre_fact
			Registres_b[ind][i>>6]=(Registres_b[ind][i>>6]&(~(((long)1)<<(i&63))))^(((long)(mult&1))<<(i&63));  // On met mult mod 2 dans registre_b

			ind++; // Recommencer avec le prochain premier de Base_fact
		}

		if (mpz_cmp_d(T,1)==0){mpz_clear(T);return 1;} // Si T=1 c'est que Q est complètement factorisé dans Base_fact
		else {mpz_clear(T);return 0;}
		}
