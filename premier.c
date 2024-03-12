#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gmp.h>
#include "premier.h"
#include "macro.h"

int legendre(int p, mpz_t N)
/*Revoit le symbole de legendre de N sur p*/
{
	uint t = 1;  
    //p mod 8
	uint r; 
    //N mod p
	uint n = mpz_fdiv_ui(N,p);

	while(n != 0)
	{
		while(n%2 == 0)
		{
			n/=2;
			r = p%8;
			if((r == 3)||(r == 5)){t = -t;}
		}
		n=n-p;p=n+p;n=p-n;
		if(((n%4)==3)&&((p%4)==3)){t = -t;}
		n = n%p;
	}
	if(p==1){return t;}
	else{return 0;}
}

uint premier(uint longueur)
//Ecrit les nombres premiers strictements inferieurs à N dans Liste_Premier.txt.
//Revoit le nombre de premiers trouvés.
{

    //------------------------------------------------------------------------------------------------------
    //Initialisation

    //compteur
    uint count = 0;
    
    //Création de Liste_Premier.txt
    FILE * l;                              
    l = fopen("Liste_Premier.txt","w");

    //Création d'un tableau de N case, 1 si la case est premier, 0 sinon
    char* T;                                 
    T = (char *) malloc(longueur * sizeof(char)) ;
    T[0] = 0; T[1] = 0;                     
    for(uint i = 2; i<longueur; i++){T[i] = 1;}

    //------------------------------------------------------------------------------------------------------
    //Traitement

    //Parcourt du tableau, si la case est un premier, on met 0 dans tout ses multiples
    for(uint i=2; i<(longueur); i++)                  
    {
        if(T[i]==1)                         
        {
            //Ecriture dans le .txt
            fprintf(l,"%d ",i);             
            count++;
            unsigned long d = 2*i;
            while(d<longueur)
            {
                T[d] = 0;
                d += i;
            }
        }
    }

    //------------------------------------------------------------------------------------------------------
    //Fermeture

    fclose(l);
    free(T);
    return count;

    //------------------------------------------------------------------------------------------------------
}

uint premier_rec(uint* List, uint longueur, mpz_t N)
/*Recupère les nombres premiers enregistrés dans Liste_Premier.txt, tel que le symbole de legendre soit 1,
et les range dans un tableau List
arguments :   Liste, tableau a remplir, de longueur "longueur"
              longueur, longueur de Liste
              N nombre a factorisé*/
{
    //------------------------------------------------------------------------------------------------------
    //Lit Liste_Premier.txt et le met dans le string text

    FILE * file;
    char* text;
    long numbytes;

    file = fopen("Liste_Premier.txt", "r");
    if(file == NULL){printf("error: Liste_Premier.txt non existant, Use '-pl uint size' to create\n");return 0;}
    
    fseek(file, 0L, SEEK_END);
    numbytes = ftell(file);
    fseek(file, 0L, SEEK_SET);  
 
    text = (char*)calloc(numbytes, sizeof(char));   
    if(text == NULL){printf("error: allocation memoire premier_rec\n");return 0;}
 
    fread(text, sizeof(char), numbytes, file);
    fclose(file);

    //------------------------------------------------------------------------------------------------------
    //Decoupe les première valeur de text et les met dans List

    char * token = strtok(text, " ");

    List[0]=1;
    uint temp;
    uint i=1;
    while(i<longueur) 
    {
        if(token == NULL){printf("Liste de Premiers trop courte\n");return 0;}
        temp = atoi(token); 
        if(legendre(temp,N) == 1){List[i] = temp;i++;}
        token = strtok(NULL, " ");
    }

    //------------------------------------------------------------------------------------------------------
    //Fermeture
    free(text);
    return 1;
}