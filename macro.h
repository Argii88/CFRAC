//MACRO
#define COORD(i,j) ((Registres_b[i][j>>6]>>(j&63))&1)
#define SETMODIF(i) Set[i>>6]^=(((long)1)<<(i&63))
#define SOLMODIF(i) Sol[i>>6]^=(((long)1)<<(i&63))
#define SETBIT(i) ((Set[i>>6]>>(i&63))&1)
#define SOLBIT(i) ((Sol[i>>6]>>(i&63))&1)

//type
typedef unsigned int uint;
