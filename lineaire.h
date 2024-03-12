uint reduit(long** Registres_b, long* Set, uint ligne, uint ligne_l);
void resoud(long** Registres_b, long* Set, long* Sol, uint ligne, uint ligne_l,long seed);
uint test_sol(uint* Base_fact, uint** Registres_fact, mpz_t* Liste_A, long* Sol, uint* exp, mpz_t p, mpz_t X, mpz_t Y, mpz_t N, uint ligne, uint ligne_l);