#include "integer_common.h"

//a and b are nonnegative
void xgcd_partial(integer& u, integer& v, integer& a, integer& b, const integer& L) {
    fmpz_t f_u; fmpz_init(f_u);
    fmpz_t f_v; fmpz_init(f_v);
    fmpz_t f_a; fmpz_init(f_a);
    fmpz_t f_b; fmpz_init(f_b);
    fmpz_t f_L; fmpz_init(f_L);

    fmpz_set_mpz(f_a, a.impl);
    fmpz_set_mpz(f_b, b.impl);
    fmpz_set_mpz(f_L, L.impl);

    fmpz_xgcd_partial(f_u, f_v, f_a, f_b, f_L);

    fmpz_get_mpz(u.impl, f_u);
    fmpz_get_mpz(v.impl, f_v);
    fmpz_get_mpz(a.impl, f_a);
    fmpz_get_mpz(b.impl, f_b);

    fmpz_clear(f_u);
    fmpz_clear(f_v);
    fmpz_clear(f_a);
    fmpz_clear(f_b);
    fmpz_clear(f_L);
}

void inject_error(mpz_struct* i) {
    if (!enable_random_error_injection) {
        return;
    }

    mark_vdf_test();

    double v=rand_integer(32).to_vector()[0]/double(1ull<<32);

    if (v<random_error_injection_rate) {
        print( "injected random error" );

        int pos=int(rand_integer(31).to_vector()[0]);
        pos%=mpz_sizeinbase(i, 2);
        mpz_combit(i, pos);
    }
}