#ifndef _CURVES_H_
#define _CURVES_H_
#define ANOMALOUS_11_SIZE 6
/*
    Contains: 
    - The representation of a point P=(x:y:z)
*/
struct Point{
    mpz_t x;
    mpz_t y;
    int z;
};
typedef struct Point Point;


/*
    Contains: 
    - The representation of a curve in the short weierstrass form: E: y^2 = x^3 + ax + b
    - A point P ∈ E.
*/
struct Curve{
    mpz_t a;
    mpz_t b;
    Point point;
};
typedef struct Curve Curve;


char *CURVES_ANOMALOUS_11 [ANOMALOUS_11_SIZE][2][4] = {
    { {"-1056","13552"}, {"33","1", "121","1"} },
    {{"-127776","-18037712"},{"1417675406433","3406706689", "183812129731156729", "198839249316863"}},
    {{"-31944","-2254714"},{"3489035","4489", "6326406713","300763"}},
    {{"-264", "1694"},{"11","1", "11","1"}},
    {{"-3194400","-2254714000"},{"10841625811568724573222992364076595009","1097994948369262177403982660144025", "35065604066149379334762170587524440281066037914774322473","36383168104020123152509793517255545827474838376125"}},
    {{"-26400", "1694000"},{"-2431","25", "228503","125"}},
};

#endif