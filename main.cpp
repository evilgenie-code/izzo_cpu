#include <cmath>
#include <cstdio>
#include <cstdlib>


void vers(const double *V_in, double *Ver_out)
{
    double v_mod = 0;
    int i;

    for (i = 0; i < 3; i++)
    {
        v_mod += V_in[i] * V_in[i];
    }

    double sqrtv_mod = sqrt(v_mod);

    for (i = 0; i < 3; i++)
    {
        Ver_out[i] = V_in[i] / sqrtv_mod;
    }
}

void vett(const double *vet1, const double *vet2, double *prod)
{
    prod[0] = (vet1[1] * vet2[2] - vet1[2] * vet2[1]);
    prod[1] = (vet1[2] * vet2[0] - vet1[0] * vet2[2]);
    prod[2] = (vet1[0] * vet2[1] - vet1[1] * vet2[0]);
}

double x2tof(const double &x, const double &s, const double &c, const int lw, const int m)
{
    double am, a, alfa, beta;

    am = s / 2;
    a = am / (1 - x * x);

    if (x < 1)//ellpise
    {
        beta = 2 * asin(sqrt((s - c) / (2 * a)));
        if (lw) beta = -beta;
        alfa = 2 * acos(x);
    }
    else
    {
        alfa = 2 * acosh(x);
        beta = 2 * asinh(sqrt((s - c) / (-2 * a)));
        if (lw) beta = -beta;
    }

    if (a > 0)
    {
        return (a * sqrt(a)* ((alfa - sin(alfa)) - (beta - sin(beta)) + 2 * acos(-1.0) * m));
    }
    else
    {
        return (-a * sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta)));
    }

}

void lambert(const float *r0, const float *rk, float t, int lw, int revs, float mu)
{
    double R1 = 0.0, R2 = 0.0, dotProd = 0.0, mcrsProd, theta, c, s, T, q,
            T0 = 0.0, V, am, lambda, x1, x2, y1, y2, xNew, yNew, x;
    double  r1Unit[3], r2Unit[3], crsProd[3], ucrsProd[3], th1Unit[3], th2Unit[3],
            r1Double[3], r2Double[3], Tprod0[4] = { T0, 0.0, 0.0, 0.0 }, Td, phr, v1[3], v2[3];
    int leftBranch = 0;


    if (t <= 0)
    {
        printf("time is not valid");
        return;
    }

    for (int j = 0; j < 3; j++)
    {
        r1Double[j] = r0[j];
        r2Double[j] = rk[j];
    }

    for (int j = 0; j < 3; j++)
    {
        R1 += pow(r1Double[j], 2);
        R2 += pow(r2Double[j], 2);
    }

    R1 = sqrt(R1);
    R2 = sqrt(R2);

    V = sqrt(mu / R1);
    T = R1 / V;
    t /= T;

    for (int j = 0; j < 3; j++)
    {
        dotProd += (r1Double[j] * r2Double[j]);
    }

    theta = acos(dotProd / R2);

    if (lw)
    {
        theta = 2 *  acos(-1.0) - theta;
    }

    c = sqrt(1 + R2 * (R2 - 2.0 * cos(theta)));
    s = (1 + R2 + c) / 2.0;
    am = s / 2.0;
    lambda = sqrt(R2) * cos(theta / 2.0) / s;

    double inn1, inn2;

    const double tolerance = 1e-11;

    int err = 1;
    int iterate = 0;
    int imax = 30;

    if (revs == 0)
    {
        x1 = log(0.4767);
        x2 = log(1.5233);
        y1 = log(x2tof(-.5233, s, c, lw, revs)) - std::log(t);
        y2 = log(x2tof(.5233, s,c,lw, revs)) - std::log(t);

        while (err > tolerance && y1 != y2)
        {
            iterate++;
            xNew = (x1 * y2 - y1 * x2)/ (y2 - y2);
            yNew = log(x2tof(exp(xNew) - 1, s, c, lw, revs)) - std::log(t);
            x1 = x2;
            y1 = y2;
            x2 = xNew;
            y2 = yNew;
            x = exp(xNew) - 1;
        }
    }
    else
    {
        inn1 = 0.7234;
        inn2 = 0.5234;

        if (leftBranch)
        {
            inn1 = -0.5234;
            inn2 = -0.2234;
        }

        x1 = tan(inn1 * acos(-1.0) / 2);
        x2 = tan(inn2*acos(-1.0) / 2);
        y1 = x2tof(inn1, s, c, lw, revs) - t;
        y2 = x2tof(inn2, s, c, lw, revs) - t;

        while ((err > tolerance) && (y1 != y2) && iterate < imax)
        {
            iterate++;
            xNew = (x1*y2 - y1 * x2) / (y2 - y1);
            yNew = x2tof(atan(xNew) * 2 / acos(-1.0), s, c, lw, revs) - t;
            x1 = x2;
            y1 = y2;
            x2 = xNew;
            y2 = yNew;
            err = std::abs(x1 - xNew);
        }

        x = atan(xNew) * 2 / acos(-1.0);
    }

    double a = am / (1 - x * x); // solution semimajor axis
    double beta, alfa, psi, eta2, eta;

    if (x < 1)// ellipse
    {
        beta = 2 * asin(sqrt((s - c) / (2 * a)));
        if (lw) beta = -beta;
        alfa = 2 * acos(x);
        psi = (alfa - beta) / 2;
        eta2 = 2 * a * pow(sin(psi), 2) / s;
        eta = sqrt(eta2);
    }
    else
    {
        beta = 2 * asinh(sqrt((c - s) / (2 * a)));
        if (lw) beta = -beta;
        alfa = 2 * acosh(x);
        psi = (alfa - beta) / 2;
        eta2 = -2 * a * pow(sinh(psi), 2) / s;
        eta = sqrt(eta2);
    }

    double p = R2 / (am * eta2) * pow(sin(theta / 2), 2);
    double sigma1 = (1 / (eta * sqrt(am))) * (2 * lambda * am - lambda + x * eta);
    double ihDum[3], ih[3], dum[3];
    vett(r1Double, r2Double, ihDum);
    vers(ihDum, ih);

    if (lw)
    {
        for (int i = 0; i < 3; ++i)
        {
            ih[i] = -ih[i];
        }
    }

    double vr1 = sigma1;
    double vt1 = sqrt(p);
    vett(ih, r1Double, dum);


}

int main() {
    double AU = 1.49597870691e8;
    double fMSun = 1.32712440018e11;             // km^3/sec^2
    double fME = 398600.4415;                    // km^3/sec^2
    double RE = 6371;                            // km
    double GaccE = 9.80665 * 1.e-3;                // km/sec^2

    double UnitR = AU;
    double UnitV = sqrt(fMSun / UnitR);          // km/sec
    double UnitT = (UnitR / UnitV) / 86400;         // day
    double UnitA = fMSun / (UnitR * UnitR);

    float unitT = 100.0 / UnitT;
    float mu = 1.0;
    int lw = 1.0, revs = 0.0;
    float r1[3] = {-7.8941608095246896e-01, -6.2501194900473045e-01, 3.5441335698377735e-05};
    float r2[3] = {1.3897892184188783e+00, 1.3377137029002054e-01, -3.1287386211010106e-02};

    lambert(r1, r2, unitT, lw, revs, mu);
//    return 0;
}
