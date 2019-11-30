#ifndef LMOMENTS_GEV_H
#define LMOMENTS_GEV_H

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "utils.h"

namespace lmoments {

template<typename T>
class GEV : public distribution<T> {
  protected:
    T u;
    T a;
    T g;
    bool converged_m = true;

  public:
    explicit inline GEV(T u_p = 0, T a_p = 0, T g_p = 0) : u(u_p), a(a_p), g(g_p) {}

    template<typename DataVector = std::vector<T>>
    static inline GEV from_data(const DataVector& x) {
        GEV res;
        res.fit_data(x);
        return res;
    }

    template<typename DataVector = std::vector<T>>
    inline void fit_data(const DataVector& x) {
        std::vector<T> xmom(3);
        lmoment_ratios_unbiased(x, xmom);
        set_lmoment_ratios(xmom);
    }

    inline bool converged() const {
        // if false: iteration has not converged, results may be unreliable
        return converged_m;
    }

    inline void set_lmoment_ratios(const std::vector<T>& xmom) override {
        constexpr int maxit = 20;
        constexpr T eu = .57721566;   // Euler's Constant
        constexpr T dl2 = .69314718;  // log(2)
        constexpr T dl3 = 1.0986123;  // log(3)
        constexpr T a0 = .2837753;
        constexpr T a1 = -1.21096399;
        constexpr T a2 = -2.50728214;
        constexpr T a3 = -1.13455566;
        constexpr T a4 = -.07138022;
        constexpr T b1 = 2.06189696;
        constexpr T b2 = 1.31912239;
        constexpr T b3 = .25077104;
        constexpr T c1 = 1.59921491;
        constexpr T c2 = -.48832213;
        constexpr T c3 = .01573152;
        constexpr T d1 = -.64363929;
        constexpr T d2 = .08985247;
        constexpr T small = 1e-5;
        constexpr T eps = 1e-6;

        const T t3 = xmom[2];
        if (xmom[1] <= 0 || std::abs(t3) >= 1) {
            throw std::invalid_argument("lmoments invalid");
        }
        converged_m = true;
        if (t3 > 0) {
            const T z = 1 - t3;
            g = (-1 + z * (c1 + z * (c2 + z * c3))) / (1 + z * (d1 + z * d2));
            if (std::abs(g) < small) {
                g = 0;  // == GUM
                a = xmom[1] / dl2;
                u = xmom[0] - eu * a;
                return;
            }
        } else {
            g = (a0 + t3 * (a1 + t3 * (a2 + t3 * (a3 + t3 * a4)))) / (1 + t3 * (b1 + t3 * (b2 + t3 * b3)));
            if (t3 < -0.8) {
                if (t3 <= -0.97) {
                    g = 1 - std::log1p(t3) / dl2;
                }
                const T t0 = (t3 + 3) / 2;
                converged_m = false;
                for (int it = 0; it < maxit; ++it) {
                    const T x2 = std::pow(2, -g);
                    const T x3 = std::pow(3, -g);
                    const T xx2 = 1 - x2;
                    const T xx3 = 1 - x3;
                    const T t = xx3 / xx2;
                    const T deriv = (xx2 * x3 * dl3 - xx3 * x2 * dl2) / (xx2 * xx2);
                    const T gold = g;
                    g -= (t - t0) / deriv;
                    if (std::abs(g - gold) <= eps * g) {
                        converged_m = true;
                        break;
                    }
                }
            }
        }
        const T gam = std::exp(dlgama(g + 1));
        a = xmom[1] * g / (gam * (1 - std::pow(2, -g)));
        if (a <= 0) {
            throw std::invalid_argument("lmoments invalid");
        }
        u = xmom[0] - a * (1 - gam) / g;
    }

    inline void get_lmoment_ratios(std::vector<T>& xmom) const override {
        const int nmom = xmom.size();
        constexpr T zmom[20] = {.577215664901532861,   .693147180559945309,   .169925001442312363,   .150374992788438185,   .0558683500577583138,
                                .0581100239999710876,  .0276242584297309125,  .0305563766579053126,  .0164650282258328802,  .0187846624298170912,
                                .0109328215063027148,  .012697312667632953,   .00778982818057231804, .00914836179621999726, .00583332389328363588,
                                .00690104287590348154, .00453267970180679549, .00538916811326595459, .0036240776777236879,  .00432387608605538096};
        constexpr T small = 1e-6;
        if (g <= -1) {
            throw std::invalid_argument("parameters invalid");
        }
        if (nmom > 20) {
            throw std::invalid_argument("too many lmoment ratios requested");
        }
        if (std::abs(g) <= small) {
            xmom[0] = u;
            if (nmom > 1) {
                xmom[1] = a * zmom[0];
                if (nmom > 2) {
                    for (int i__ = 2; i__ < nmom; ++i__) {
                        xmom[i__] = zmom[i__ - 1];
                    }
                }
            }
        } else {
            const T gam = std::exp(dlgama(1 + g));
            xmom[0] = u + a * (1 - gam) / g;
            if (nmom > 1) {
                const T xx2 = 1 - std::pow(2, -g);
                xmom[1] = a * xx2 * gam / g;
                if (nmom > 2) {
                    T z0 = 1;
                    for (int j = 2; j < nmom; ++j) {
                        const T dj = j + 1;
                        const T beta = (1 - std::pow(dj, -g)) / xx2;
                        z0 = z0 * (4 * dj - 6) / dj;
                        T z__ = z0 * 3 * (dj - 1) / (dj + 1);
                        T sum = z0 * beta - z__;
                        if (j > 2) {
                            for (int i__ = 2; i__ <= j - 1; ++i__) {
                                const T di = i__;
                                z__ = z__ * (di + di + 1) * (dj - di) / ((di + di - 1) * (dj + di));
                                sum -= z__ * xmom[i__];
                            }
                        }
                        xmom[j] = sum;
                    }
                }
            }
        }
    }

    inline T cdf(T x) const override {
        constexpr T small = 1e-15;
        T y = (x - u) / a;
        if (g == 0.0) {
            return std::exp(-std::exp(-y));
        }
        const T arg = 1 - g * y;
        if (arg > small) {
            y = -std::log(arg) / g;
            return std::exp(-std::exp(-y));
        }
        if (g < 0) {
            return 0;
        }
        return 1;
    }

    inline T quantile(T f) const override {
        if (f <= 0.0 || f >= 1.0) {
            if ((f == 0 && g < 0) || (f == 1 && g > 0)) {
                return u + a / g;
            }
            throw std::invalid_argument("argument out of range");
        }
        T y = -std::log(-std::log(f));
        if (g != 0) {
            y = (1 - std::exp(-g * y)) / g;
        }
        return u + a * y;
    }

    inline std::vector<T> get_parameters() const override { return std::vector<T>({u, a, g}); }
};
}  // namespace lmoments

#endif
