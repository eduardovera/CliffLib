#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <metric.h>
#include <functional>
#include <stdexcept>
#include <vector>

namespace CliffLib {

    const int N_DIMS = 4;

    typedef std::bitset<1 << N_DIMS> mask;

    template<std::size_t N>
    struct Comparer {
        bool operator () (const std::bitset<N> &x, const std::bitset<N> &y) const {
            return x.to_ulong() < y.to_ulong();
        }
    };


    template<class coeff_type>
    class multivector {

        private:
            std::map<mask, coeff_type, Comparer<1 << N_DIMS>> M;

            std::string get_base(std::bitset<1 << N_DIMS> m) const {
                std::string ret = "";
                mask binary = m;
                for (int i = 1; i <= N_DIMS; i++) {
                    mask b = binary >> (i);
                    if (b.test(0)) {
                        ret.append("(e" + std::to_string(i) + ")");
                        b.flip(0);
                        if (b.any()) {
                            ret.append("^");
                        }
                    }
                }
                return ret;
            }

        public:

            template<class T> friend multivector<T> e(int);
            template<class T> friend multivector<T> scalar(T);
            template<class T> friend multivector<T> take_grade(const multivector<T> &, const int);
            template<class T> friend int get_grade(const multivector<T> &);
            template<class T> friend multivector<T> REVERSE(const multivector<T> &);
            template<class T> friend multivector<T> INV(const multivector<T> &, Metric<T> &);
            template<class T> friend T SQR_NORM_REVERSE(const multivector<T> &, Metric<T> &);

            template<class T> friend std::ostream& operator << (std::ostream &, const multivector<T> &);

            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator + (const multivector<T> &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator - (const multivector<T> &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator * (const T &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator ^ (const multivector<T> &, const multivector<U> &);

            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> generic_call (const multivector<T> &m1, const multivector<U> &m2, auto f );

            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> RP (const multivector<T> &, const multivector<U> &, int);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> GP (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> IGP (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> LCONT (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> RCONT (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> DELTA (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend typename std::common_type<T, U>::type  SCP (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);

            template<class T> friend std::vector<multivector<T>> FACTORIZE(const multivector<T> &, Metric<T> &, T &);

            template<class T> friend multivector<T> DUAL(const multivector<T> &, int DIMS, Metric<T> &);
            template<class T, class U> friend multivector<T> DUAL(const multivector<T> &, const multivector<U> &, Metric<Metric<typename std::common_type<T, U>::type>> &);

            template<class T> friend multivector<T> UNDUAL(const multivector<T> &, int DIMS, Metric<T> &);
            template<class T, class U> friend multivector<T> UNDUAL(const multivector<T> &, const multivector<U> &, Metric<Metric<typename std::common_type<T, U>::type>> &);

            template<class T, class U> friend std::vector<multivector<typename std::common_type<T, U>::type>> MEET_AND_JOIN(const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);

    };

    template<class T>
    int get_grade(const multivector<T> &m) {
        auto it = m.M.begin();
        int grade = it->first.count();
        while (it != m.M.end()) {
            if ((int)it->first.count() != grade) {
                return -1;
            }
            ++it;
        }
        return grade;
    }

    template<class T>
    std::ostream& operator << (std::ostream &os, const multivector<T> &m) {
        for (auto it = m.M.begin(); it != m.M.end(); ++it) {
            os << (it->second > 0.0 ? "+" : "") << it->second << m.get_base(it->first);
        }
        return os;
    }

    template<class T = double>
    multivector<T> scalar(T value) {
        multivector<T> r;
        r.M[mask(0)] = value;
        return r;
    }


    template<class T = double>
    multivector<T> e(int index) {
        multivector<T> r;
        r.M[mask(1 << index)] = 1;
        return r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator - (const multivector<T> &m1, const multivector<U> &m2) {
        return m1 + (-1 * m2);
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator + (const multivector<T> &m1, const multivector<U> &m2) {
        multivector<typename std::common_type<T, U>::type> multivector_r;
        auto it1 = m1.M.begin();
        auto it2 = m2.M.begin();
        while (it1 != m1.M.end() && it2 != m2.M.end()) {
            if (it1->first.to_ullong() < it2->first.to_ullong()) {
                multivector_r.M[it1->first] = it1->second;
                it1++;
            } else if (it1->first.to_ullong() > it2->first.to_ullong()) {
                multivector_r.M[it2->first] = it2->second;
                it2++;
            } else {
                multivector_r.M[it2->first] = it1->second + it2->second;
                it1++;
                it2++;
            }
        }

        while (it1 != m1.M.end()) {
            multivector_r.M[it1->first] = it1->second;
            it1++;
        }

        while (it2 != m2.M.end()) {
            multivector_r.M[it2->first] = it2->second;
            it2++;
        }
        return multivector_r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator + (const T m1, const multivector<U> &m2) {
        return scalar(m1) + m2;
    }


    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator + (const multivector<T> &m1, const U m2) {
        return m1 + scalar(m2);
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator - (const T m1, const multivector<U> &m2) {
        return scalar(m1) + (-1 * m2);
    }


    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator - (const multivector<T> &m1, const U m2) {
        return m1 + scalar(-1 * m2);
    }


    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator * (const T &s, const multivector<U> &m) {
        multivector<typename std::common_type<T, U>::type> r;
        for (auto it = m.M.begin(); it != m.M.end(); ++it) {
            r.M[(*it).first] = s * (*it).second;
        }
        return r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator * (const multivector<U> &m, const T &s) {
        return s * m;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> generic_call (const multivector<T> &m1, const multivector<U> &m2, auto f ) {
        multivector<typename std::common_type<T, U>::type> result;
        for (auto it1 = m1.M.begin(); it1 != m1.M.end(); ++it1) {
            for (auto it2 = m2.M.begin(); it2 != m2.M.end(); ++it2) {
                result = result + f(it1->second, it1->first, it2->second, it2->first);
            }
        }
        return result;
    }

    int canonical_sort(mask mask_1, mask mask_2) {
        int swaps = 0;
        mask mask_it = mask_1 >> 1;
        while (mask_it.count() != 0) {
            swaps += (mask_it & mask_2).count();
            mask_it = mask_it >> 1;
        }
        if (!(mask(swaps) & mask(1)).any()) {
            return 1;
        }
        return - 1;
    }


    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> operator ^ (const multivector<T> &m1, const multivector<U> &m2) {
        auto OP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            if (!(mask1 & mask2).any()) {
                multivector_r.M[(mask1 | mask2)] = canonical_sort(mask1, mask2) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, OP);
    }

    template<class T>
    multivector<T> take_grade(const multivector<T> &m, const int grade) {
        multivector<T> multivector_r;
        for (auto it = m.M.begin(); it != m.M.end(); ++it) {
            if ((int)(it->first).count() == grade) {
                multivector_r.M[(*it).first] = (*it).second;
            }
        }
        return multivector_r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> RP(const multivector<T> &m1, const multivector<U> &m2, int NDIMS) {
        auto RP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_r = mask1 & mask2;
            if ((int)(mask1.count() + mask2.count() - mask_r.count()) == NDIMS) {
                multivector_r.M[mask_r] = canonical_sort(mask1 ^ mask_r, mask2 ^ mask_r) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, RP);
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> GP (const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        auto GP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_r = mask1 ^ mask2;
            multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1 & mask2).to_ulong(), (mask1 & mask2).to_ulong()) * coef1 * coef2;
            return multivector_r;
        };
        return generic_call(m1, m2, GP);
    }

    template<class T, class U>
    typename std::common_type<T, U>::type  SCP (const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        auto SCP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_r = mask1 ^ mask2;
            if (mask_r.count() == 0) {
                multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1 & mask2).to_ullong(), (mask1 & mask2).to_ulong()) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, SCP).M[0];
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> LCONT (const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        auto LCONT = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_r = mask1 ^ mask2;
            if (mask_r.count() == mask2.count() - mask1.count()) {
                multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1 & mask2).to_ullong(), (mask1 & mask2).to_ulong()) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, LCONT);
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> RCONT (const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        auto RCONT = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_r = mask1 ^ mask2;
            if (mask_r.count() == mask1.count() - mask2.count()) {
                multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1 & mask2).to_ullong(), (mask1 & mask2).to_ulong()) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, RCONT);
    }

    template<class T>
    multivector<T> REVERSE(const multivector<T> &m) {
        // assert that m is a blade
        multivector<T> multivector_r;
        auto it = m.M.begin();
        int g = get_grade(m);
        int ex = ((int)(g * (g-1))) >> 1;
        multivector_r.M[it->first] = std::pow(-1, ex) * m.M.at(it->first);
        return multivector_r;
    }

    template<class T>
    T SQR_NORM_REVERSE(const multivector<T> &m, Metric<T> &metric) {
        // assert that m is a blade
        return SCP(m, REVERSE(m), metric);
    }

    template<class T>
    multivector<T> INV(const multivector<T> &m, Metric<T> &metric) {
        // assert that m is a blade
        multivector<T> multivector_r;
        multivector<T> temp = REVERSE(m);
        T norm = SQR_NORM_REVERSE(m, metric);
        auto it = temp.M.begin();
        multivector_r.M[it->first] = it->second / norm;
        return multivector_r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> IGP(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        multivector<typename std::common_type<T, U>::type> multivector_r = GP(m1, INV(m2, metric), metric);
        return multivector_r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> DELTA(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        multivector<typename std::common_type<T, U>::type> temp = GP(m1, m2, metric);
        unsigned long long max_grade = 0;
        for (auto it = temp.M.begin(); it != temp.M.end(); ++it) {
            if (it->first.count() > max_grade) {
                max_grade = it->first.count();
            }
        }
        return take_grade(temp, max_grade);
    }

    template<class T, class U>
    multivector<T> PROJECTION(const multivector<T> &m1, const multivector<U> &m2, Metric<T> &metric) {
        return LCONT(LCONT(m1, INV(m2, metric), metric), m2, metric);
    }

    template<class T>
    std::vector<multivector<T>> FACTORIZE(const multivector<T> &m, Metric<T> &metric, T &scaling) {
        // assert it is a blade
        scaling = sqrt(SQR_NORM_REVERSE(m, metric));
        double inv_scaling = 1.0 / scaling;

        auto maxIt = m.M.begin();

        for (auto it = m.M.begin(); it != m.M.end(); ++it) {
            if (maxIt->second < it->second) {
                maxIt = it;
            }
        }

        mask masc = maxIt->first;

        multivector<T> temp = m * inv_scaling;
        std::vector<multivector<T>> r;

        for (int i = 1; i < 1 << N_DIMS; i++) {
            if (masc.test(i)) {
                mask e_ = masc.none();
                e_[i] = 1;
                multivector<T> e_j;
                e_j.M[1 << i] = 1;

                multivector<T> proj = PROJECTION(e_j, temp, metric);
                multivector<T> factor_j = proj * inv_scaling;
                r.push_back(factor_j);
                temp = LCONT(INV(factor_j, metric), temp, metric);
            }
        }

        return r;
    }

    template<class T>
    multivector<T> DUAL(const multivector<T> &m, int DIMS, Metric<T> &metric) {
        mask masc;
        for (int i = 1; i <= DIMS; i++) {
            masc[i] = 1;
        }
        multivector<T> I;
        I.M[masc] = 1;

        return LCONT(m, INV(I, metric), metric);
    }

    template<class T, class U>
    multivector<T> DUAL(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        return LCONT(m1, INV(m2, metric), metric);
    }

    template<class T>
    multivector<T> UNDUAL(const multivector<T> &m, int DIMS, Metric<T> &metric) {
        mask masc;
        for (int i = 1; i <= DIMS; i++) {
            masc[i] = 1;
        }
        multivector<T> I;
        I.M[masc] = 1;

        return LCONT(m, I, metric);
    }

    template<class T, class U>
    multivector<T> UNDUAL(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        return LCONT(m1, m2, metric);
    }


    template<class T, class U>
    std::vector<multivector<typename std::common_type<T, U>::type>> MEET_AND_JOIN(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {

        multivector<T> A = m1;
        multivector<T> B = m2;

        int r = get_grade(m1);
        int s = get_grade(m2);

        if (r > s) {
            multivector<T> A = m2;
            multivector<T> B = m1;
        }

        multivector<typename std::common_type<T, U>::type> delta = DELTA(A, B, metric);
        int t = (r + s - get_grade(delta)) >> 1;


    }
}





#endif // MULTIVECTOR_H
