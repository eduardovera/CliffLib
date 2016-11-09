#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <utils.h>
#include <metric.h>
#include <functional>

#define MAX_DIMENSIONS 10
#ifdef __MINGW32__
typedef uint8_t uint;
#endif

template<int K = MAX_DIMENSIONS>
int canonical_sort(uint, uint);


template<class coeff_type, int k = MAX_DIMENSIONS>
class multivector {

    private:
        std::map<uint, coeff_type> M;

        std::string get_base(coeff_type index) const {
            std::string ret = "";
            std::bitset<1 << k> binary(index);
            for (int i = 1; i <= k; i++) {
                std::bitset<1 << k> b = binary >> (i);
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

        template<class T, int K> friend multivector<T, K> e(int);
        template<class T, int K> friend multivector<T, K> scalar(T);
        template<class T, int K> friend multivector<T, K> take_grade(const multivector<T, K> &, const int);
        template<class T, int K> friend int check_grade(const multivector<T, K> &);
        template<class T, int K> friend multivector<T, K> REVERSE(const multivector<T, K> &);
        template<class T, int K> friend multivector<T, K> INV(const multivector<T, K> &, Metric<T> &);
        template<class T, int K> friend T SQR_NORM_REVERSE(const multivector<T, K> &, Metric<T> &);

        template<class T, int K> friend std::ostream& operator << (std::ostream &, const multivector<T, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator + (const multivector<T, K> &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator - (const multivector<T, K> &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator * (const T &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator ^ (const multivector<T, K> &, const multivector<U, K> &);

        template<class T, class U, int K> friend typename std::common_type<T, U>::type  SCP (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> RP (const multivector<T, K> &, const multivector<U, K> &, int N_DIMS);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> GP (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> IGP (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> LCONT (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> RCONT (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> generic_call (const multivector<T, K> &m1, const multivector<U, K> &m2, auto f );
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> DELTA (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);

};

template<class coeff_type, int K = MAX_DIMENSIONS>
std::ostream& operator << (std::ostream &os, const multivector<coeff_type, K> &m) {
    for (auto it = m.M.begin(); it != m.M.end(); ++it) {
        os << ((*it).second > 0.0 ? "+" : "") << (*it).second << m.get_base((*it).first);
    }
    return os;
}

template<class coeff_type = double, int K = MAX_DIMENSIONS>
multivector<coeff_type, K> scalar(coeff_type value) {
    multivector<coeff_type, K> r;
    r.M[0] = value;
    return r;
}


template<class coeff_type = double, int K = MAX_DIMENSIONS>
multivector<coeff_type, K> e(int index) {
    multivector<coeff_type, K> r;
    r.M[1 << index] = 1;
    return r;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator - (const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2) {
    return m1 + (-1 * m2);
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator + (const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2) {
    multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> multivector_r;
    auto it1 = m1.M.begin();
    auto it2 = m2.M.begin();

    while (it1 != m1.M.end() && it2 != m2.M.end()) {
        if ((*it1).first < (*it2).first) {
            multivector_r.M[(*it1).first] = (*it1).second;
            it1++;
        } else if ((*it1).first > (*it2).first) {
            multivector_r.M[(*it2).first] = (*it2).second;
            it2++;
        } else {
            multivector_r.M[(*it2).first] = (*it1).second + (*it2).second;
            it1++;
            it2++;
        }
    }

    while (it1 != m1.M.end()) {
        multivector_r.M[(*it1).first] = (*it1).second;
        it1++;
    }

    while (it2 != m2.M.end()) {
        multivector_r.M[(*it2).first] = (*it2).second;
        it2++;
    }
    return multivector_r;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator + (const coeff_type1 m1, const multivector<coeff_type2, K> &m2) {
    return scalar(m1) + m2;
}


template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator + (const multivector<coeff_type1, K> &m1, const coeff_type2 m2) {
    return m1 + scalar(m2);
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator - (const coeff_type1 m1, const multivector<coeff_type2, K> &m2) {
    return scalar(m1) + (-1 * m2);
}


template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator - (const multivector<coeff_type1, K> &m1, const coeff_type2 m2) {
    return m1 + scalar(-1 * m2);
}


template<class coeff_type, class U, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type, U>::type, K> operator * (const coeff_type &s, const multivector<U, K> &m) {
    multivector<typename std::common_type<coeff_type, U>::type, K> r;
    for (auto it = m.M.begin(); it != m.M.end(); ++it) {
        r.M[(*it).first] = s * (*it).second;
    }
    return r;
}

template<class coeff_type, class U, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type, U>::type, K> operator * (const multivector<U, K> &m, const coeff_type &s) {
    return s * m;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> generic_call (const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2, auto f ) {
    multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> result;
    for (auto it1 = m1.M.begin(); it1 != m1.M.end(); ++it1) {
        for (auto it2 = m2.M.begin(); it2 != m2.M.end(); ++it2) {
            result = result + f(it1->second, it1->first, it2->second, it2->first);
        }
    }
    return result;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator ^ (const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2) {
    auto OP = [&] (const coeff_type1 coef1, uint mask1, const coeff_type2 coef2, uint mask2) {
        multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> multivector_r;
        if (!(mask1 & mask2)) {
            multivector_r.M[(mask1 | mask2)] = canonical_sort(mask1, mask2) * coef1 * coef2;
            return multivector_r;
        }
    };
    return generic_call(m1, m2, OP);
}

template<int K = MAX_DIMENSIONS>
int canonical_sort(uint mask_1, uint mask_2) {
    int swaps = 0;
    mask_1 = mask_1 >> 1;
    while (hamming_weight(mask_1) != 0) {
        swaps += hamming_weight(mask_1 & mask_2);
        mask_1 = mask_1 >> 1;
    }
    if (hamming_weight(swaps & 1) == 0) {
        return 1;
    }
    return - 1;
}


/**
 *
 * CANDIDATO A REMOÇÃO
 *
 */
template<class coeff_type, int K = MAX_DIMENSIONS>
int check_grade(const multivector<coeff_type, K> &m) {
    auto it = m.M.begin();
    int count = hamming_weight((*it).first);
    ++it;
    while (it != m.M.end()) {
        if (hamming_weight((*it).first) != count) {
            return -1;
        }
        ++it;
    }
    return count;
}

int check_grade(int mask) {
    return hamming_weight(mask);
}

template<class coeff_type, int K = MAX_DIMENSIONS>
multivector<coeff_type, K> take_grade(const multivector<coeff_type, K> &m, const int grade) {
    multivector<coeff_type, K> multivector_r;
    for (auto it = m.M.begin(); it != m.M.end(); ++it) {
        if (check_grade((*it).first) == grade) {
            multivector_r.M[(*it).first] = (*it).second;
        }
    }
    return multivector_r;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> RP(const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2, int N_DIMS) {
    auto RP = [&] (const coeff_type1 coef1, uint mask1, const coeff_type2 coef2, uint mask2) {
        multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> multivector_r;
        uint mask_r = mask1 & mask2;
        if (check_grade(mask1) + check_grade(mask2) - check_grade(mask_r) == N_DIMS) {
            multivector_r.M[mask_r] = canonical_sort(mask1 ^ mask_r, mask2 ^ mask_r) * coef1 * coef2;
            return multivector_r;
        }
    };
    return generic_call(m1, m2, RP);
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> GP (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    auto GP = [&] (const T coef1, uint mask1, const U coef2, uint mask2) {
        multivector<typename std::common_type<T, U>::type, K> multivector_r;
        uint mask_r = mask1 ^ mask2;
        multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1, mask2)) * coef1 * coef2;
        return multivector_r;
    };
    return generic_call(m1, m2, GP);
}

template<class T, class U, int K>
typename std::common_type<T, U>::type  SCP (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    auto SCP = [&] (const T coef1, uint mask1, const U coef2, uint mask2) {
        multivector<typename std::common_type<T, U>::type, K> multivector_r;
        uint mask_r = mask1 ^ mask2;
        if (check_grade(mask_r) == 0) {
            multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1, mask2)) * coef1 * coef2;
        }
        return multivector_r;
    };
    return generic_call(m1, m2, SCP).M[0];
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> LCONT (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    auto LCONT = [&] (const T coef1, uint mask1, const U coef2, uint mask2) {
        multivector<typename std::common_type<T, U>::type, K> multivector_r;
        uint mask_r = mask1 ^ mask2;
        if (check_grade(mask_r) == (check_grade(mask2) - check_grade(mask1))) {
            multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1, mask2)) * coef1 * coef2;
        }
        return multivector_r;
    };
    return generic_call(m1, m2, LCONT);
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> RCONT (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    auto RCONT = [&] (const T coef1, uint mask1, const U coef2, uint mask2) {
        multivector<typename std::common_type<T, U>::type, K> multivector_r;
        uint mask_r = mask1 ^ mask2;
        if (check_grade(mask_r) == (check_grade(mask1) - check_grade(mask2))) {
            multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1, mask2)) * coef1 * coef2;
        }
        return multivector_r;
    };
    return generic_call(m1, m2, RCONT);
}

template<class coeff_type, int K = MAX_DIMENSIONS>
multivector<coeff_type, K> REVERSE(const multivector<coeff_type, K> &m) {
    // assert that m is a blade
    multivector<coeff_type, K> multivector_r;
    auto it = m.M.begin();
    int g = check_grade(m);
    int ex = ((int)(g * (g-1))) >> 1;
    multivector_r.M[it->first] = std::pow(-1, ex) * m.M.at(it->first);
    return multivector_r;
}

template<class coeff_type, int K = MAX_DIMENSIONS>
coeff_type SQR_NORM_REVERSE(const multivector<coeff_type, K> &m, Metric<coeff_type> &metric) {
    // assert that m is a blade
    return SCP(m, REVERSE(m), metric);
}

template<class coeff_type, int K = MAX_DIMENSIONS>
multivector<coeff_type, K> INV(const multivector<coeff_type, K> &m, Metric<coeff_type> &metric) {
    // assert that m is a blade
    multivector<coeff_type, K> multivector_r;
    multivector<coeff_type, K> temp = REVERSE(m);
    coeff_type norm = SQR_NORM_REVERSE(m, metric);
    auto it = temp.M.begin();
    multivector_r.M[it->first] = it->second / norm;
    return multivector_r;
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> IGP (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    multivector<typename std::common_type<T, U>::type, K> multivector_r = GP(m1, INV(m2, metric), metric);
    return multivector_r;
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K>   DELTA (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    multivector<typename std::common_type<T, U>::type, K> temp = GP(m1, m2, metric);
    std::cout << check_grade(m1) << std::endl;
    std::cout << check_grade(m2) << std::endl;
    auto it = temp.M.end();
    --it;
    int max = check_grade(it->first);
    std::cout << max << std::endl;
    return take_grade(temp, max);
}


#endif // MULTIVECTOR_H
