#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <utils.h>
#include <metric.h>

#define MAX_DIMENSIONS 10
#ifdef __MINGW32__
typedef uint8_t uint;
#endif

template<int K = MAX_DIMENSIONS>
int canonical_sort(const int, const int);


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
        template<class T, int K> friend int take_grade(const multivector<T, K> &);
        template<class T, int K> friend multivector<T, K> get_element_with_grade(const multivector<T, K> &, const int);

        template<class T, int K> friend std::ostream& operator << (std::ostream &, const multivector<T, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator + (const multivector<T, K> &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator - (const multivector<T, K> &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator * (const T &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator ^ (const multivector<T, K> &, const multivector<U, K> &);

        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> RP (const multivector<T, K> &, const multivector<U, K> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> GP (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> SCP (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> LCONT (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);
        template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> RCONT (const multivector<T, K> &, const multivector<U, K> &, Metric<typename std::common_type<T, U>::type> &);


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
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator ^ (const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2) {
    multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> multivector_r;

    for (auto it1 = m1.M.begin(); it1 != m1.M.end(); ++it1) {
        for (auto it2 = m2.M.begin(); it2 != m2.M.end(); ++it2) {
            if (!((*it1).first & (*it2).first)) {
                multivector_r.M[((*it1).first | (*it2).first)] = canonical_sort((*it1).first, (*it2).first) * (*it1).second * (*it2).second;
            }
        }
    }
    return multivector_r;
}

template<int K = MAX_DIMENSIONS>
int canonical_sort(int mask_1, int mask_2) {
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
int take_grade(const multivector<coeff_type, K> &m) {
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

int take_grade(int mask) {
    return hamming_weight(mask);
}

template<class coeff_type, int K = MAX_DIMENSIONS>
multivector<coeff_type, K> get_element_with_grade(const multivector<coeff_type, K> &m, const int grade) {
    multivector<coeff_type, K> multivector_r;
    for (auto it = m.M.begin(); it != m.M.end(); ++it) {
        if (take_grade((*it).first) == grade) {
            multivector_r.M[(*it).first] = (*it).second;
        }
    }
    return multivector_r;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> RP(const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2) {
    multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> multivector_r;
    int n = std::max(take_grade(m1), take_grade(m2));
    for (auto it1 = m1.M.begin(); it1 != m1.M.end(); ++it1) {
        for (auto it2 = m2.M.begin(); it2 != m2.M.end(); ++it2) {
            int mask_r = (*it1).first & (*it2).first;
            if ((take_grade((*it1).first) + take_grade((*it2).first) - take_grade(mask_r)) == n) {
                multivector_r.M[mask_r] += canonical_sort((*it1).first ^ mask_r, (*it2).first ^ mask_r) * (*it1).second * (*it2).second;
            }
        }
    }
    return multivector_r;
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> GP (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    multivector<typename std::common_type<T, U>::type, K> multivector_r;
    for (auto it1 = m1.M.begin(); it1 != m1.M.end(); ++it1) {
        for (auto it2 = m2.M.begin(); it2 != m2.M.end(); ++it2) {
            int mask_r = (*it1).first ^ (*it2).first;
            multivector_r.M[mask_r] += canonical_sort((*it1).first, (*it2).first) * metric.factor((*it1).first & (*it2).first) * (*it1).second * (*it2).second;
        }
    }
    return multivector_r;
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> SCP (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    multivector<typename std::common_type<T, U>::type, K> multivector_r;
    multivector<typename std::common_type<T, U>::type, K> temp = GP(m1, m2, metric);
    multivector_r.M[0] = temp.M[0];
    return multivector_r;
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> LCONT (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    multivector<typename std::common_type<T, U>::type, K> multivector_r;
    multivector<typename std::common_type<T, U>::type, K> temp = GP(m1, m2, metric);
    if (take_grade(m1) <= take_grade(m2)) {
        multivector_r = get_element_with_grade(temp, take_grade(m2) - take_grade(m1));
    }
    return multivector_r;
}

template<class T, class U, int K>
multivector<typename std::common_type<T, U>::type, K> RCONT (const multivector<T, K> &m1, const multivector<U, K> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
    multivector<typename std::common_type<T, U>::type, K> multivector_r;
    multivector<typename std::common_type<T, U>::type, K> temp = GP(m1, m2, metric);
    if (take_grade(m1) >= take_grade(m2)) {
        multivector_r = get_element_with_grade(temp, take_grade(m1) - take_grade(m2));
    }
    return multivector_r;
}


#endif // MULTIVECTOR_H
