#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>
#include <type_traits>

#define MAX_DIMENSIONS 6

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

        multivector() {
        }


    template<class T, int K> friend std::ostream& operator << (std::ostream &, const multivector<T, K> &);
    template<class T, int K> friend multivector<T, K> e(int);
    template<class T, int K> friend multivector<T, K> scalar(T);

    template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator + (const multivector<T, K> &, const multivector<U, K> &);
    template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator - (const multivector<T, K> &, const multivector<U, K> &);
    template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator * (const T &, const multivector<U, K> &);
    template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> operator ^ (const multivector<T, K> &, const multivector<U, K> &);

    template<class T, class U, int K> friend multivector<typename std::common_type<T, U>::type, K> canonical_form (const multivector<T, K> &, const multivector<U, K> &);
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
int canonical_form (multivector<coeff_type1, K> &m1, multivector<coeff_type2, K> &m2) {
    return m2;
}

template<class coeff_type1, class coeff_type2, int K = MAX_DIMENSIONS>
multivector<typename std::common_type<coeff_type1, coeff_type2>::type, K> operator ^ (const multivector<coeff_type1, K> &m1, const multivector<coeff_type2, K> &m2) {
    return m2;
}

#endif // MULTIVECTOR_H
