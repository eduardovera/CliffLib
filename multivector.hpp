#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <functional>
#include <stdexcept>
#include <vector>
#include <assert.h>
#include <string>
#define MAX_BITS 25
typedef std::bitset<MAX_BITS> mask;
#include <metric.h>

#define TOLERANCE 0.00001



bool operator < (const std::bitset<MAX_BITS> &a, const std::bitset<MAX_BITS> &b) {
    return a.to_string().compare(b.to_string()) < 0;
}

bool operator >(const std::bitset<MAX_BITS> &a, const std::bitset<MAX_BITS> &b) {
    return a.to_string().compare(b.to_string()) > 0;
}


namespace CliffLib {

    template<class coeff_type>
    class multivector;
    static void build_masks();

    static int N_DIMS = -1;
    static std::map<int, multivector<double>> MASKS;

    enum multivector_type {
        BLADE = 0,
        VERSOR = 1,
        NO_MEANING = 2
    };


    template<std::size_t N>
    struct Comparer {
        bool operator () (const std::bitset<N> &x, const std::bitset<N> &y) const {
            return x.to_string().compare(y.to_string()) < 0;

//            std::bitset<MAX_BITS> a_c = x;
//            std::bitset<MAX_BITS> b_c = y;

//            for (int i = 0; i <= MAX_BITS; i++) {
//                if (a_c.test(i)) {
//                    a_c.flip(i);
//                }
//                if (b_c.test(i)) {
//                    b_c.flip(i);
//                }
//                if (a_c.count() == 0 && b_c.count() != 0) {
//                    return true;
//                }
//                if (b_c.count() == 0) {
//                    return false;
//                }
//            }

//            size_t y_first = y._Find_first();
//            return (x._Find_first() < y_first) && (y_first != y.size());
//            for (int i = N-1; i >= 0; i--) {
//                if (x[i] ^ y[i]) return y[i];
//            }
//            std::cout << x._Find_first() << std::endl;
//            std::cout << y._Find_first() << std::endl;

//            std::cout << x << std::endl;
//            std::cout << y << std::endl;

//            return x.to_ullong() < y.to_ullong();
        }
    };

    static std::map<mask, std::map<mask, multivector<double>, Comparer<MAX_BITS>>, Comparer<MAX_BITS>> LOOKUP_TABLE;


    template<class coeff_type>
    class multivector {

        private:

        public:

            std::map<mask, coeff_type, Comparer<MAX_BITS>> M;

            coeff_type getScalar() {
                auto it = M.find(mask(0));
                if (it == M.end()) {
                    return 0;
                }
                return it->second;
            }

            multivector() {

            }

            multivector<coeff_type>(const multivector<coeff_type> &m) {
                this->M = m.M;
            }

            std::string get_base(std::bitset<MAX_BITS> m) const {
                std::string ret = "";
                mask binary = m;
                for (int i = 1; i < MAX_BITS; i++) {
                    if (binary.test(i)) {
                        ret.append("(e" + std::to_string(i) + ")");
                        binary.flip(i);
                        if (binary.any()) {
                            ret.append("^");
                        }
                    }
                }
                return ret;
            }

            multivector<coeff_type> get_portion(long);
            template<class T> friend multivector<T> e(int);
            template<class T> friend multivector<T> scalar(T);
            template<class T> friend multivector<T> take_grade(const multivector<T> &, const int);
            template<class T> friend int get_grade(const multivector<T> &);
            template<class T> friend multivector<T> REVERSE(const multivector<T> &);
            template<class T> friend multivector<T> INV(const multivector<T> &, Metric<double> &);
            template<class T> friend T SQR_NORM_REVERSE(const multivector<T> &, Metric<double> &);

            template<class T> friend std::ostream& operator << (std::ostream &, const multivector<T> &);
            template<class T, class U> friend bool operator == (const multivector<T> &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator + (const multivector<T> &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator - (const multivector<T> &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator * (const T &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> operator ^ (const multivector<T> &, const multivector<U> &);

            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> generic_call (const multivector<T> &m1, const multivector<U> &m2, auto f );

            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> RP (const multivector<T> &, const multivector<U> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> GP (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> IGP (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> LCONT (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> RCONT (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend multivector<typename std::common_type<T, U>::type> DELTA (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);
            template<class T, class U> friend typename std::common_type<T, U>::type  SCP (const multivector<T> &, const multivector<U> &, Metric<typename std::common_type<T, U>::type> &);

            template<class T> friend std::vector<multivector<T>> FACTORIZE(const multivector<T> &, Metric<double> &, double &);

            template<class T> friend multivector<T> DUAL(const multivector<T> &, Metric<double> &);
            template<class T, class U> friend multivector<T> DUAL(const multivector<T> &, const multivector<U> &, Metric<Metric<typename std::common_type<T, U>::type>> &);

            template<class T> friend multivector<T> UNDUAL(const multivector<T> &, Metric<double> &);
            template<class T, class U> friend multivector<T> UNDUAL(const multivector<T> &, const multivector<U> &, Metric<Metric<typename std::common_type<T, U>::type>> &);

            template<class T, class U> friend std::vector<multivector<typename std::common_type<T, U>::type>> MEET_AND_JOIN(const multivector<T> &, const multivector<U> &, int, Metric<typename std::common_type<T, U>::type> &);

            template<class T> friend multivector<T> PSEUDOSCALAR();
            template<class T> friend multivector<T> SCALAR();

            template<class T> friend multivector<T> GRADE_INVOL(const multivector<T> &m);

            bool isZero() const {
                for (auto it = M.begin(); it != M.end(); ++it) {
                    if (it->second != 0) {
                        return false;
                    }
                }
                return true;
            }

            void handle_numeric_error() {
                for (auto it = this->M.begin(); it != this->M.end(); ++it) {
                    if (std::abs(it->second) < TOLERANCE) {
                        this->M.erase(it);
                    }
                }
            }

            mask getMask() {
                return M.begin()->first;
            }

            multivector_type get_type(Metric<coeff_type> &metric) const {
                if (SQR_NORM_REVERSE(*this, metric) < TOLERANCE) {
                    return NO_MEANING;
                }
                auto grade_preserving = [] (const multivector<coeff_type> &involution, const multivector<coeff_type> &reverse, Metric<coeff_type> &metric) {
                    for (int i = 1; i <= N_DIMS; i++) {
                        mask masc(1 << i);
                        multivector<coeff_type> e_;
                        e_.M[masc] = 1;
                        std::cout << "e_" << i << ": " << e_ << std::endl;
                        multivector<coeff_type> M = GP(GP(involution, e_, metric), reverse, metric);
                        M.handle_numeric_error();
                        std::cout << "M: " << M << std::endl;
                        if (get_grade(M) != 1) {
                            return false;
                        }
                    }
                    return true;
                };

                multivector<coeff_type> involution = GRADE_INVOL(*this);
                multivector<coeff_type> inverse = INV(*this, metric);

                multivector<coeff_type> m = GP(involution, inverse, metric);
                m.handle_numeric_error();
                multivector<coeff_type> diff = m - (GP(inverse, involution, metric));
                diff.handle_numeric_error();

                std::cout << "m: " << m << std::endl;
                std::cout << "diff: " << diff << std::endl;
                std::cout << "grade(m): " << get_grade(m) << std::endl;
                std::cout << "diff is zero? " << diff.isZero() << std::endl;


                if (get_grade(m) == 0 && diff.isZero()) {
                    if (grade_preserving(involution, REVERSE(*this), metric)) {
                        if (get_grade(*this) != -1) {
                            return BLADE;
                        } else {
                            return VERSOR;
                        }
                    } else {
                        return NO_MEANING;
                    }
                } else {
                    return NO_MEANING;
                }
            }


    };

    template<class T>
    int get_grade(const multivector<T> &m) {
        if (m.isZero()) {
            return 0;
        }
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
    multivector<T> PSEUDOSCALAR() {
        mask masc;
        for (int i = 1; i <= N_DIMS; i++) {
            masc[i] = 1;
        }
        multivector<T> I;
        I.M[masc] = 1;
        return I;
    }

    template<class T>
    multivector<T> SCALAR() {
        mask masc(0);
        multivector<T> s;
        s.M[masc] = 1;
        return s;
    }

    template<class T>
    std::ostream& operator << (std::ostream &os, const multivector<T> &m) {
        for (auto it = m.M.begin(); it != m.M.end(); ++it) {
            os << (it->second > 0.0 ? "+" : "") << it->second << m.get_base(it->first);
        }
        return os;
    }

    template<class T, class U>
    bool operator == (const multivector<T> &m1, const multivector<U> &m2) {
        auto it1 = m1.M.begin();
        auto it2 = m2.M.begin();
        while (it1 != m1.M.end() || it2 != m2.M.end()) {
            if (it1->first != it2->first || it1->second != it2->second) {
                return false;
            }
            ++it1;
            ++it2;
        }
        if (it1 == m1.M.end() && it2 != m2.M.end()) {
            return false;
        } else if (it2 == m2.M.end() && it1 != m1.M.end()) {
            return false;
        }
        return true;
    }


    template<class T = double>
    multivector<T> scalar(T value) {
        multivector<T> r;
        r.M[mask(0)] = value;
        return r;
    }


    template<class T = double>
    multivector<T> e(int index) {
        if (N_DIMS == -1) {
            throw std::invalid_argument("Please define the dimensionality of your subspace by setting CliffLib::N_DIMS param.");
        }
        if (index >= MAX_BITS) {
            std::string message = "Can't represent e(" + std::to_string(index) + ") with this amount of bits. Please select a higher value for MAX_BITS param by using #define MAX_BITS directive before including multivector.hpp";
            throw std::invalid_argument(message);
        }
        multivector<T> r;
        mask k(1);
        r.M[k << index] = 1;
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
            if (it1->first < it2->first) {
                multivector_r.M[it1->first] = it1->second;
                it1++;
            } else if (it1->first > it2->first) {
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
        multivector<typename std::common_type<T, U>::type> partial;
        for (auto it1 = m1.M.begin(); it1 != m1.M.end(); ++it1) {
            for (auto it2 = m2.M.begin(); it2 != m2.M.end(); ++it2) {
                partial = f(it1->second, it1->first, it2->second, it2->first);
                result = result + partial;
            }
        }

        result.handle_numeric_error();
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
    multivector<typename std::common_type<T, U>::type> RP(const multivector<T> &m1, const multivector<U> &m2) {
        auto RP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_r = mask1 & mask2;
            if ((int)(mask1.count() + mask2.count() - mask_r.count()) == N_DIMS) {
                multivector_r.M[mask_r] = canonical_sort(mask1 ^ mask_r, mask2 ^ mask_r) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, RP);
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> GP (const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {

        auto GP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            mask mask_and = mask1 & mask2;
            mask mask_xor = mask1 ^ mask2;
            multivector<typename std::common_type<T, U>::type> multivector_r;
            multivector_r.M[mask_xor] = canonical_sort(mask1, mask2) * metric.factorByMask((mask_and), (mask_and)) * coef1 * coef2;
            return multivector_r;
        };


//        auto GP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
////            std::cout << mask1 << std::endl;
////            std::cout << mask2 << std::endl;
//            multivector<typename std::common_type<T, U>::type> multivector_r;
//            multivector_r.M[mask_xor] = coef1 * coef2;
//            if (mask1.any() && mask2.any()) {
//                return CliffLib::LOOKUP_TABLE[mask1][mask2] * coef1 * coef2;
//            } else {
//                mask mask_xor = mask1 ^ mask2;
//                if (mask1.any()) {
////                } else {
////                    multivector_r.M[mask2] = coef1 * coef2;
//                }
//                return multivector_r;
//            }
//        };
        return generic_call(m1, m2, GP);
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> GP (T s, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        return GP(scalar(s), m2, metric);
    }


    template<class T, class U>
    typename std::common_type<T, U>::type  SCP (const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        auto SCP = [&] (const T coef1, mask mask1, const U coef2, mask mask2) {
            multivector<typename std::common_type<T, U>::type> multivector_r;
            mask mask_and = mask1 & mask2;

            mask mask_xor = mask1 ^ mask2;
//            std::cout << mask_xor << std::endl;
            if (mask_xor.count() == 0) {
                multivector_r.M[mask_xor] = canonical_sort(mask1, mask2) * metric.factorByMask(mask1, mask2) * coef1 * coef2;
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
                multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1 & mask2).to_ullong(), (mask1 & mask2).to_ullong()) * coef1 * coef2;
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
                multivector_r.M[mask_r] = canonical_sort(mask1, mask2) * metric.factor((mask1 & mask2).to_ullong(), (mask1 & mask2).to_ullong()) * coef1 * coef2;
            }
            return multivector_r;
        };
        return generic_call(m1, m2, RCONT);
    }

    template<class T>
    multivector<T> REVERSE(const multivector<T> &m) {
        multivector<T> multivector_r;
        for (auto it = m.M.begin(); it != m.M.end(); ++it) {
            int g = it->first.count();
            int ex = (g * (g-1)) >> 1 ;
            multivector_r.M[it->first] = std::pow(-1, ex) * it->second;
        }
        return multivector_r;
    }

    template<class T>
    T SQR_NORM_REVERSE(const multivector<T> &m, Metric<double> &metric) {
        return SCP(m, REVERSE(m), metric);
    }

    template<class T>
    multivector<T> INV(const multivector<T> &m, Metric<double> &metric) {
        double norm = SQR_NORM_REVERSE(m, metric);
        if (std::abs(norm) < TOLERANCE) {
            throw std::invalid_argument("Multivector not invertible");
        }
        double inv_norm = 1.0 / norm;
        multivector<T> r = (REVERSE(m) * inv_norm);
        r.handle_numeric_error();
        return r;
    }

    template<class T, class U>
    multivector<typename std::common_type<T, U>::type> IGP(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        multivector<typename std::common_type<T, U>::type> multivector_r = GP(m1, INV(m2, metric), metric);
        multivector_r.handle_numeric_error();
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
    multivector<T> ORTHO_PROJECT(const multivector<T> &m1, const multivector<U> &m2, Metric<double> &metric) {
        return LCONT(LCONT(m1, INV(m2, metric), metric), m2, metric);
    }

    template<class T, class U>
    multivector<T> ORTHO_REJECT(const multivector<T> &m1, const multivector<U> &m2, Metric<double> &metric) {
        return RCONT((m1^m2), INV(m2, metric), metric);
    }

    template<class T>
    multivector<T> GRADE_INVOL(const multivector<T> &m) {
        multivector<T> invol = m;
        for (auto it = invol.M.begin(); it != invol.M.end(); ++it) {
            invol.M[it->first] = pow(-1, it->first.count()) * it->second;
        }
        return invol;
    }

    template<class T>
    std::vector<multivector<T>> FACTORIZE(const multivector<T> &m, Metric<double> &metric, double &scaling) {
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

        for (int i = 1; i < MAX_BITS; i++) {
            if (masc.test(i)) {
                mask e_ = masc.none();
                e_[i] = 1;
                multivector<T> e_j;
                e_j.M[1 << i] = 1;

                multivector<T> proj = ORTHO_PROJECT(e_j, temp, metric);
                multivector<T> factor_j = proj * (1.0 / sqrt(SQR_NORM_REVERSE(proj, metric)));
                r.push_back(factor_j);
                temp = LCONT(INV(factor_j, metric), temp, metric);
            }
        }

        return r;
    }

    template<class T>
    multivector<T> DUAL(const multivector<T> &m, Metric<double> &metric) {
        return LCONT(m, INV(PSEUDOSCALAR<T>(), metric), metric);
    }

    template<class T, class U>
    multivector<T> DUAL(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        return LCONT(m1, INV(m2, metric), metric);
    }

    template<class T>
    multivector<T> UNDUAL(const multivector<T> &m, Metric<double> &metric) {
        return LCONT(m, PSEUDOSCALAR<T>(), metric);
    }

    template<class T, class U>
    multivector<T> UNDUAL(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        return LCONT(m1, m2, metric);
    }

    template<class T, class U>
    std::vector<multivector<typename std::common_type<T, U>::type>> MEET_AND_JOIN(const multivector<T> &m1, const multivector<U> &m2, Metric<typename std::common_type<T, U>::type> &metric) {
        assert(m1.get_type(metric) == BLADE);
        assert(m2.get_type(metric) == BLADE);
        std::vector<multivector<typename std::common_type<T, U>::type>> output(2);

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
        double scaling;
        std::vector<multivector<typename std::common_type<T, U>::type>> factors = FACTORIZE(DUAL(delta, metric), metric, scaling);

        multivector<typename std::common_type<T, U>::type> meet = SCALAR<typename std::common_type<T, U>::type>();
        multivector<typename std::common_type<T, U>::type> join = PSEUDOSCALAR<typename std::common_type<T, U>::type>();

        for (auto f_j : factors) {
            auto proj = ORTHO_PROJECT(f_j, A, metric);
            if (!proj.isZero()) {
                meet = meet ^ proj;
                if (get_grade(meet) == t) {
                    join = RCONT(A, INV(meet, metric), metric)^B;
                    continue;
                }
            }
            auto rej = ORTHO_REJECT(f_j, A, metric);
            if (!rej.isZero()) {
                join = LCONT(rej, join, metric);
                if (get_grade(join) == (r+s-t)) {
                    meet = UNDUAL(DUAL(B, metric) ^ DUAL(A, metric), metric);
                    continue;
                }
            }
        }

        if (r > s) {
            output[0] = meet * pow(-1, (r-t)*(s-t));
            output[1] = join * pow(-1, (r-t)*(s-t));
        } else {
            output[0] = meet;
            output[1] = join;
        }
        return output;
    }

    static void build_lookup_table(OrthonormalMetric<double> m) {

        if (N_DIMS == -1) {
            throw std::invalid_argument("Please define the dimensionality of your subspace by setting CliffLib::N_DIMS param.");
        }
        auto GP_MOD = [&] (const double coef1, mask mask1, const double coef2, mask mask2) {
            mask mask_and = mask1 & mask2;
            mask mask_xor = mask1 ^ mask2;
            multivector<double> multivector_r;
            multivector_r.M[mask_xor] = canonical_sort(mask1, mask2) * m.factorByMask((mask_and), (mask_and));
//            std::cout << multivector_r << std::endl;
            return multivector_r;
        };

//        #pragma omp parallel for
        for (int j = 0; j < MAX_BITS; j++) {
            for (int i = 0; i < MAX_BITS; i++) {
                LOOKUP_TABLE[e(j).getMask()][e(i).getMask()] = generic_call(e(j), e(i), GP_MOD);
            }
        }
    }

//#pragma omp parallel num_threads(30)
    static void build_masks() {
        if (N_DIMS == -1) {
            throw std::invalid_argument("Please define the dimensionality of your subspace by setting CliffLib::N_DIMS param.");
        }
//        #pragma omp parallel for
        for (int i = 0; i <= N_DIMS; i++) {
            MASKS[i] = e(i);
        }
    }

    template<class T>
    multivector<T> multivector<T>::get_portion(long index) {
        auto it = M.find(mask(1) << index);
        if (it == M.end()) {
            return 0 * e<T>(index);
        }
        return e<T>(index);
    }
}

#endif // MULTIVECTOR_H
